//
//  Curvature.cpp
//  OpenFlipper
//
//  Created by 饒亢 on 2015/1/7.
//
//

#include "PlushPatternGenerator.hh"

#include <QTextStream>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Monge_via_jet_fitting.h>
#include <CGAL/Eigen_svd.h>

/** Expand selection by n-ring connectivity **/
void expandVerticeSelection(TriMesh *mesh, std::set<VertexHandle> &verticesSelection, int n) {
    for (int i = 0; i < n; i++) {
        std::vector<VertexHandle> ring;
        for (std::set<VertexHandle>::iterator it = verticesSelection.begin(); it != verticesSelection.end(); it++) {
            for (TriMesh::VertexVertexIter vv_it = mesh->vv_iter(*it); vv_it++;) {
                ring.push_back(*vv_it);
            }
        }
        verticesSelection.insert(ring.begin(), ring.end());
    }
}

class Vertex_distance
{
public:
    VertexHandle vertex;
    double distance;

    Vertex_distance(VertexHandle v, double d) : vertex(v), distance(d) {}
};

bool operator< (const Vertex_distance &A, const Vertex_distance &B) {
    return A.distance < B.distance;
}

bool PlushPatternGenerator::calcCurvature() {
    typedef double                   DFT;
    typedef CGAL::Simple_cartesian<DFT>     Data_Kernel;
    typedef CGAL::Simple_cartesian<DFT>     Local_Kernel;
    typedef CGAL::Eigen_svd SvdTraits;
    typedef Data_Kernel::Point_3     DPoint;
    typedef CGAL::Monge_via_jet_fitting<Data_Kernel,Local_Kernel,SvdTraits> My_Monge_via_jet_fitting;
    typedef My_Monge_via_jet_fitting::Monge_form     My_Monge_form;

    size_t d_fitting = 2;
    size_t d_monge = 2;
    
    isJobCanceled = false;
    
    std::vector<int> selectedVertices;
    // If no one is selected, show all of them
    if (selectedVertices.size() == 0) {
        for (VertexIter v_it = m_mesh->vertices_begin(); v_it != m_mesh->vertices_end(); v_it++) {
            selectedVertices.push_back(v_it->idx());
        }
    }
    
    // for each selected vertex, calculate curvature
    for (size_t i = 0; i < selectedVertices.size(); i++) {
        if (isJobCanceled) {
            emit log(LOGINFO, "Curvature calculation canceled.");
            return false;
        }
        
        VertexHandle centerV = m_mesh->vertex_handle(selectedVertices[i]);

        std::set<VertexHandle> neighboringVertices;
        neighboringVertices.insert(centerV);

        // Expand selection until we got enough samples
        // Or terminate if selection doesn't grow anymore
        size_t requiredSamples = (d_fitting+1)*(d_fitting+2)/2;
        size_t nNeighbors;
        do {
            nNeighbors = neighboringVertices.size();
            expandVerticeSelection(m_mesh, neighboringVertices, 1);
        } while (nNeighbors != neighboringVertices.size() && neighboringVertices.size() < requiredSamples*3);
        
        if (nNeighbors == neighboringVertices.size() && nNeighbors < requiredSamples) {
            QString str = QString("Not enough points for fitting for mesh: %1").arg(m_meshName);
            emit log(LOGERR, str);
            continue;
        }
        
        std::vector<DPoint> in_points;
        for (std::set<VertexHandle>::iterator it = neighboringVertices.begin(); it != neighboringVertices.end(); it++)
        {
            TriMesh::Point p = m_mesh->point(*it);
            in_points.push_back(DPoint(p[0],p[1],p[2]));
        }
        
        My_Monge_form monge_form;
        My_Monge_via_jet_fitting monge_fit;
        monge_form = monge_fit(in_points.begin(), in_points.end(), d_fitting, d_monge);
        // fix with normal
        TriMesh::Normal normal = m_mesh->normal(centerV);
        CGAL::Vector_3<Data_Kernel> normalACGL(normal[0], normal[1], normal[2]);
        monge_form.comply_wrt_given_normal(normalACGL);
        
        // Save result
        OpenMesh::Vec3d d1(monge_form.maximal_principal_direction()[0],
                           monge_form.maximal_principal_direction()[1],
                           monge_form.maximal_principal_direction()[2]);
        OpenMesh::Vec3d d2(monge_form.minimal_principal_direction()[0],
                           monge_form.minimal_principal_direction()[1],
                           monge_form.minimal_principal_direction()[2]);
        m_mesh->property(maxCurvatureDirectionHandle, centerV) = d1;
        m_mesh->property(minCurvatureDirectionHandle, centerV) = d2;
        m_mesh->property(maxCurvatureHandle, centerV) = monge_form.principal_curvatures(0);
        m_mesh->property(minCurvatureHandle, centerV) = monge_form.principal_curvatures(1);
        
        int status = (int)((double)i / selectedVertices.size() * 100);
        emit setJobState(status);
    }
    
    // Prepare file for saving data
    QFile file(m_meshName + "_curvature.txt");
    file.open(QIODevice::WriteOnly | QIODevice::Text);
    QTextStream out(&file);

    for (VertexIter v_it = m_mesh->vertices_begin(); v_it != m_mesh->vertices_end(); v_it++)
    {
        OpenMesh::Vec3d d1 = m_mesh->property(maxCurvatureDirectionHandle, *v_it);
        OpenMesh::Vec3d d2 = m_mesh->property(minCurvatureDirectionHandle, *v_it);
        out << d1[0] << " " << d1[1] << " " << d1[2] << " ";
        out << d2[0] << " " << d2[1] << " " << d2[2] << " ";
        out << m_mesh->property(maxCurvatureHandle, *v_it) << " ";
        out << m_mesh->property(minCurvatureHandle, *v_it) << "\n";
    }
    file.close();
    
    return true;
}


void PlushPatternGenerator::loadCurvature() {
    QString curvatureFilename = m_meshName + "_curvature.txt";
    QFile file(curvatureFilename);
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
        emit log(LOGERR, QString("Unable to read file %1").arg(curvatureFilename));
        return;
    }
    
    QTextStream in(&file);
    
    double minK1 = 1e9, minK2 = 1e9;
    double d1x, d1y, d1z, d2x, d2y, d2z, k1, k2;
    for (VertexIter v_it = m_mesh->vertices_begin(); v_it != m_mesh->vertices_end(); v_it++) {
        in >> d1x >> d1y >> d1z >> d2x >> d2y >> d2z >> k1 >> k2;
        m_mesh->property(maxCurvatureDirectionHandle, *v_it) = OpenMesh::Vec3d(d1x, d1y, d1z);
        m_mesh->property(minCurvatureDirectionHandle, *v_it) = OpenMesh::Vec3d(d2x, d2y, d2z);
        m_mesh->property(maxCurvatureHandle, *v_it) = k1;
        m_mesh->property(minCurvatureHandle, *v_it) = k2;
    }
}
