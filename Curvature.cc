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
    
    double avgEdgeLength = 0;
    for (EdgeHandle eh : m_mesh->edges()) {
        avgEdgeLength += m_mesh->calc_edge_length(eh);
    }
    avgEdgeLength /= m_mesh->n_edges();

    // for each selected vertex, calculate curvature
    int count = 0;
    for (VertexHandle centerV : m_mesh->vertices()) {
        if (isJobCanceled) {
            emit log(LOGINFO, "Curvature calculation canceled.");
            return false;
        }
        
        std::set<VertexHandle> neighboringVertices;
        neighboringVertices.insert(centerV);

        // Expand selection until we got enough samples
        // Or terminate if selection doesn't grow anymore
        double maxDistance = avgEdgeLength * 2;
        size_t requiredSamples = (d_fitting+1)*(d_fitting+2)/2;
        while (neighboringVertices.size() < requiredSamples*3) {
            size_t nNeighbors = neighboringVertices.size();
            expandVertice(m_mesh, centerV, neighboringVertices, 1, maxDistance);
            
            // We can't expand selection, try large area
            if (nNeighbors == neighboringVertices.size()) {
                maxDistance *= 2;
                if (maxDistance < avgEdgeLength*100) {
                    continue;
                } else {
                    assert("Not enough sample points for fitting. And it seems there are no possible ways to expand selection.");
                }
            }
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
        m_mesh->property(maxCurvatureHandle, centerV) = monge_form.principal_curvatures(0);
        
        int status = (int)((double)++count / m_mesh->n_vertices() * 100);
        emit setJobState(status);
    }
    
    // Prepare file for saving data
    QFile file(m_meshName + "_curvature.txt");
    file.open(QIODevice::WriteOnly | QIODevice::Text);
    QTextStream out(&file);

    for (VertexIter v_it = m_mesh->vertices_begin(); v_it != m_mesh->vertices_end(); v_it++)
    {
        out << m_mesh->property(maxCurvatureHandle, *v_it) << " ";
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
    double k1;
    for (VertexIter v_it = m_mesh->vertices_begin(); v_it != m_mesh->vertices_end(); v_it++) {
        in >> k1;
        m_mesh->property(maxCurvatureHandle, *v_it) = k1;
    }
}
