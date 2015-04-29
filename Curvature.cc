//
//  Curvature.cpp
//  OpenFlipper
//
//  Created by 饒亢 on 2015/1/7.
//
//

#include "PlushPatternGenerator.hh"

#include <MeshTools/Curvature.hh>

#include <QTextStream>

bool PlushPatternGenerator::calcCurvature() {
    isJobCanceled = false;
    
    // Gause curvature
//    for (VertexHandle v : m_mesh->vertices()) {
//        m_mesh->property(maxCurvatureHandle, v) =  curvature::gauss_curvature(*m_mesh, v);
//    }

    // Mean
    for (VertexHandle v : m_mesh->vertices()) {
        ACG::Vec3d curva(0.0,0.0,0.0);
        double area = 0.0;
        curvature::discrete_mean_curv_op<TriMesh,ACG::Vec3d,double>(*m_mesh, v, curva, area);
        double curv  = curva.norm();
        
        if ((curva | m_mesh->normal(v)) < 0)
            curv = -curv;
        
        m_mesh->property(maxCurvatureHandle, v) =  curv;
    }

    // Prepare file for saving data
    QFile file(m_meshName + "_curvature.txt");
    file.open(QIODevice::WriteOnly | QIODevice::Text);
    QTextStream out(&file);

    for (VertexHandle v : m_mesh->vertices()) {
        out << m_mesh->property(maxCurvatureHandle, v) << " ";
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
    
    double mean_curvature;
    for (VertexHandle v : m_mesh->vertices()) {
        in >> mean_curvature;
        m_mesh->property(maxCurvatureHandle, v) = mean_curvature;
    }
}
