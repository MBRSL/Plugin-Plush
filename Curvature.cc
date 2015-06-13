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
    
    // Gaussian curvature
    for (VertexHandle v : m_mesh->vertices()) {
        m_mesh->property(gaussianCurvatureHandle, v) =  curvature::gauss_curvature(*m_mesh, v);
    }

    // Mean
    for (VertexHandle v : m_mesh->vertices()) {
        ACG::Vec3d curva(0.0,0.0,0.0);
        double area = 0.0;
        curvature::discrete_mean_curv_op<TriMesh,ACG::Vec3d,double>(*m_mesh, v, curva, area);
        double curv  = curva.norm();
        
        if ((curva | m_mesh->normal(v)) < 0)
            curv = -curv;
        
        m_mesh->property(meanCurvatureHandle, v) =  curv;
    }

    // Prepare file for saving data
    QFile mean_file(m_meshName + "_mean_curvature.txt");
    mean_file.open(QIODevice::WriteOnly | QIODevice::Text);
    QTextStream mean_out(&mean_file);

    for (VertexHandle v : m_mesh->vertices()) {
        mean_out << m_mesh->property(meanCurvatureHandle, v) << " ";
    }
    mean_file.close();

    QFile gaussian_file(m_meshName + "_gaussian_curvature.txt");
    gaussian_file.open(QIODevice::WriteOnly | QIODevice::Text);
    QTextStream gaussian_out(&gaussian_file);
    
    for (VertexHandle v : m_mesh->vertices()) {
        gaussian_out << m_mesh->property(meanCurvatureHandle, v) << " ";
    }
    gaussian_file.close();

    return true;
}


void PlushPatternGenerator::loadCurvature() {
    QFile mean_file(m_meshName + "_mean_curvature.txt");
    if (!mean_file.open(QIODevice::ReadOnly | QIODevice::Text)) {
        emit log(LOGERR, QString("Unable to read file %1").arg(m_meshName + "_mean_curvature.txt"));
        return;
    }
    QTextStream mean_in(&mean_file);
    double mean_curvature;
    for (VertexHandle v : m_mesh->vertices()) {
        mean_in >> mean_curvature;
        m_mesh->property(gaussianCurvatureHandle, v) = mean_curvature;
    }
    mean_file.close();
    
    QFile gaussian_file(m_meshName + "_gaussian_curvature.txt");
    if (!gaussian_file.open(QIODevice::ReadOnly | QIODevice::Text)) {
        emit log(LOGERR, QString("Unable to read file %1").arg(m_meshName + "_gaussian_curvature.txt"));
        return;
    }
    QTextStream gaussian_in(&gaussian_file);
    double gaussian_curvature;
    for (VertexHandle v : m_mesh->vertices()) {
        gaussian_in >> gaussian_curvature;
        m_mesh->property(gaussianCurvatureHandle, v) = gaussian_curvature;
    }
    gaussian_file.close();
}
