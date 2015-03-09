//
//  Skeleton.cpp
//  OpenFlipper
//
//  Created by 饒亢 on 2015/1/31.
//
//

#include "PlushPatternGenerator.hh"
#include "SuperDeform/Skeleton.hh"
#include "SuperDeform/Weight.hh"

bool PlushPatternGenerator::calcSkeletonWeight() {
    Skeleton *skeleton = m_mesh->property(skeletonHandle);
    if (!skeleton) {
        loadSkeleton();
        if (!skeleton) {
            return false;
        }
        skeleton = m_mesh->property(skeletonHandle);
    }
    
    Weight weightGenerator;
    m_mesh->request_face_normals();
    weightGenerator.computeBoneWeight(m_mesh, skeleton);
    m_mesh->release_face_normals();
    
    Eigen::MatrixXd weight = weightGenerator.getWeight();
    for (VertexIter v_it = m_mesh->vertices_begin(); v_it != m_mesh->vertices_end(); v_it++) {
        double *w = new double[skeleton->bones.size()];
        for (size_t i = 0; i < skeleton->bones.size(); i++) {
            w[i] = weight(i, v_it->idx());
        }
        if (m_mesh->property(bonesWeightHandle, *v_it) != NULL) {
            delete[] m_mesh->property(bonesWeightHandle, *v_it);
        }
        m_mesh->property(bonesWeightHandle, *v_it) = w;
    }
    return true;
}

bool PlushPatternGenerator::loadBoneWeight() {
    Skeleton *skeleton = m_mesh->property(skeletonHandle);
    if(skeleton == NULL)
    {
        // No skeleton yet, try to load skeleton
        if (!loadSkeleton()) {
            return false;
        }
        skeleton = m_mesh->property(skeletonHandle);
    }
    
    QString weightFilename = m_meshName + ".skeleton.weight";
    QFile file(weightFilename);
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text))
    {
        emit log(LOGERR, QString("Error opening file: %1").arg(weightFilename));
        return false;
    }

    QTextStream fin(&file);
    // Assign bone weight to vertices
    for(VertexIter v_it = m_mesh->vertices_begin(); v_it != m_mesh->vertices_end(); v_it++)
    {
        double *weights = new double[skeleton->bones.size()];
        for(size_t j = 0; j < skeleton->bones.size(); j++)
        {
            double w;
            fin >> w;
            weights[j] = w;
        }
        if (m_mesh->property(bonesWeightHandle, *v_it) != NULL) {
            delete[] m_mesh->property(bonesWeightHandle, *v_it);
        }
        m_mesh->property(bonesWeightHandle, *v_it) = weights;
    }
    file.close();

    return true;
}

bool PlushPatternGenerator::saveBoneWeight() {
    QString weightFilename = m_meshName + ".skeleton.weight";
    QFile file(weightFilename);
    if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) {
        emit log(LOGERR, QString("Error opening file: %1").arg(weightFilename));
        return false;
    }
    
    QTextStream out(&file);
    int nBones = m_mesh->property(skeletonHandle)->bones.size();
    for (VertexIter v_it = m_mesh->vertices_begin(); v_it != m_mesh->vertices_end(); v_it++) {
        double *w = m_mesh->property(bonesWeightHandle, *v_it);
        for (int i = 0; i < nBones; i++) {
            out << w[i] << " ";
        }
        out << endl;
    }
    file.close();
    return true;
}

bool PlushPatternGenerator::loadSkeleton() {
    QString skeletonFilename = m_meshName + ".skeleton";
    
    // Add skeleton to mesh property
    Skeleton *skeleton = new Skeleton;
    if(!skeleton->build(skeletonFilename.toLocal8Bit().data()))
    {
        emit log(LOGERR, QString("Error loading file: %1").arg(skeletonFilename));
        delete skeleton;
        return false;
    }
    
    if (m_mesh->property(skeletonHandle) != NULL) {
        delete m_mesh->property(skeletonHandle);
    }
    m_mesh->property(skeletonHandle) = skeleton;
    
    return true;
}
