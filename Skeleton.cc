//
//  Skeleton.cpp
//  OpenFlipper
//
//  Created by 饒亢 on 2015/1/31.
//
//

#include "PlushPlugin.hh"
#include "SuperDeform/Skeleton.hh"
#include "SuperDeform/Weight.hh"

bool PlushPlugin::loadBoneWeight(int meshId) {
    BaseObjectData *obj;
    PluginFunctions::getObject(meshId, obj);
    if (!obj->dataType(DATA_TRIANGLE_MESH)) {
        emit log(LOGERR, QString("Not a valid TriMesh of object %1").arg(QString::number(meshId)));
        return false;
    }
    
    TriMesh *mesh = PluginFunctions::triMesh(obj);
    QString meshName = QFileInfo(obj->name()).baseName();

    QString weightFilename = meshName+".skeleton.weight";

    if(mesh->property(skeletonHandle) == NULL)
    {
        // No skeleton yet, try to load skeleton
        if (!loadSkeleton(meshId)) {
            return false;
        }
    }
    Skeleton *skeleton = mesh->property(skeletonHandle);
    
    ifstream boneWeightFile(weightFilename.toLocal8Bit().data());
    if(!boneWeightFile.is_open())
    {
        emit log(LOGERR, QString("Error opening file: %1").arg(weightFilename));
        return false;
    }

    // Assign bone weight to vertices
    for(VertexIter v_it = mesh->vertices_begin(); v_it != mesh->vertices_end(); v_it++)
    {
        double *weights = new double[skeleton->bones.size()];
        for(size_t j = 0; j < skeleton->bones.size(); j++)
        {
            double w;
            boneWeightFile >> w;
            weights[j] = w;
        }
        if (mesh->property(bonesWeightHandle, *v_it) != NULL) {
            delete[] mesh->property(bonesWeightHandle, *v_it);
        }
        mesh->property(bonesWeightHandle, *v_it) = weights;
    }
    boneWeightFile.close();

    return true;
}

void PlushPlugin::saveBoneWeight(int meshId) {
    BaseObjectData *obj;
    PluginFunctions::getObject(meshId, obj);
    if (!obj->dataType(DATA_TRIANGLE_MESH)) {
        emit log(LOGERR, QString("Not a valid TriMesh of object %1").arg(QString::number(meshId)));
        return;
    }
    
    TriMesh *mesh = PluginFunctions::triMesh(obj);
    QString meshName = QFileInfo(obj->name()).baseName();

    // Prepare file for saving data
    QFile file(meshName+".skeleton.weight");
    file.open(QIODevice::WriteOnly | QIODevice::Text);
    QTextStream out(&file);
    
    int nBones = mesh->property(skeletonHandle)->bones.size();
    for (VertexIter v_it = mesh->vertices_begin(); v_it != mesh->vertices_end(); v_it++) {
        double *w = mesh->property(bonesWeightHandle, *v_it);
        for (int i = 0; i < nBones; i++) {
            out << w[i] << " ";
        }
        out << endl;
    }
    file.close();
}

bool PlushPlugin::loadSkeleton(int meshId) {
    BaseObjectData *obj;
    PluginFunctions::getObject(meshId, obj);
    if (!obj->dataType(DATA_TRIANGLE_MESH)) {
        emit log(LOGERR, QString("Not a valid TriMesh of object %1").arg(QString::number(meshId)));
        return false;
    }
    
    TriMesh *mesh = PluginFunctions::triMesh(obj);
    QString meshName = QFileInfo(obj->name()).baseName();

    QString skeletonFilename = meshName+".skeleton";
    
    // Add skeleton to mesh property
    Skeleton *skeleton = new Skeleton;
    if(!skeleton->build(skeletonFilename.toLocal8Bit().data()))
    {
        emit log(LOGERR, QString("Error loading file: %1").arg(skeletonFilename));
        delete skeleton;
        return false;
    }
    
    if (mesh->property(skeletonHandle) != NULL) {
        delete mesh->property(skeletonHandle);
    }
    mesh->property(skeletonHandle) = skeleton;
    
    return true;
}