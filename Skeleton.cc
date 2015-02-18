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
#include <ObjectTypes/PolyLine/PolyLine.hh>
//#include <assimp/Importer.hpp> // C++ importer interface
//#include <assimp/scene.h> // Output data structure
//#include <assimp/postprocess.h> // Post processing flags

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
    
//    // -------------------------------------------------------------------
//    /** Describes a bone weight on a vertex */
//    struct Weight
//    {
//        unsigned int mBone; ///< Index of the bone
//        float mWeight;      ///< Weight of that bone on this vertex
//        Weight() { }
//        Weight( unsigned int pBone, float pWeight)
//        {
//            mBone = pBone;
//            mWeight = pWeight;
//        }
//        
//        /** Comparision operator to sort bone weights by descending weight */
//        bool operator < (const Weight& pWeight) const
//        { 
//            return mWeight > pWeight.mWeight;
//        }
//    };
//
//    QString meshFilename = meshName+"_skeleton.dae";
//    QByteArray meshSkeletonFilenameByteArray = meshFilename.toLocal8Bit();
//    const char *meshSkeletonFilenameCStr = meshSkeletonFilenameByteArray.data();
//    Assimp::Importer imp;
//    const aiScene* scene = imp.ReadFile(meshSkeletonFilenameCStr,
//                                        0);
//    if (!scene) {
//        emit log(LOGERR, QString("assimp info: Unable to load input file %1\n").arg(meshFilename));
//        return;
//    }
// 
//    // Find corresponding mesh
//    aiNode *meshNode = NULL;
//    aiNode *armatureNode = NULL;
//    std::vector<std::pair<aiNode*, aiMatrix4x4> > queue;
//    aiMatrix4x4 identity;
//    queue.push_back(std::make_pair(scene->mRootNode, identity));
//    while (queue.size() > 0) {
//        aiNode *currentNode = queue[queue.size()-1].first;
//        aiMatrix4x4 accTransform = queue[queue.size()-1].second;
//        queue.pop_back();
//        if (strcmp(currentNode->mName.C_Str(), meshName.toLocal8Bit().data()) == 0) {
//            meshNode = currentNode;
//        }
//        if (strcmp(currentNode->mName.C_Str(), "Armature") == 0) {
//            armatureNode = currentNode;
//        }
//        
//        if (meshNode && armatureNode) {
////            break;
//        }
//        
//        for (unsigned int i = 0; i < currentNode->mNumChildren; ++i ) {
//            queue.push_back(std::make_pair(currentNode->mChildren[i], accTransform * currentNode->mTransformation));
//        }
//    }
//    
//    if (!meshNode) {
//        emit log(LOGERR, QString("Unable to read file %1").arg(meshFilename));
//        return;
//    }
//    if (!armatureNode) {
//        emit log(LOGERR, QString("Skeleton not found in %1, only supports .fbx generated from Blender.").arg(meshFilename));
//        return;
//    }
//
////    std::map<aiBone*, OpenMesh::Vec3d> transforms;
//    std::map<QString, aiMatrix4x4> transforms;
//    for (unsigned int meshNo = 0; meshNo < meshNode->mNumMeshes; meshNo++) {
//        aiMesh *pMesh = scene->mMeshes[meshNo];
//        for (unsigned int i = 0; i < pMesh->mNumBones; i++) {
//            aiBone *bone = pMesh->mBones[i];
//            if (transforms.find(QString(bone->mName.C_Str())) == transforms.end()) {
//                
//                queue.clear();
//                queue.push_back(std::make_pair(armatureNode, armatureNode->mTransformation.Inverse()));
//                while (queue.size() > 0) {
//                    aiNode *currentNode = queue[queue.size()-1].first;
//                    aiMatrix4x4 accTransform = queue[queue.size()-1].second;
//                    queue.pop_back();
//                    if (strcmp(currentNode->mName.C_Str(), bone->mName.C_Str()) == 0) {
//                        transforms.insert(std::make_pair(QString(bone->mName.C_Str()), accTransform * currentNode->mTransformation * bone->mOffsetMatrix));
//                    } else {
//                        for (unsigned int i = 0; i < currentNode->mNumChildren; ++i ) {
//                            queue.push_back(std::make_pair(currentNode->mChildren[i], accTransform * currentNode->mTransformation));
//                        }
//                    }
//                }
//            }
//        }
//    }
//    
//    // Draw joint orientation
//    for (std::map<QString, aiMatrix4x4>::iterator it = transforms.begin(); it != transforms.end(); it++) {
//        QString boneName = it->first;
//        aiMatrix4x4 transform = it->second;
//        
//        aiVector3t<float> scaling;
//        aiQuaterniont<float> rotation;
//        aiVector3t<float> position;
//        transform.Decompose(scaling, rotation, position);
//        rotation.Normalize();
//        
//        double s = sqrt(1 - rotation.w * rotation.w);
//        OpenMesh::Vec3d direction(rotation.x / s,
//                                  rotation.y / s,
//                                  rotation.z / s);
//        
//        emit log(LOGINFO, QString("(%1: %2 %3 %4)").arg(boneName, 12).arg(   QString::number(direction[0], 'g', 2),
//                                                                             QString::number(direction[1], 'g', 2),
//                                                                             QString::number(direction[2], 'g', 2)));
//        
//        int objId;
//        emit addEmptyObject(DATA_POLY_LINE, objId);
//        PolyLineObject *object = 0;
//        PluginFunctions::getObject(objId, object);
//        object->setName(boneName);
//        PolyLine *polyLine = object->line();
//        
//        polyLine->add_point(PolyLine::Point(scaling.x * position.x,
//                                            scaling.y * position.y,
//                                            scaling.z * position.z));
//        polyLine->add_point(PolyLine::Point(scaling.x * (position.x + direction[0]),
//                                            scaling.y * (position.y + direction[1]),
//                                            scaling.z * (position.z + direction[2])));
//    }
}