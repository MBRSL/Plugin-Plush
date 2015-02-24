//
//  Geodesic.cc
//  OpenFlipper
//
//  Created by 饒亢 on 2015/1/13.
//
//

#include <queue>
#include "PlushPlugin.hh"
#include "CGAL_Polyhedron_builder.hh"

/**
 @brief Detect if there are any intersection between two paths.
 Note that overlapping (paths with identical edges) is not considered intersection here.
 
 @param mesh input mesh
 @param pathA vertex ids list of path one
 @param pathB vertex ids list of path two
 @return is the two paths are interected or not
 */
bool PlushPlugin::isIntersected(TriMesh *mesh, IdList pathA, IdList pathB) {
    OpenMesh::Vec3d prevCrossVec(0,0,0);
    // first check if there are any same vertex (says v) in both path
    for (int i = 0; i < (int)pathA.size(); i++) {
        VertexHandle vA = mesh->vertex_handle(pathA[i]);
        for (int j = 0; j < (int)pathB.size(); j++) {
            VertexHandle vB = mesh->vertex_handle(pathB[j]);
            if (vA == vB) {
                // index out of range, which means v1 or v2 is at the end.
                // consider no crossover at this vertex
                if (i-1 < 0 || i+1 >= (int)pathA.size()
                ||  j-1 < 0 || j+1 >= (int)pathB.size()) {
                    continue;
                }
                // then compare the sign of the following cross products
                // (edge_before_v_on_path1 x edge_before_v_on_path2)
                // and
                // (edge_after_v_on_path1 x edge_after_v_on_path2)
                VertexHandle vA0 = mesh->vertex_handle(pathA[i-1]);
                VertexHandle vA1 = mesh->vertex_handle(pathA[i]);
                VertexHandle vA2 = mesh->vertex_handle(pathA[i+1]);
                
                VertexHandle vB0 = mesh->vertex_handle(pathB[j-1]);
                VertexHandle vB1 = mesh->vertex_handle(pathB[j]);
                VertexHandle vB2 = mesh->vertex_handle(pathB[j+1]);
                
                OpenMesh::Vec3d vecA01 = mesh->point(vA0) - mesh->point(vA1);
                OpenMesh::Vec3d vecB01 = mesh->point(vB0) - mesh->point(vB1);

                OpenMesh::Vec3d vecA21 = mesh->point(vA2) - mesh->point(vA1);
                OpenMesh::Vec3d vecB21 = mesh->point(vB2) - mesh->point(vB1);

                OpenMesh::Vec3d cross1 = vecA01 % vecB01;
                OpenMesh::Vec3d cross2 = vecA21 % vecB21;
                
                // if one cross is 0, means it's overlapping
                if ((cross1 | cross2) == 0) {
                    // record non 0 vector
                    if (cross1.norm() != 0) {
                        prevCrossVec = cross1;
                    } else if (cross2.norm() == 0) {
                        prevCrossVec = cross2;
                    }
                    continue;
                }
                // if their signs are not the same same, no crossover
                if ((cross1 | cross2) < 0) {
                    continue;
                }
                // otherwise, it must be crossover
                else {
                    return true;
                }
            }
        }
    }
    return false;
}

/**
 *  @brief Generate spanning tree (steiner tree indeed) from given vertices.
 *  This function will calculate shortest paths between all vertex pairs and then
 *  tries to connect them to form a minimal spanning tree.
 *
 *  @param _jobId Thread identifier
 *  @param meshId Target mesh id
 *  @param spanningTree Vector of paths with cost in this spanning tree. Path are composed with vertex ids.
 *  @param selectedVertices Resulting tree will span through these vertices.
 *  @param limitNum (Debugging) If set to postive number, only first num paths are used for spanningTree. 0 means no limitation.
 *  @param allPath (Debugging) If false, only used (limitNum)-th path for spanningTree.
 *  @return False if error occured.
 */
bool PlushPlugin::calcSpanningTree(QString _jobId,
                                   int meshId,
                                   std::vector<EdgeHandle> &spanningTree,
                                   IdList selectedVertices,
                                   int limitNum = 0,
                                   bool allPath = true) {
    BaseObjectData *obj;
    PluginFunctions::getObject(meshId, obj);
    if (!obj->dataType(DATA_TRIANGLE_MESH)) {
        emit log(LOGERR, QString("Not a valid TriMesh of object %1").arg(QString::number(meshId)));
        return false;
    }
    
    isJobCanceled = false;
    
    TriMesh *mesh = PluginFunctions::triMesh(obj);
    QString meshName = QFileInfo(obj->name()).baseName();
    emit setJobDescription(_jobId, QString("Calculating curvature: %1").arg(meshName));
    
    std::map<std::pair<VertexHandle, VertexHandle>, double> &geodesicDistance = mesh->property(geodesicDistanceHandle);
    std::map<std::pair<VertexHandle, VertexHandle>, IdList> &geodesicPath = mesh->property(geodesicPathHandle);
    
    // convert mesh from OpenFlipper to CGAL
    Polyhedron P;
    std::map<int, boost_vertex_descriptor> verticesMapping;
    CGAL_Polyhedron_helper::convert_OpenMesh_to_CGAL(mesh, P, verticesMapping);
    
    // create a new graph with selectedVerices, calculate all paths between them
    std::vector< std::pair<double, IdList> > result;
    for (size_t i = 0; i < selectedVertices.size(); i++) {
        for (size_t j = i+1; j < selectedVertices.size(); j++) {
            if (isJobCanceled) {
                emit log(LOGINFO, "Geodesic calculation canceled.");
                return false;
            }
            
            int sourceIdx = selectedVertices[i];
            int destIdx = selectedVertices[j];
            VertexHandle sourceHandle = mesh->vertex_handle(sourceIdx);
            VertexHandle destHandle = mesh->vertex_handle(destIdx);
            
            // if result is not found, caculate now
            std::map<std::pair<VertexHandle, VertexHandle>, double>::iterator distanceFound =
            geodesicDistance.find(std::make_pair(sourceHandle, destHandle));
            std::map<std::pair<VertexHandle, VertexHandle>, IdList>::iterator pathFound =
            geodesicPath.find(std::make_pair(sourceHandle, destHandle));
            
            if (distanceFound == geodesicDistance.end() || pathFound == geodesicPath.end()) {
                calcGeodesic(mesh, P, verticesMapping, sourceHandle, selectedVertices);
                distanceFound = geodesicDistance.find(std::make_pair(sourceHandle, destHandle));
                pathFound = geodesicPath.find(std::make_pair(sourceHandle, destHandle));
            }
            double cost = distanceFound->second;
            IdList path = pathFound->second;
            
            // assign edge & weight
            result.push_back(std::make_pair(cost, path));
        }
        // most of the time is spent on geodesic calculation
        int status = (double)(i+1)/selectedVertices.size() * 100;
        emit setJobState(_jobId, status);
    }
    
    // remove paths which intersect with previous (lower-cost) path.
    std::sort(result.begin(), result.end());
    bool modified = false;
    do {
        modified = false;
        for (std::vector<std::pair<double, IdList> >::iterator it = result.begin(); it != result.end() && ! modified; it++) {
            for (std::vector<std::pair<double, IdList> >::iterator it2 = result.begin()+1; it2 != result.end() && !modified; it2++) {
                if (isIntersected(mesh, it->second, it2->second)) {
                    result.erase(it2);
                    modified = true;
                    break;
                }
            }
        }
    } while (modified);
    
    // insert edges into spanning tree
    int count = 0;
    for (std::vector<std::pair<double, IdList> >::iterator it = result.begin(); it != result.end() && ! modified; it++, count++) {
        // Break if we reach limitNum
        if (count >= limitNum && limitNum != 0) {
            break;
        }
        // Skip if this path is not what we want
        if (!allPath && count != limitNum-1) {
            continue;
        }
        IdList path = it->second;
        
        // There should be at least two vertices on a path.
        assert(path.size() > 1);
        
        for (int i = 1; i < (int)path.size(); i++) {
            EdgeHandle eh;
            assert(getEdge(mesh, eh, path[i-1], path[i]));
            spanningTree.push_back(eh);
        }
        
        QString msg = QString("Total weight of path #%1: %2").arg(count+1).arg(it->first);
        emit log(LOGINFO, msg);
    }
    
    return true;
}
