//
//  Geodesic.cc
//  OpenFlipper
//
//  Created by 饒亢 on 2015/1/13.
//
//

#include "PlushPatternGenerator.hh"

#include <MeshTools/MeshSelectionT.hh>
#include <queue>

/**
 @brief Detect if there are any intersections between two paths.
 This function returns true if
 1. One path contains both the starting/ending points of another path.
 2. One path crossover another path.
 
 Note that overlapping (paths with identical edges) is not considered intersection here.
 
 @param pathA vertex ids list of path one
 @param pathB vertex ids list of path two
 @return is the two paths are interected or not
 */
bool PlushPatternGenerator::isIntersected(std::vector<VertexHandle> pathA, std::vector<VertexHandle> pathB) {
    OpenMesh::Vec3d prevCrossVec(0,0,0);
    
    bool AcontainsStartingPointB = false;
    bool AcontainsEndingPointB = false;
    bool BcontainsStartingPointA = false;
    bool BcontainsEndingPointA = false;
    // first check if there are any same vertex (says v) in both path
    for (int i = 0; i < (int)pathA.size(); i++) {
        VertexHandle vA = pathA[i];
        for (int j = 0; j < (int)pathB.size(); j++) {
            VertexHandle vB = pathB[j];
            if (vA == vB) {
                // index out of range, which means v1 or v2 is at the end.
                if (i-1 < 0 || i+1 >= (int)pathA.size()
                ||  j-1 < 0 || j+1 >= (int)pathB.size()) {
                    if (i-1 < 0) {
                        BcontainsStartingPointA = true;
                    }
                    if (i+1 >= (int)pathA.size()) {
                        BcontainsEndingPointA = true;
                    }
                    if (j-1 < 0) {
                        AcontainsStartingPointB = true;
                    }
                    if (j+1 >= (int)pathB.size()) {
                        AcontainsEndingPointB = true;
                    }
                    // Detecting case 1
                    if ((AcontainsStartingPointB && AcontainsEndingPointB)
                    ||  (BcontainsStartingPointA && BcontainsEndingPointA)) {
                        return true;
                    }
                    // If not case 1, consider no crossover on these vertices
                    continue;
                }
                // then compare the sign of the following cross products
                // (edge_before_v_on_path1 x edge_before_v_on_path2)
                // and
                // (edge_after_v_on_path1 x edge_after_v_on_path2)
                VertexHandle vA0 = pathA[i-1];
                VertexHandle vA1 = pathA[i];
                VertexHandle vA2 = pathA[i+1];
                
                VertexHandle vB0 = pathB[j-1];
                VertexHandle vB1 = pathB[j];
                VertexHandle vB2 = pathB[j+1];
                
                OpenMesh::Vec3d vecA01 = m_mesh->point(vA0) - m_mesh->point(vA1);
                OpenMesh::Vec3d vecB01 = m_mesh->point(vB0) - m_mesh->point(vB1);

                OpenMesh::Vec3d vecA21 = m_mesh->point(vA2) - m_mesh->point(vA1);
                OpenMesh::Vec3d vecB21 = m_mesh->point(vB2) - m_mesh->point(vB1);

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
 @brief Generate circular seams that are not possible to get with calcSeams.
 This function find vertices with high curvature and try to form rings of seams.
 The result is append to seamsHandle property of m_mesh.
 @param <#parameter#>
 @return <#retval#>
 @retval <#meaning#>
 */
bool PlushPatternGenerator::calcCircularSeams(TriMesh *mesh) {
    // Find the most distorted vertex, also calculate mean & std
    double minAreaDiff = 1e9;
    double sumCurvature = 0;
    double sumCurvatureSqr = 0;
    VertexHandle minV;
    OpenMesh::VPropHandleT<VertexHandle> inverseMapping = getInverseMappingHandle(mesh);
    for (VertexIter v_it = mesh->vertices_begin(); v_it != mesh->vertices_end(); v_it++) {
        VertexHandle originalV = mesh->property(inverseMapping, *v_it);
        double avgAreaDiff = m_mesh->property(distortionVHandle, originalV);
        
        if (avgAreaDiff < minAreaDiff) {
            minAreaDiff = avgAreaDiff;
            minV = *v_it;
        }
        
        double curvature = m_mesh->property(maxCurvatureHandle, originalV);
        sumCurvature += curvature;
        sumCurvatureSqr += curvature * curvature;
    }
    double mean = sumCurvature/mesh->n_vertices();
    double std = sqrt(sumCurvatureSqr/mesh->n_vertices() - mean*mean);

    minV = *mesh->vertices_begin();
    // Try to search boundary with high curvature from minV
    std::set<VertexHandle> visited;
    std::queue<VertexHandle> queue;
    queue.push(minV);
    
    while (queue.size() > 0) {
        VertexHandle v = queue.front();
        queue.pop();
        visited.insert(v);
        
        for (TriMesh::ConstVertexOHalfedgeIter cvoh_it = mesh->cvoh_iter(v); cvoh_it; cvoh_it++) {
            VertexHandle neighborV = mesh->to_vertex_handle(*cvoh_it);
            VertexHandle originalNeighborV = mesh->property(inverseMapping, neighborV);
            
            // outward expanding
            if (visited.find(neighborV) == visited.end()) {
                visited.insert(neighborV);
                
                double curvature = m_mesh->property(maxCurvatureHandle, originalNeighborV);
                // Only continue exploring from low curvature vertex
                if (curvature < mean + 1.5*std) {
                    queue.push(neighborV);
                }
            }
        }
    }

    // visualize vertices beyond threshold
//    std::vector<int> verticesId;
//    for (VertexIter v_it = mesh->vertices_begin(); v_it != mesh->vertices_end(); v_it++) {
//        VertexHandle originalV = mesh->property(inverseMapping, *v_it);
//        if (m_mesh->property(maxCurvatureHandle, originalV) >= mean + std) {
//            verticesId.push_back(originalV.idx());
//        }
//    }
//    MeshSelection::selectVertices(m_mesh, verticesId);
    
    // Now visited vector should contains some part of the mesh and stopped at high curvature area
    // Extract boundary vertices from visited vector
    std::set<VertexHandle> finalBorderVertices;
    for (auto v_it = visited.begin(); v_it != visited.end(); v_it++) {
        // If this vertex has a non-visited neighbor, it's boundary vertex
        for (TriMesh::ConstVertexVertexIter cvv_it = mesh->cvv_iter(*v_it); cvv_it; cvv_it++) {
            if (visited.find(*cvv_it) == visited.end()) {
                finalBorderVertices.insert(*v_it);
                break;
            }
        }
    }
    
    // Form edges from vertices
    visited.clear();
    OpenMesh::MPropHandleT< std::vector<EdgeHandle> > seamsHandle = getSeamsHandle(m_mesh);
    std::vector<EdgeHandle> &seams = m_mesh->property(seamsHandle);
    for (auto v_it = finalBorderVertices.begin(); v_it != finalBorderVertices.end(); v_it++) {
        VertexHandle startingV = *v_it;
        if (visited.find(startingV) != visited.end()) {
            continue;
        }
        queue.push(startingV);
        visited.insert(startingV);
        
        std::vector<HalfedgeHandle> subBorderEdges;
        while (queue.size() > 0) {
            VertexHandle v = queue.front();
            queue.pop();
            
            for (TriMesh::ConstVertexOHalfedgeIter cvoh_it = mesh->cvoh_iter(v); cvoh_it; cvoh_it++) {
                VertexHandle neighborV = mesh->to_vertex_handle(*cvoh_it);
                
                if (finalBorderVertices.find(neighborV) != finalBorderVertices.end()) {
                    if (visited.find(neighborV) == visited.end()) {
                        visited.insert(neighborV);
                        queue.push(neighborV);
                        
                        subBorderEdges.push_back(*cvoh_it);
                        break;
                    } else if (neighborV == startingV
                           && std::find(subBorderEdges.begin(), subBorderEdges.end(),
                                        mesh->opposite_halfedge_handle(*cvoh_it)) == subBorderEdges.end()) {
                       subBorderEdges.push_back(*cvoh_it);
                       break;
                    }
                }
            }
        }
        
        if (subBorderEdges.size() > 1) {
            VertexHandle endingV = mesh->to_vertex_handle(subBorderEdges[subBorderEdges.size()-1]);
            // Test if this is a ring using last halfedge
            if (startingV == endingV) {
                for (auto he_it = subBorderEdges.begin(); he_it != subBorderEdges.end(); he_it++) {
                    VertexHandle originalV1 = mesh->property(inverseMapping, mesh->from_vertex_handle(*he_it));
                    VertexHandle originalV2 = mesh->property(inverseMapping, mesh->to_vertex_handle(*he_it));
                    EdgeHandle originalEh;
                    assert(getEdge(m_mesh, originalEh, originalV1, originalV2));
                    seams.push_back(originalEh);
                }
            }
        }
    }
    return true;
}

/**
 *  @brief Generate seams by connecting given vertices. Seams are actually steiner tree.
 *  This function will calculate shortest paths between all vertex pairs and then
 *  tries to connect them to form a minimal spanning tree.
 *  This function requires geodesicDistance & geodesicPath of given vertices. They should be calculated before calling this function.
 *  The result is stored in seamsHandle property of m_mesh.
 *
 *  @param selectedVertices Resulting tree will span through these vertices.
 *  @param limitNum (Debugging) If set to postive number, only first num paths are used for spanningTree. 0 means no limitation.
 *  @param elimination (Debugging) If set to true, it will to remove paths which intersected with lower-cost paths.
 *  @param allPath (Debugging) If false, only used (limitNum)-th path for spanningTree.
 *  @return False if error occured.
 */
bool PlushPatternGenerator::calcSeams(std::vector<VertexHandle> selectedVertices,
                                             int limitNum,
                                             bool elimination,
                                             bool allPath) {
    std::map<std::pair<VertexHandle, VertexHandle>, double> &geodesicDistance = m_mesh->property(geodesicDistanceHandle);
    std::map<std::pair<VertexHandle, VertexHandle>, std::vector<VertexHandle> > &geodesicPath = m_mesh->property(geodesicPathHandle);
    
    // create a new graph with selectedVerices, calculate all paths between them
    std::vector< std::pair<double, std::vector<VertexHandle> > > pathCandidates;
    for (size_t i = 0; i < selectedVertices.size(); i++) {
        VertexHandle sourceHandle = selectedVertices[i];

        for (size_t j = i+1; j < selectedVertices.size(); j++) {
            VertexHandle destHandle = selectedVertices[j];
            
            // if result is not found, caculate now
            std::map<std::pair<VertexHandle, VertexHandle>, double>::iterator distanceFound =
            geodesicDistance.find(std::make_pair(sourceHandle, destHandle));
            std::map<std::pair<VertexHandle, VertexHandle>, std::vector<VertexHandle> >::iterator pathFound =
            geodesicPath.find(std::make_pair(sourceHandle, destHandle));

            if (distanceFound == geodesicDistance.end() || pathFound == geodesicPath.end()) {
                emit log(LOGERR, QString("Unreachable from vertex %1 to %2.").arg(sourceHandle.idx(), destHandle.idx()));
                emit log(LOGERR, "Do you have geodesic calculated before?");
                return false;
            }
            double cost = distanceFound->second;
            std::vector<VertexHandle> path = pathFound->second;
            
            // assign edge & weight
            pathCandidates.push_back(std::make_pair(cost, path));
        }
    }

    // Sort candidates from lost cost to high cost
    std::sort(pathCandidates.begin(), pathCandidates.end());
    
    // Choose paths which do not intersect with previous (lower-cost) path.
    std::vector< std::pair<double, std::vector<VertexHandle> > > result;
    if (elimination) {
        for (std::vector<std::pair<double, std::vector<VertexHandle> > >::iterator it = pathCandidates.begin(); it != pathCandidates.end(); it++) {
            bool noIntersection = true;
            std::vector<VertexHandle> path1 = it->second;
            
            for (std::vector<std::pair<double, std::vector<VertexHandle> > >::iterator it2 = result.begin(); it2 != result.end(); it2++) {
                std::vector<VertexHandle> path2 = it2->second;
                if (isIntersected(path1, path2)) {
                    noIntersection = false;
                    break;
                }
            }
            if (noIntersection) {
                result.push_back(*it);
            }
        }
    } else {
        result = pathCandidates;
    }
    
    // insert edges into seams
    OpenMesh::MPropHandleT< std::vector<EdgeHandle> > seamsHandle = getSeamsHandle(m_mesh);
    std::vector<EdgeHandle> &seams = m_mesh->property(seamsHandle);
    seams.clear();

    int count = 0;
    for (std::vector<std::pair<double, std::vector<VertexHandle> > >::iterator it = result.begin(); it != result.end(); it++, count++) {
        // Break if we reach limitNum
        if (count >= limitNum && limitNum != 0) {
            break;
        }
        // Skip if this path is not what we want
        if (!allPath && count != limitNum-1) {
            continue;
        }
        std::vector<VertexHandle> path = it->second;
        
        // There should be at least two vertices on a path.
        assert(path.size() > 1);
        
        for (size_t i = 1; i < path.size(); i++) {
            EdgeHandle eh;
            assert(getEdge(m_mesh, eh, path[i-1], path[i]));
            seams.push_back(eh);
        }
        
        QString msg = QString("Weight of path #%1 from %2 to %3: %4").arg(count+1).arg((path.begin())->idx()).arg((path.end()-1)->idx()).arg(it->first);
        emit log(LOGINFO, msg);
    }
    
    return true;
}
