#include "PlushPatternGenerator.hh"

#include <queue>

OpenMesh::VPropHandleT<double> PlushPatternGenerator::minCurvatureHandle;
OpenMesh::VPropHandleT<double> PlushPatternGenerator::maxCurvatureHandle;
OpenMesh::VPropHandleT<OpenMesh::Vec3d> PlushPatternGenerator::minCurvatureDirectionHandle;
OpenMesh::VPropHandleT<OpenMesh::Vec3d> PlushPatternGenerator::maxCurvatureDirectionHandle;

OpenMesh::EPropHandleT<double> PlushPatternGenerator::edgeWeightHandle;
OpenMesh::MPropHandleT< std::map<std::pair<VertexHandle, VertexHandle>, double> > PlushPatternGenerator::geodesicDistanceHandle;
OpenMesh::MPropHandleT< std::map<std::pair<VertexHandle, VertexHandle>, std::vector<VertexHandle> > > PlushPatternGenerator::geodesicPathHandle;

OpenMesh::MPropHandleT<Skeleton*> PlushPatternGenerator::skeletonHandle;
OpenMesh::VPropHandleT<double*> PlushPatternGenerator::bonesWeightHandle;

OpenMesh::MPropHandleT< std::vector<TriMesh> > PlushPatternGenerator::flattenedMeshesHandle;
OpenMesh::FPropHandleT<double> PlushPatternGenerator::distortionFHandle;
OpenMesh::VPropHandleT<double> PlushPatternGenerator::distortionVHandle;

OpenMesh::VPropHandleT<VertexHandle> PlushPatternGenerator::getInverseMappingHandle(TriMesh *mesh) {
    OpenMesh::VPropHandleT<VertexHandle> inverseMappingHandle;
    if (!mesh->get_property_handle(inverseMappingHandle, "inverseMapping")) {
        mesh->add_property(inverseMappingHandle, "inverseMapping");
    }
    return inverseMappingHandle;
}

OpenMesh::MPropHandleT< std::vector<EdgeHandle> > PlushPatternGenerator::getSeamsHandle(TriMesh *mesh) {
    OpenMesh::MPropHandleT< std::vector<EdgeHandle> > seamsHandle;
    if(!mesh->get_property_handle(seamsHandle, "seams")) {
        mesh->add_property(seamsHandle, "seams");
    }
    return seamsHandle;
}

PlushPatternGenerator::PlushPatternGenerator(TriMesh *mesh, QString meshName) : m_mesh(mesh), m_meshName(meshName) {
    // convert m_mesh from OpenFlipper to CGAL
    CGAL_Polyhedron_helper::convert_OpenMesh_to_CGAL(m_mesh, m_polyhedron, m_verticesMapping);
    
    // setup properties
    initProperties();
    
    loadCurvature();
    
    loadSkeleton();
    
    loadBoneWeight();
    
    loadGeodesic();
}

PlushPatternGenerator::~PlushPatternGenerator() {
    uninitProperties();
}

void PlushPatternGenerator::initProperties() {
    m_mesh->add_property(minCurvatureHandle, "Min Curvature");
    m_mesh->add_property(maxCurvatureHandle, "Max Curvature");
    m_mesh->add_property(minCurvatureDirectionHandle, "Min curvature direction");
    m_mesh->add_property(maxCurvatureDirectionHandle, "Max curvature direction");
    
    m_mesh->add_property(edgeWeightHandle, "Edge weight");
    m_mesh->add_property(geodesicDistanceHandle, "Geodesic distance between vertices pair");
    m_mesh->add_property(geodesicPathHandle, "Geodesic path between vertices pair");
    
    m_mesh->add_property(skeletonHandle, "Skeleton");
    m_mesh->add_property(bonesWeightHandle, "Bone weights for each vertex");
    
    m_mesh->add_property(flattenedMeshesHandle, "flattened meshes");

    m_mesh->add_property(distortionVHandle, "Flattening distortion of vertices");
    m_mesh->add_property(distortionFHandle, "Flattening distortion of faces");
    
    for (EdgeIter e_it = m_mesh->edges_begin(); e_it != m_mesh->edges_end(); e_it++) {
        m_mesh->property(edgeWeightHandle, *e_it) = -1;
    }
}

void PlushPatternGenerator::uninitProperties() {
    m_mesh->remove_property(minCurvatureHandle);
    m_mesh->remove_property(maxCurvatureHandle);
    m_mesh->remove_property(minCurvatureDirectionHandle);
    m_mesh->remove_property(maxCurvatureDirectionHandle);
    
    m_mesh->remove_property(edgeWeightHandle);
    m_mesh->remove_property(geodesicDistanceHandle);
    m_mesh->remove_property(geodesicPathHandle);
    
    if (m_mesh->property(skeletonHandle) != NULL) {
        delete m_mesh->property(skeletonHandle);
    }
    for (VertexIter v_it = m_mesh->vertices_begin(); v_it != m_mesh->vertices_end(); v_it++) {
        if (m_mesh->property(bonesWeightHandle, *v_it) != NULL) {
            delete[] m_mesh->property(bonesWeightHandle, *v_it);
        }
    }
    m_mesh->remove_property(skeletonHandle);
    m_mesh->remove_property(bonesWeightHandle);
    
    OpenMesh::MPropHandleT< std::vector<EdgeHandle> > seamsHandle;
    if(!m_mesh->get_property_handle(seamsHandle, "seams")) {
        m_mesh->remove_property(seamsHandle);
    }
    m_mesh->remove_property(flattenedMeshesHandle);

    m_mesh->remove_property(distortionVHandle);
    m_mesh->remove_property(distortionFHandle);
    
}

void PlushPatternGenerator::cancelJob() {
    isJobCanceled = true;
}

bool PlushPatternGenerator::getHalfedge(const TriMesh *mesh, HalfedgeHandle &heh, VertexHandle from, VertexHandle to) {
    for (TriMesh::ConstVertexEdgeIter cve_it = mesh->cve_iter(from); cve_it; cve_it++)
    {
        // Return halfedge with correct direction
        HalfedgeHandle _heh = mesh->halfedge_handle(*cve_it, 0);
        if (mesh->from_vertex_handle(_heh) == from
            &&  mesh->to_vertex_handle(_heh) == to)
        {
            heh = _heh;
            return true;
        } else if (mesh->from_vertex_handle(_heh) == to
                   &&  mesh->to_vertex_handle(_heh) == from) {
            heh = mesh->halfedge_handle(*cve_it, 1);
            return true;
        }
    }
    return false;
}

bool PlushPatternGenerator::getHalfedge(const TriMesh *mesh, HalfedgeHandle &heh, int fromNo, int toNo) {
    VertexHandle from = mesh->vertex_handle(fromNo);
    VertexHandle to = mesh->vertex_handle(toNo);
    return getHalfedge(mesh, heh, from, to);
}

bool PlushPatternGenerator::getEdge(const TriMesh *mesh, EdgeHandle &eh, VertexHandle v1, VertexHandle v2) {
    HalfedgeHandle _heh;
    bool result = getHalfedge(mesh, _heh, v1, v2);
    if (result) {
        eh = mesh->edge_handle(_heh);
        return true;
    } else {
        return false;
    }
}

bool PlushPatternGenerator::getEdge(const TriMesh *mesh, EdgeHandle &eh, int v1No, int v2No) {
    VertexHandle v1 = mesh->vertex_handle(v1No);
    VertexHandle v2 = mesh->vertex_handle(v2No);
    return getEdge(mesh, eh, v1, v2);
}

bool PlushPatternGenerator::getFace(const TriMesh *mesh, FaceHandle &fh, VertexHandle v1, VertexHandle v2, VertexHandle v3) {
    for (TriMesh::ConstVertexFaceIter cvf_it = mesh->cvf_begin(v1); cvf_it; cvf_it++) {
        HalfedgeHandle heh = mesh->halfedge_handle(*cvf_it);
        VertexHandle _v1 = mesh->from_vertex_handle(heh);
        VertexHandle _v2 = mesh->to_vertex_handle(heh);
        VertexHandle _v3 = mesh->to_vertex_handle(mesh->next_halfedge_handle(heh));

        if ((_v1 == v1 && _v2 == v2 && _v3 == v3)
        ||  (_v1 == v1 && _v2 == v3 && _v3 == v2)
        ||  (_v1 == v2 && _v2 == v1 && _v3 == v3)
        ||  (_v1 == v2 && _v2 == v3 && _v3 == v1)
        ||  (_v1 == v3 && _v2 == v1 && _v3 == v2)
        ||  (_v1 == v3 && _v2 == v2 && _v3 == v1)) {
            fh = *cvf_it;
            return true;
        }
    }
    
    // Not found
    return false;
}

std::vector<TriMesh>* PlushPatternGenerator::getFlattenedMeshes() {
    if (!flattenedMeshesHandle.is_valid() || m_mesh->property(flattenedMeshesHandle).size() == 0) {
        return NULL;
    } else {
        return &m_mesh->property(flattenedMeshesHandle);
    }
}

std::vector<EdgeHandle>* PlushPatternGenerator::getSeams() {
    OpenMesh::MPropHandleT< std::vector<EdgeHandle> > seamsHandle = getSeamsHandle(m_mesh);
    if (m_mesh->property(seamsHandle).size() == 0) {
        return NULL;
    } else {
        return &m_mesh->property(seamsHandle);
    }
}

double PlushPatternGenerator::getSumInnerAngle(const TriMesh *mesh, HalfedgeHandle heh1, HalfedgeHandle heh2) {
    assert(mesh->to_vertex_handle(heh1) == mesh->from_vertex_handle(heh2));
    assert(!mesh->is_boundary(heh1) && !mesh->is_boundary(heh2));

    double sumInnerAngle = 0;
    
    // Loop through inner edges and sum up inner angles
    HalfedgeHandle current_heh = heh1;
    while (mesh->edge_handle(current_heh) != mesh->edge_handle(heh2)) {
        HalfedgeHandle next_heh = mesh->opposite_halfedge_handle(
                                                                 mesh->next_halfedge_handle(current_heh));
        sumInnerAngle += mesh->calc_sector_angle(current_heh);
        
        current_heh = next_heh;
    }
    return sumInnerAngle;
}

/**
 @brief Return rings of boundary halfedges. May contains multiple rings
 @param boundary Result will be stored here, Ordered halfedges forms a ring.
 @param getInteriorHalfedge If true, return the opposite halfedge of boundary halfedge. These halfedge are not boundary halfedge.
 @retval false if it contains no boundaries.
 */
bool PlushPatternGenerator::getBoundaryOfOpenedMesh(std::vector< std::vector<HalfedgeHandle> > &boundaries, const TriMesh *mesh, bool getInteriorHalfedge) {
    std::set<HalfedgeHandle> visited;
    for (HalfedgeIter he_it = mesh->halfedges_begin(); he_it != mesh->halfedges_end(); he_it++) {
        if (!mesh->is_boundary(*he_it)
        ||  visited.find(*he_it) != visited.end()) {
            continue;
        }
        
        std:vector<HalfedgeHandle> boundary;
        std::queue<HalfedgeHandle> queue;
        queue.push(*he_it);
        visited.insert(*he_it);
        
        while (!queue.empty()) {
            HalfedgeHandle current_heh = queue.front();
            queue.pop();
            
            boundary.push_back(current_heh);
            
            HalfedgeHandle next_heh = mesh->next_halfedge_handle(current_heh);
            if (visited.find(next_heh) == visited.end()) {
                queue.push(next_heh);
                visited.insert(next_heh);
            }
        }
        
        assert(boundary.size() > 0);
        if (mesh->from_vertex_handle(*he_it) != mesh->to_vertex_handle(boundary[boundary.size()-1])) {
            // Not a loop
            continue;
        }
        
        // The lopp of opposite halfedges is in reversed direction
        if (getInteriorHalfedge) {
            std::reverse(boundary.begin(), boundary.end());
            for (size_t i = 0; i < boundary.size(); i++) {
                boundary[i] = mesh->opposite_halfedge_handle(boundary[i]);
            }
        }
        
        boundaries.push_back(boundary);
    }
    
    
    return boundaries.empty();
}

/**
 * @brief This function returns a set of rings of halfedges. These rings are separated by given edges.
 * @param loops Result loops will be stored here.
 * @param jointVertices (This parameter should be compute internally) The set of intersection point of m_spanningTree.
 */
//void PlushPatternGenerator::getBoundariesByEdges(std::vector< std::vector<HalfedgeHandle> > *boundaries, std::vector<EdgeHandle> *separator) {
//    std::set<HalfedgeHandle> isVisited;
//    for (size_t i = 0; i < separator->size(); i++) {
//        HalfedgeHandle heh = m_mesh->halfedge_handle(separator->at(i), 0);
//        VertexHandle startV = m_mesh->from_vertex_handle(heh);
//        
//        for (TriMesh::VertexOHalfedgeIter voh_it = m_mesh->voh_iter(startV); voh_it; voh_it++) {
//            // Skip visited edges and non-spanning-tree edges
//            if (std::find(isVisited.begin(), isVisited.end(), *voh_it) != isVisited.end()
//            ||  std::find(separator->begin(), separator->end(), m_mesh->edge_handle(*voh_it)) == separator->end()) {
//                continue;
//            }
//            
//            // We got a non-visited edge, lets start searching from this edge
//            HalfedgeHandle start_heh = *voh_it;
//            HalfedgeHandle current_heh = start_heh;
//            isVisited.insert(start_heh);
//            
//            std::vector<HalfedgeHandle> boundary;
//            boundary.push_back(current_heh);
//            do {
//                assert(!m_mesh->is_boundary(current_heh));
//                
//                HalfedgeHandle next_heh = current_heh;
//                // Find next edge which is in spanning tree
//                do {
//                    // While there are possibly many adjacent edges to choose, we use next->opposite to iterator through them
//                    // so that we will visit them all clockwise without missing any one.
//                    next_heh = m_mesh->opposite_halfedge_handle(
//                                                                m_mesh->next_halfedge_handle(next_heh));
//                    
//                    // This means there are no other edges except current_heh (which is where you came from)
//                    // In this case, just go back and treat this path as a dart
//                    if (next_heh == current_heh) {
//                        // Do nothing, it will be flipped later.
//                    }
//                } while (std::find(separator->begin(), separator->end(), m_mesh->edge_handle(next_heh)) == separator->end());
//                
//                // Flip it to point to the next vertex
//                current_heh = m_mesh->opposite_halfedge_handle(next_heh);
//                isVisited.insert(current_heh);
//                
//                boundary.push_back(current_heh);
//            } while (current_heh != start_heh);
//            
//            // Delete last duplicated halfedge
//            boundary.pop_back();
//            
//            boundaries->push_back(boundary);
//        }
//    }
//}

/**
 Split existing m_mesh into sub meshes with boundary as separator.
 @param loops The given set of boundaries. Each boundary should be a loop.
 @param subMeshes Resulting sub meshes will be stored here.
 */
bool PlushPatternGenerator::splitWithBoundary(std::vector<TriMesh> *subMeshes)
{
    // construct set for efficient searching
    OpenMesh::MPropHandleT< std::vector<EdgeHandle> > seamsHandle = getSeamsHandle(m_mesh);
    std::set<EdgeHandle> seamsEdgeSet(m_mesh->property(seamsHandle).begin(), m_mesh->property(seamsHandle).end());
    assert(!seamsEdgeSet.empty() && "seams should not be empty.");
    
    std::set<VertexHandle> seamsVertexSet;
    for (auto e_it = m_mesh->property(seamsHandle).begin(); e_it != m_mesh->property(seamsHandle).end(); e_it++) {
        HalfedgeHandle heh = m_mesh->halfedge_handle(*e_it, 0);
        seamsVertexSet.emplace(m_mesh->from_vertex_handle(heh));
        seamsVertexSet.emplace(m_mesh->to_vertex_handle(heh));
    }

    // Segment mesh into sub meshes using BFS.
    std::set<FaceHandle> visitedF;
    for (FaceIter f_it = m_mesh->faces_begin(); f_it != m_mesh->faces_end(); f_it++) {
        if (visitedF.find(*f_it) != visitedF.end()) {
            continue;
        }
        
        TriMesh newMesh;
        OpenMesh::VPropHandleT<VertexHandle> inverseMapping = getInverseMappingHandle(&newMesh);

        // For non-boundary vertices, we construct a one-to-one mapping from old to new
        std::map<VertexHandle, VertexHandle> vertexToVertexMapping;
        // But for boundary vertices, one old vertex maps to two new vertices. This is the problem.
        // So we use halfedge-vertex mapping because two old halfedges can map to two new vertices.
        std::map<HalfedgeHandle, VertexHandle> halfedgeToVertexMapping;
        
        std::queue<FaceHandle> queue;
        std::set<FaceHandle> connectedFaces;

        queue.push(*f_it);
        visitedF.insert(*f_it);
        connectedFaces.insert(*f_it);
        
        // Explore & insert new vertices
        while (!queue.empty()) {
            FaceHandle f = queue.front();
            queue.pop();
            connectedFaces.insert(f);
            
            HalfedgeHandle heh[3];
            heh[0] = m_mesh->halfedge_handle(f);
            heh[1] = m_mesh->next_halfedge_handle(heh[0]);
            heh[2] = m_mesh->next_halfedge_handle(heh[1]);
            
            for (int i = 0; i < 3; i++) {
                VertexHandle vi = m_mesh->from_vertex_handle(heh[i]);
                // Insert vertices into new mesh
                // If this is a boundary vertex...
                if (seamsVertexSet.find(vi) != seamsVertexSet.end()) {
                    // And only add new vertex if this is a boundary halfedge
                    if (seamsEdgeSet.find(m_mesh->edge_handle(heh[i])) != seamsEdgeSet.end()) {
                        VertexHandle newV = newMesh.add_vertex(m_mesh->point(vi));
                        newMesh.property(inverseMapping, newV) = vi;
                        halfedgeToVertexMapping.emplace(heh[i], newV);
                    }
                }
                // If not a boundary vertex, we have to check if this vertex is already added
                else if (vertexToVertexMapping.find(vi) == vertexToVertexMapping.end()) {
                    VertexHandle newV = newMesh.add_vertex(m_mesh->point(vi));
                    newMesh.property(inverseMapping, newV) = vi;
                    vertexToVertexMapping.emplace(vi, newV);
                }
                
                // Explore other faces through non-boundary halfedges
                FaceHandle neighborFi = m_mesh->opposite_face_handle(heh[i]);
                if (seamsEdgeSet.find(m_mesh->edge_handle(heh[i])) == seamsEdgeSet.end()
                &&  visitedF.find(neighborFi) == visitedF.end()) {
                    queue.push(neighborFi);
                    visitedF.insert(neighborFi);
                }
            }
        }
        
        // Insert new faces
        newMesh.request_face_colors();
        for (std::set<FaceHandle>::iterator f_it = connectedFaces.begin(); f_it != connectedFaces.end(); f_it++) {
            HalfedgeHandle heh[3];
            heh[0] = m_mesh->halfedge_handle(*f_it);
            heh[1] = m_mesh->next_halfedge_handle(heh[0]);
            heh[2] = m_mesh->next_halfedge_handle(heh[1]);
            
            // Search for corresponding new vertex
            VertexHandle newV[3];
            for (int i = 0; i < 3; i++) {
                VertexHandle vi = m_mesh->from_vertex_handle(heh[i]);
                // Non-boundary vertex
                if (seamsVertexSet.find(vi) == seamsVertexSet.end()) {
                    newV[i] = vertexToVertexMapping[vi];
                }
                // boundary vertex on boundary halfedge
                else if (seamsEdgeSet.find(m_mesh->edge_handle(heh[i])) != seamsEdgeSet.end()) {
                    newV[i] = halfedgeToVertexMapping[heh[i]];
                }
                // boundary vertex on non-boundary halfedge, we need to search for neighbor boundary halfedge
                else {
                    // loop until we find a boundary halfedge
                    HalfedgeHandle current_heh = m_mesh->next_halfedge_handle(m_mesh->opposite_halfedge_handle(heh[i]));
                    for (unsigned int count = 0; seamsEdgeSet.find(m_mesh->edge_handle(current_heh)) == seamsEdgeSet.end(); count++) {
                        current_heh = m_mesh->next_halfedge_handle(m_mesh->opposite_halfedge_handle(current_heh));
                        assert(count < m_mesh->valence(vi));
                    }
                    newV[i] = halfedgeToVertexMapping[current_heh];
                }
            }
            
            std::vector<VertexHandle> face_vhs;
            face_vhs.push_back(newV[0]);
            face_vhs.push_back(newV[1]);
            face_vhs.push_back(newV[2]);
            FaceHandle newF = newMesh.add_face(face_vhs);
            newMesh.set_color(newF, m_mesh->color(*f_it));
        }
        subMeshes->push_back(newMesh);
    }
    return true;
}
