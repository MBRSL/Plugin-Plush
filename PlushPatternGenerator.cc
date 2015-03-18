#include "PlushPatternGenerator.hh"

OpenMesh::VPropHandleT<double> PlushPatternGenerator::minCurvatureHandle;
OpenMesh::VPropHandleT<double> PlushPatternGenerator::maxCurvatureHandle;
OpenMesh::VPropHandleT<OpenMesh::Vec3d> PlushPatternGenerator::minCurvatureDirectionHandle;
OpenMesh::VPropHandleT<OpenMesh::Vec3d> PlushPatternGenerator::maxCurvatureDirectionHandle;

OpenMesh::EPropHandleT<double> PlushPatternGenerator::edgeWeightHandle;
OpenMesh::MPropHandleT< std::map<std::pair<VertexHandle, VertexHandle>, double> > PlushPatternGenerator::geodesicDistanceHandle;
OpenMesh::MPropHandleT< std::map<std::pair<VertexHandle, VertexHandle>, std::vector<VertexHandle> > > PlushPatternGenerator::geodesicPathHandle;

OpenMesh::MPropHandleT<Skeleton*> PlushPatternGenerator::skeletonHandle;
OpenMesh::VPropHandleT<double*> PlushPatternGenerator::bonesWeightHandle;

OpenMesh::VPropHandleT<double> PlushPatternGenerator::distortionHandle;

OpenMesh::VPropHandleT<VertexHandle> PlushPatternGenerator::inverseMapping;

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
    
    m_mesh->add_property(distortionHandle, "Flattening distortion");
    
    for (EdgeIter e_it = m_mesh->edges_begin(); e_it != m_mesh->edges_end(); e_it++) {
        m_mesh->property(edgeWeightHandle, *e_it) = -1;
    }
}

void PlushPatternGenerator::uninitProperties() {
    if (m_mesh->property(skeletonHandle) != NULL) {
        delete m_mesh->property(skeletonHandle);
    }
    for (VertexIter v_it = m_mesh->vertices_begin(); v_it != m_mesh->vertices_end(); v_it++) {
        if (m_mesh->property(bonesWeightHandle, *v_it) != NULL) {
            delete[] m_mesh->property(bonesWeightHandle, *v_it);
        }
    }
}

void PlushPatternGenerator::cancelJob() {
    isJobCanceled = true;
}

bool PlushPatternGenerator::getHalfedge(TriMesh *mesh, HalfedgeHandle &heh, VertexHandle from, VertexHandle to) {
    for (TriMesh::VertexEdgeIter ve_it = mesh->ve_iter(from); ve_it; ve_it++)
    {
        // Return halfedge with correct direction
        HalfedgeHandle _heh = mesh->halfedge_handle(*ve_it, 0);
        if (mesh->from_vertex_handle(_heh) == from
            &&  mesh->to_vertex_handle(_heh) == to)
        {
            heh = _heh;
            return true;
        } else if (mesh->from_vertex_handle(_heh) == to
                   &&  mesh->to_vertex_handle(_heh) == from) {
            heh = mesh->halfedge_handle(*ve_it, 1);
            return true;
        }
    }
    return false;
}

bool PlushPatternGenerator::getHalfedge(TriMesh *mesh, HalfedgeHandle &heh, int fromNo, int toNo) {
    VertexHandle from = mesh->vertex_handle(fromNo);
    VertexHandle to = mesh->vertex_handle(toNo);
    return getHalfedge(mesh, heh, from, to);
}

bool PlushPatternGenerator::getEdge(TriMesh *mesh, EdgeHandle &eh, VertexHandle v1, VertexHandle v2) {
    HalfedgeHandle _heh;
    bool result = getHalfedge(mesh, _heh, v1, v2);
    if (result) {
        eh = mesh->edge_handle(_heh);
        return true;
    } else {
        return false;
    }
}

bool PlushPatternGenerator::getEdge(TriMesh *mesh, EdgeHandle &eh, int v1No, int v2No) {
    VertexHandle v1 = mesh->vertex_handle(v1No);
    VertexHandle v2 = mesh->vertex_handle(v2No);
    return getEdge(mesh, eh, v1, v2);
}

OpenMesh::Vec3d PlushPatternGenerator::getVector(TriMesh *mesh, EdgeHandle &_eh) {
    HalfedgeHandle heh = mesh->halfedge_handle(_eh, 0);
    return getVector(mesh, heh);
}

OpenMesh::Vec3d PlushPatternGenerator::getVector(TriMesh *mesh, HalfedgeHandle &_heh) {
    TriMesh::Point p1 = mesh->point(mesh->from_vertex_handle(_heh));
    TriMesh::Point p2 = mesh->point(mesh->to_vertex_handle(_heh));
    return p2-p1;
}

double PlushPatternGenerator::getSumInnerAngle(TriMesh *mesh, HalfedgeHandle heh1, HalfedgeHandle heh2) {
    double sumInnerAngle = 0;

    if (mesh->is_boundary(heh1) || mesh->is_boundary(heh2)) {
        HalfedgeHandle tmp = heh1;
        heh1 = mesh->opposite_halfedge_handle(heh2);
        heh2 = mesh->opposite_halfedge_handle(tmp);
    }
    // Loop through inner edges and sum up inner angles
    HalfedgeHandle current_heh = heh1;
    while (mesh->edge_handle(current_heh) != mesh->edge_handle(heh2)) {
        HalfedgeHandle next_heh = mesh->opposite_halfedge_handle(
                                  mesh->next_halfedge_handle(current_heh));
        OpenMesh::Vec3d vec1 = PlushPatternGenerator::getVector(mesh, current_heh);
        OpenMesh::Vec3d vec2 = PlushPatternGenerator::getVector(mesh, next_heh);
        sumInnerAngle += acos((vec1|vec2) / (vec1.norm()*vec2.norm()));
        
        current_heh = next_heh;
    }
    return sumInnerAngle;
}

/**
 @brief Return a ring of boundary halfedges.
 It will randomly return one ring if mesh contains more than one boundary rings.
 @param boundary Result will be stored here, Ordered halfedges forms a ring.
 @retval false if it contains no boundaries.
 */
bool PlushPatternGenerator::getBoundaryOfOpenedMesh(std::vector<HalfedgeHandle> &boundary, TriMesh &mesh) {
    HalfedgeHandle heh;
    for (HalfedgeIter he_it = mesh.halfedges_begin(); he_it != mesh.halfedges_end(); he_it++) {
        if (mesh.is_boundary(*he_it)) {
            heh = *he_it;
            break;
        }
    }
    
    if (!mesh.is_valid_handle(heh)) {
        return false;
    }
    
    HalfedgeHandle current_heh = heh;
    do {
        boundary.push_back(current_heh);
        current_heh = mesh.next_halfedge_handle(current_heh);
    } while (current_heh != heh);
    
    return true;
}

/**
 * @brief This function returns a set of rings of halfedges. These rings are separated by given edges.
 * @param loops Result loops will be stored here.
 * @param jointVertices (This parameter should be compute internally) The set of intersection point of m_spanningTree.
 */
void PlushPatternGenerator::getBoundariesByEdges(std::vector< std::vector<HalfedgeHandle> > &boundaries, std::vector<EdgeHandle> &separator) {
    std::set<HalfedgeHandle> isVisited;
    for (size_t i = 0; i < separator.size(); i++) {
        HalfedgeHandle heh = m_mesh->halfedge_handle(separator[i], 0);
        VertexHandle startV = m_mesh->from_vertex_handle(heh);
        
        for (TriMesh::VertexOHalfedgeIter voh_it = m_mesh->voh_iter(startV); voh_it; voh_it++) {
            // Skip visited edges and non-spanning-tree edges
            if (std::find(isVisited.begin(), isVisited.end(), *voh_it) != isVisited.end()
            ||  std::find(separator.begin(), separator.end(), m_mesh->edge_handle(*voh_it)) == separator.end()) {
                continue;
            }
            
            // We got a non-visited edge, lets start searching from this edge
            HalfedgeHandle start_heh = *voh_it;
            HalfedgeHandle current_heh = start_heh;
            isVisited.insert(start_heh);
            
            std::vector<HalfedgeHandle> boundary;
            boundary.push_back(current_heh);
            do {
                assert(!m_mesh->is_boundary(current_heh));
                
                HalfedgeHandle next_heh = current_heh;
                // Find next edge which is in spanning tree
                do {
                    // While there are possibly many adjacent edges to choose, we use next->opposite to iterator through them
                    // so that we will visit them all clockwise without missing any one.
                    next_heh = m_mesh->opposite_halfedge_handle(
                                                                m_mesh->next_halfedge_handle(next_heh));
                    
                    // This means there are no other edges except current_heh (which is where you came from)
                    // In this case, just go back and treat this path as a dart
                    if (next_heh == current_heh) {
                        // Do nothing, it will be flipped later.
                    }
                } while (std::find(separator.begin(), separator.end(), m_mesh->edge_handle(next_heh)) == separator.end());
                
                // Flip it to point to the next vertex
                current_heh = m_mesh->opposite_halfedge_handle(next_heh);
                isVisited.insert(current_heh);
                
                boundary.push_back(current_heh);
            } while (current_heh != start_heh);
            
            // Delete last duplicated halfedge
            boundary.pop_back();
            
            boundaries.push_back(boundary);
        }
    }
}