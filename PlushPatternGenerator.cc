#include "PlushPatternGenerator.hh"

#include <queue>

OpenMesh::EPropHandleT<int> PlushPatternGenerator::segment_no_handle;

OpenMesh::MPropHandleT<int> PlushPatternGenerator::merge_iterations_handle;

OpenMesh::MPropHandleT< std::map< std::vector<HalfedgeHandle>, double> > PlushPatternGenerator::joint_boundary_area_handle;
OpenMesh::MPropHandleT< std::map< std::vector<HalfedgeHandle>, double> > PlushPatternGenerator::joint_boundary_distortion_handle;

OpenMesh::MPropHandleT< Eigen::SparseMatrix<double> > PlushPatternGenerator::conformalMatrixHandle;

OpenMesh::VPropHandleT<double> PlushPatternGenerator::maxCurvatureHandle;

OpenMesh::EPropHandleT<double> PlushPatternGenerator::edgeWeightHandle;
OpenMesh::MPropHandleT< std::map<std::pair<VertexHandle, VertexHandle>, double> > PlushPatternGenerator::geodesicDistanceHandle;
OpenMesh::MPropHandleT< std::map<std::pair<VertexHandle, VertexHandle>, std::vector<VertexHandle> > > PlushPatternGenerator::geodesicPathHandle;

OpenMesh::MPropHandleT<Skeleton*> PlushPatternGenerator::skeletonHandle;
OpenMesh::VPropHandleT<double*> PlushPatternGenerator::bonesWeightHandle;

OpenMesh::MPropHandleT< std::vector<TriMesh> > PlushPatternGenerator::flattenedMeshesHandle;
OpenMesh::FPropHandleT<double> PlushPatternGenerator::distortionFHandle;
OpenMesh::VPropHandleT<double> PlushPatternGenerator::distortionVHandle;

OpenMesh::VPropHandleT<VertexHandle> PlushPatternGenerator::getInverseMappingHandle(TriMesh *mesh) {
    OpenMesh::VPropHandleT<VertexHandle> inverse_mappingHandle;
    if (!mesh->get_property_handle(inverse_mappingHandle, "inverse_mapping")) {
        mesh->add_property(inverse_mappingHandle, "inverse_mapping");
    }
    return inverse_mappingHandle;
}

OpenMesh::MPropHandleT< std::set<EdgeHandle> > PlushPatternGenerator::getSeamsHandle(TriMesh *mesh) {
    OpenMesh::MPropHandleT< std::set<EdgeHandle> > seamsHandle;
    if(!mesh->get_property_handle(seamsHandle, "seams")) {
        mesh->add_property(seamsHandle, "seams");
    }
    return seamsHandle;
}

PlushPatternGenerator::PlushPatternGenerator(TriMesh *mesh, QString meshName) : m_mesh(mesh), m_meshName(meshName) {
    m_mesh->request_edge_colors();
    
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
    m_mesh->add_property(segment_no_handle, "Segment number");
    m_mesh->add_property(merge_iterations_handle, "Merge iterations");
    
    m_mesh->add_property(joint_boundary_area_handle, "Sum of area for joint boundary");
    m_mesh->add_property(joint_boundary_distortion_handle, "Distortion for joint boundary");
    
    m_mesh->add_property(conformalMatrixHandle, "Conformal sparse matrix");
    
    m_mesh->add_property(maxCurvatureHandle, "Max Curvature");
    
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
    m_mesh->remove_property(segment_no_handle);
    
    m_mesh->remove_property(merge_iterations_handle);

    m_mesh->remove_property(joint_boundary_area_handle);
    m_mesh->remove_property(joint_boundary_distortion_handle);

    m_mesh->remove_property(conformalMatrixHandle);
    
    m_mesh->remove_property(maxCurvatureHandle);
    
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

std::set<EdgeHandle>* PlushPatternGenerator::getSeams() {
    OpenMesh::MPropHandleT< std::set<EdgeHandle> > seamsHandle = getSeamsHandle(m_mesh);
    if (m_mesh->property(seamsHandle).size() == 0) {
        return NULL;
    } else {
        return &m_mesh->property(seamsHandle);
    }
}

double PlushPatternGenerator::getSumInnerAngle(const TriMesh *mesh, HalfedgeHandle heh1, HalfedgeHandle heh2) {
    assert(mesh->to_vertex_handle(heh1) == mesh->from_vertex_handle(heh2));

    if (mesh->is_boundary(heh1) || mesh->is_boundary(heh2)) {
        HalfedgeHandle tmp = heh1;
        heh1 = mesh->opposite_halfedge_handle(heh2);
        heh2 = mesh->opposite_halfedge_handle(tmp);
    }

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

double PlushPatternGenerator::getSumInnerAngle(const TriMesh *mesh, VertexHandle v) {
    double sumAngle = 0;
    for (HalfedgeHandle cvih : mesh->vih_range(v)) {
        if (!mesh->is_boundary(cvih)) {
            sumAngle += mesh->calc_sector_angle(cvih);
        }
    }
    return sumAngle;
}

VertexHandle PlushPatternGenerator::get_original_handle(TriMesh *mesh, const VertexHandle vh) const {
    if (mesh == m_mesh) {
        return vh;
    } else {
        OpenMesh::VPropHandleT<VertexHandle> inverse_mapping = getInverseMappingHandle(mesh);
        return mesh->property(inverse_mapping, vh);
    }
}

EdgeHandle PlushPatternGenerator::get_original_handle(TriMesh *mesh, const EdgeHandle eh) const {
    if (mesh == m_mesh) {
        return eh;
    } else {
        OpenMesh::VPropHandleT<VertexHandle> inverse_mapping = getInverseMappingHandle(mesh);
        HalfedgeHandle heh = mesh->halfedge_handle(eh, 0);
        VertexHandle v1 = mesh->property(inverse_mapping, mesh->from_vertex_handle(heh));
        VertexHandle v2 = mesh->property(inverse_mapping, mesh->to_vertex_handle(heh));
        EdgeHandle original_eh;
        bool found = getEdge(m_mesh, original_eh, v1, v2);
        assert(found);
        return original_eh;
    }
}

HalfedgeHandle PlushPatternGenerator::get_original_handle(TriMesh *mesh, const HalfedgeHandle heh) const {
    if (mesh == m_mesh) {
        return heh;
    } else {
        OpenMesh::VPropHandleT<VertexHandle> inverse_mapping = getInverseMappingHandle(mesh);
        VertexHandle v1 = mesh->property(inverse_mapping, mesh->from_vertex_handle(heh));
        VertexHandle v2 = mesh->property(inverse_mapping, mesh->to_vertex_handle(heh));
        HalfedgeHandle original_heh;
        bool found = getHalfedge(m_mesh, original_heh, v1, v2);
        assert(found);
        return original_heh;
    }
}

FaceHandle PlushPatternGenerator::get_original_handle(TriMesh *mesh, const FaceHandle fh) const {
    if (mesh == m_mesh) {
        return fh;
    } else {
        OpenMesh::VPropHandleT<VertexHandle> inverse_mapping = getInverseMappingHandle(mesh);
        TriMesh::ConstFaceVertexIter cfv_it = mesh->cfv_iter(fh);
        VertexHandle v1 = mesh->property(inverse_mapping, *cfv_it++);
        VertexHandle v2 = mesh->property(inverse_mapping, *cfv_it++);
        VertexHandle v3 = mesh->property(inverse_mapping, *cfv_it++);
        FaceHandle original_fh;
        bool found = getFace(m_mesh, original_fh, v1, v2, v3);
        assert(found);
        return original_fh;
    }
}

double PlushPatternGenerator::calcArea(OpenMesh::Vec3d p1, OpenMesh::Vec3d p2, OpenMesh::Vec3d p3) {
    OpenMesh::Vec3d v12 = p2-p1;
    OpenMesh::Vec3d v13 = p3-p1;
    return (v12 % v13).norm() / 2;
}

bool PlushPatternGenerator::is_different_texture(TriMesh *mesh, HalfedgeHandle heh) {
    return is_different_texture(mesh, mesh->edge_handle(heh));
}

bool PlushPatternGenerator::is_different_texture(TriMesh *mesh, EdgeHandle eh) {
    HalfedgeHandle heh = mesh->halfedge_handle(eh, 0);
    assert(mesh->has_face_colors());
    if (mesh->is_boundary(eh)
    ||  mesh->color(mesh->face_handle(heh)) == mesh->color(mesh->opposite_face_handle(heh))) {
        return false;
    } else {
        return true;
    }
}

/** Expand selection by n-ring connectivity **/
void PlushPatternGenerator::expandVertice(const TriMesh *mesh, const VertexHandle centerV, std::set<VertexHandle> &verticesSelection, int n, double maxDistance) {
    TriMesh::Point p = mesh->point(centerV);
    for (int i = 0; i < n; i++) {
        std::set<VertexHandle> ring;
        for (VertexHandle v : verticesSelection) {
            for (VertexHandle cvv : mesh->vv_range(v)) {
                if ((p - mesh->point(cvv)).norm() < maxDistance) {
                    ring.insert(cvv);
                }
            }
        }
        verticesSelection.insert(ring.begin(), ring.end());
    }
}

void PlushPatternGenerator::expandVertice(const TriMesh *mesh, std::set<VertexHandle> &verticesSelection, int n) {
    expandVertice(mesh, *mesh->vertices_begin(), verticesSelection, n, std::numeric_limits<double>::max());
}

void PlushPatternGenerator::shrinkVertice(const TriMesh *mesh, std::set<VertexHandle> &verticesSelection, int n) {
    for (int i = 0; i < n; i++) {
        std::set<VertexHandle> innerVertices;
        for (VertexHandle v : verticesSelection) {
            bool isInner = true;
            for (VertexHandle cvv : mesh->vv_range(v)) {
                if (verticesSelection.find(cvv) == verticesSelection.end()) {
                    isInner = false;
                    break;
                }
            }
            if (isInner) {
                innerVertices.insert(v);
            }
        }
        verticesSelection = innerVertices;
    }
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

        HalfedgeHandle start_heh = *he_it;
        HalfedgeHandle current_heh = start_heh;
        do {
            visited.insert(current_heh);
            boundary.push_back(current_heh);
            current_heh = mesh->next_halfedge_handle(current_heh);
        } while (current_heh != start_heh);
        
        assert(boundary.size() > 0);
        assert(mesh->from_vertex_handle(start_heh) == mesh->to_vertex_handle(boundary[boundary.size()-1]));
        
        // The loop of opposite halfedges is in reversed direction
        if (getInteriorHalfedge) {
            std::reverse(boundary.begin(), boundary.end());
            for (size_t i = 0; i < boundary.size(); i++) {
                boundary[i] = mesh->opposite_halfedge_handle(boundary[i]);
                assert(!mesh->is_boundary(boundary[i]));
            }
        }
        boundaries.push_back(boundary);
    }
    return boundaries.empty();
}

/**
 * @brief This function returns a set of rings of halfedges. These rings are separated by given seams.
 * @param closed_seams Result loops will be stored here.
 * @param seams given as separatror.
 */
void PlushPatternGenerator::get_closed_boundaries_of_seams(std::vector< std::vector<HalfedgeHandle> > *closed_seams,
                                                           std::set<EdgeHandle> *seams) {
    std::set<HalfedgeHandle> isVisited;
    for (EdgeHandle eh : *seams) {
        HalfedgeHandle heh = m_mesh->halfedge_handle(eh, 0);
        VertexHandle startV = m_mesh->from_vertex_handle(heh);
        
        for (TriMesh::VertexOHalfedgeIter voh_it = m_mesh->voh_iter(startV); voh_it; voh_it++) {
            // Skip visited edges and non-spanning-tree edges
            if (isVisited.find(*voh_it) != isVisited.end()
            ||  seams->find(m_mesh->edge_handle(*voh_it)) == seams->end()) {
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
                } while (std::find(seams->begin(), seams->end(), m_mesh->edge_handle(next_heh)) == seams->end());
                
                // Flip it to point to the next vertex
                current_heh = m_mesh->opposite_halfedge_handle(next_heh);
                isVisited.insert(current_heh);
                
                boundary.push_back(current_heh);
            } while (current_heh != start_heh);
            
            // Delete last duplicated halfedge
            boundary.pop_back();
            
            closed_seams->push_back(boundary);
        }
    }
}

bool PlushPatternGenerator::extract_mesh_with_boundary(TriMesh *new_mesh, FaceHandle root_face, std::vector<HalfedgeHandle> *global_seams)
{
    std::set<EdgeHandle> eh_global_seams;
    for (HalfedgeHandle heh : *global_seams) {
        eh_global_seams.insert(m_mesh->edge_handle(heh));
    }
    extract_mesh_with_boundary(new_mesh, root_face, &eh_global_seams);
}
/**
 Create a new sub mesh from m_mesh. The new sub mesh is a connected component containing root_face and bounded by given boundary.
 @param new_mesh Resulting sub meshes will be stored here.
 @param root_face The given face handle used for identify which sub mesh should be extract.
 @param global_seam The given boundary.
 */
bool PlushPatternGenerator::extract_mesh_with_boundary(TriMesh *new_mesh, FaceHandle root_face, std::set<EdgeHandle> *global_seams)
{
    assert(!global_seams->empty() && "seams should not be empty.");
    
    OpenMesh::VPropHandleT<VertexHandle> inverse_mapping = getInverseMappingHandle(new_mesh);
    
    // For non-boundary vertices, we construct a one-to-one mapping from old to new
    std::map<VertexHandle, VertexHandle> vertex_to_vertex_mapping;
    // But for boundary vertices, one old vertex maps to two new vertices. This is the problem.
    // So we use halfedge-vertex mapping because two old halfedges can map to two new vertices.
    std::map<HalfedgeHandle, VertexHandle> halfedge_to_vertex_mapping;
    
    std::set<VertexHandle> seams_vertices;
    for (EdgeHandle eh : *global_seams) {
        HalfedgeHandle heh = m_mesh->halfedge_handle(eh, 0);
        seams_vertices.insert(m_mesh->from_vertex_handle(heh));
        seams_vertices.insert(m_mesh->to_vertex_handle(heh));
    }

    std::queue<FaceHandle> queue;
    
    std::set<FaceHandle> visitedF;
    
    queue.push(root_face);
    visitedF.insert(root_face);
    
    // Explore & insert new vertices
    while (!queue.empty()) {
        FaceHandle f = queue.front();
        queue.pop();
        
        TriMesh::ConstFaceHalfedgeIter cfh_it = m_mesh->cfh_iter(f);
        HalfedgeHandle heh[3];
        heh[0] = *cfh_it++;
        heh[1] = *cfh_it++;
        heh[2] = *cfh_it++;
        
        for (int i = 0; i < 3; i++) {
            VertexHandle vi = m_mesh->from_vertex_handle(heh[i]);
            // Insert vertices into new mesh
            // If this is a boundary vertex...
            if (seams_vertices.find(vi) != seams_vertices.end()) {
                // And only add new vertex if this is a boundary halfedge
                if (global_seams->find(new_mesh->edge_handle(heh[i])) != global_seams->end()) {
                    VertexHandle new_v = new_mesh->add_vertex(m_mesh->point(vi));
                    new_mesh->property(inverse_mapping, new_v) = vi;
                    halfedge_to_vertex_mapping.emplace(heh[i], new_v);
                }
            }
            // If not a boundary vertex, we have to check if this vertex is already added
            else if (vertex_to_vertex_mapping.find(vi) == vertex_to_vertex_mapping.end()) {
                VertexHandle new_v = new_mesh->add_vertex(m_mesh->point(vi));
                new_mesh->property(inverse_mapping, new_v) = vi;
                vertex_to_vertex_mapping.emplace(vi, new_v);
            }
            
            // Explore other faces through non-boundary halfedges
            FaceHandle neighborFi = m_mesh->opposite_face_handle(heh[i]);
            if (global_seams->find(new_mesh->edge_handle(heh[i])) == global_seams->end()
            &&  visitedF.find(neighborFi) == visitedF.end()) {
                queue.push(neighborFi);
                visitedF.insert(neighborFi);
            }
        }
    }
    
    // Insert new faces
    new_mesh->request_face_colors();
    for (FaceHandle f : visitedF) {
        TriMesh::ConstFaceHalfedgeIter cfh_it = m_mesh->cfh_iter(f);
        HalfedgeHandle heh[3];
        heh[0] = *cfh_it++;
        heh[1] = *cfh_it++;
        heh[2] = *cfh_it++;
        
        // Search for corresponding new vertex
        VertexHandle new_v[3];
        for (int i = 0; i < 3; i++) {
            VertexHandle vi = m_mesh->from_vertex_handle(heh[i]);
            // Non-boundary vertex
            if (seams_vertices.find(vi) == seams_vertices.end()) {
                new_v[i] = vertex_to_vertex_mapping[vi];
            }
            // boundary vertex on boundary halfedge
            else if (global_seams->find(new_mesh->edge_handle(heh[i])) != global_seams->end()) {
                new_v[i] = halfedge_to_vertex_mapping[heh[i]];
            }
            // boundary vertex on non-boundary halfedge, we need to search for neighbor boundary halfedge
            else {
                // loop until we find a boundary halfedge
                HalfedgeHandle current_heh = m_mesh->next_halfedge_handle(m_mesh->opposite_halfedge_handle(heh[i]));
                for (unsigned int count = 0; global_seams->find(new_mesh->edge_handle(current_heh)) == global_seams->end(); count++) {
                    current_heh = m_mesh->next_halfedge_handle(m_mesh->opposite_halfedge_handle(current_heh));
                    assert(count < m_mesh->valence(vi));
                }
                new_v[i] = halfedge_to_vertex_mapping[current_heh];
            }
        }
        
        std::vector<VertexHandle> face_vhs;
        face_vhs.push_back(new_v[0]);
        face_vhs.push_back(new_v[1]);
        face_vhs.push_back(new_v[2]);
        FaceHandle newF = new_mesh->add_face(face_vhs);
        new_mesh->set_color(newF, m_mesh->color(f));
    }
    return true;
}

/**
 Split existing m_mesh into sub meshes with boundary as seams.
 @param loops The given set of boundaries. Each boundary should be a loop.
 @param subMeshes Resulting sub meshes will be stored here.
 */
bool PlushPatternGenerator::splitWithBoundary(std::vector<TriMesh> *subMeshes, std::set<EdgeHandle> *seams)
{
    // For each closed boundary, generate a set of closed halfedge
    std::vector< std::vector<HalfedgeHandle> > closed_seams;
    get_closed_boundaries_of_seams(&closed_seams, seams);

    std::set<FaceHandle> visited;
    for (std::vector<HalfedgeHandle> closed_seam : closed_seams) {
        // Check if the sub mesh bounded by this boundary is already calculated
        bool is_already_calculated = true;
        for (HalfedgeHandle heh : closed_seam) {
            if (visited.find(m_mesh->face_handle(heh)) == visited.end()) {
                is_already_calculated = false;
                break;
            }
        }
        
        if (!is_already_calculated) {
            subMeshes->push_back(TriMesh());
            TriMesh *new_mesh = &subMeshes->back();
            extract_mesh_with_boundary(new_mesh, m_mesh->face_handle(*closed_seam.begin()), seams);
            
            for (FaceHandle f : new_mesh->faces()) {
                visited.insert(get_original_handle(new_mesh, f));
            }
        }
    }
    assert(visited.size() == m_mesh->n_faces());
    return true;
}

void PlushPatternGenerator::show_intersection_points() {
    std::set<VertexHandle> intersection_points;
    std::set<EdgeHandle> *seams = getSeams();
    get_intersection_points(seams, intersection_points);
    for (VertexHandle v : intersection_points) {
        m_mesh->status(v).set_feature(true);
    }
}
