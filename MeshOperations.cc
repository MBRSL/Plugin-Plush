#include "PlushPatternGenerator.hh"
#include "FilteredTriMesh.hh"

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/copy.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/dynamic_bitset.hpp>

#include <ACG/Utils/ColorCoder.hh>

#include <queue>
#include <unordered_set>

double PlushPatternGenerator::getSumInnerAngle(const TriMesh *mesh, VertexHandle v) {
    double sumAngle = 0;
    for (HalfedgeHandle cvih : mesh->vih_range(v)) {
        if (!mesh->is_boundary(cvih)) {
            sumAngle += mesh->calc_sector_angle(cvih);
        }
    }
    return sumAngle;
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
    
    if (!mesh->has_face_colors()) {
        return false;
    }
    
    HalfedgeHandle heh = mesh->halfedge_handle(eh, 0);
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

std::set<VertexHandle> PlushPatternGenerator::get_intersection_points(std::set<EdgeHandle> &seams) {
    std::set<VertexHandle> intersection_points;
    
    // Find intersection points with valance != 2
    for (EdgeHandle eh : seams) {
        HalfedgeHandle heh = m_mesh->halfedge_handle(eh, 0);
        VertexHandle v[2];
        v[0] = m_mesh->from_vertex_handle(heh);
        v[1] = m_mesh->to_vertex_handle(heh);
        
        for (int i = 0; i < 2; i++) {
            int valance = 0;
            for (HalfedgeHandle heh : m_mesh->voh_range(v[i])) {
                if (seams.find(m_mesh->edge_handle(heh)) != seams.end()) {
                    valance++;
                }
            }
            if (valance != 2) {
                intersection_points.insert(v[i]);
            }
        }
    }
    return intersection_points;
}

std::vector< std::vector<EdgeHandle> >
PlushPatternGenerator::get_segments_from_seams(std::set<EdgeHandle> &seams) {
    std::vector< std::vector< HalfedgeHandle> > heh_segments = get_halfedge_segments_from_seams(seams);

    std::vector< std::vector<EdgeHandle> > segments;
    for (auto heh_segment : heh_segments) {
        std::vector<EdgeHandle> segment;
        for (HalfedgeHandle heh : heh_segment) {
            segment.push_back(m_mesh->edge_handle(heh));
        }
        segments.push_back(segment);
    }
    return segments;
}

std::vector< std::vector<HalfedgeHandle> >
PlushPatternGenerator::get_halfedge_segments_from_seams(std::set<EdgeHandle> &seams) {
    std::set<VertexHandle> intersection_points = get_intersection_points(seams);
    
    std::vector< std::vector<HalfedgeHandle> > segments = get_halfedge_segments_from_seams(intersection_points, seams);
    
    return segments;
}
/**
 Extract oriented path segments from seams. segments are divided by seam points with valance != 2.
 The unused (isolated) intersection_points will be removed.
 @param <#parameter#>
 @return <#retval#>
 @retval <#meaning#>
 */
std::vector< std::vector<HalfedgeHandle> >
PlushPatternGenerator::get_halfedge_segments_from_seams(std::set<VertexHandle> intersection_points,
                                                        std::set<EdgeHandle> &seams) {
    std::vector< std::vector<HalfedgeHandle> > segments;
    
    // Use DFS to find segments
    auto hasher = [](EdgeHandle eh) -> int {
        return eh.idx();
    };
    std::unordered_set<EdgeHandle, decltype(hasher)> visited(0, hasher);
    
    
    for (auto start_v_it = intersection_points.begin(); start_v_it != intersection_points.end(); start_v_it++) {
        // Loop until no new segment for this intersection point
        std::vector<HalfedgeHandle> segment;
        do {
            segment.clear();
            VertexHandle v = *start_v_it;
            bool is_isolated_point = true;
            for (const HalfedgeHandle cvoh : m_mesh->voh_range(v)) {
                EdgeHandle cve = m_mesh->edge_handle(cvoh);
                if (seams.find(cve) != seams.end()) {
                    is_isolated_point = false;
                    break;
                }
            }
            // Remove isolated points
            if (is_isolated_point) {
                start_v_it = intersection_points.erase(start_v_it);
                continue;
            } else {
                do {
                    for (const HalfedgeHandle cvoh : m_mesh->voh_range(v)) {
                        EdgeHandle cve = m_mesh->edge_handle(cvoh);
                        if (seams.find(cve) != seams.end()
                            &&  visited.find(cve) == visited.end()) {
                            visited.insert(cve);
                            
                            segment.push_back(cvoh);
                            
                            VertexHandle neighbor_v = m_mesh->to_vertex_handle(cvoh);
                            v = neighbor_v;
                            is_isolated_point = false;
                            break;
                        }
                    }
                } while (intersection_points.find(v) == intersection_points.end());
            }
            // When reaching intersection point, create a new segment
            if (!segment.empty()) {
                segments.push_back(segment);
            }
        } while (!segment.empty());
    }
    
    // For the rest of the seams, they are individual loops.
    for (EdgeHandle eh : seams) {
        std::vector<HalfedgeHandle> segment;
        
        HalfedgeHandle heh = m_mesh->halfedge_handle(eh, 0);
        VertexHandle v = m_mesh->from_vertex_handle(heh);
        
        bool is_all_neighbors_visited = true;
        do {
            is_all_neighbors_visited = true;
            for (const HalfedgeHandle cvoh : m_mesh->voh_range(v)) {
                EdgeHandle cve = m_mesh->edge_handle(cvoh);
                VertexHandle neighbor_v = m_mesh->to_vertex_handle(cvoh);
                if (seams.find(cve) != seams.end()
                    &&  visited.find(cve) == visited.end()) {
                    visited.insert(cve);
                    
                    segment.push_back(cvoh);
                    
                    v = neighbor_v;
                    is_all_neighbors_visited = false;
                    break;
                }
            }
        } while (!is_all_neighbors_visited);

        if (!segment.empty()) {
            segments.push_back(segment);
        }
    }
    
    // Check if result is valid
    assert(visited.size() == seams.size());
    
#ifndef NDEBUG
    std::set<EdgeHandle> duplication_test;
    for (auto segment : segments) {
        for (HalfedgeHandle heh : segment) {
            EdgeHandle eh = m_mesh->edge_handle(heh);
            if (duplication_test.find(eh) == duplication_test.end()) {
                duplication_test.insert(eh);
            } else {
                assert("Duplicated edge handle");
            }
        }
    }
#endif
    
    return segments;
}

/**
 Find the closed boundary containing root_heh.
 @param root_heh The halfedge for identifying which closed boundary we need.
 @param seams The set containing desired boundary.
 @param boundary The resulting boundary will be stored here.
 @return <#retval#>
 @retval <#meaning#>
 */
bool PlushPatternGenerator::get_adjacent_boundary(HalfedgeHandle root_heh,
                                                  std::set<EdgeHandle> &seams,
                                                  std::set<HalfedgeHandle> &boundary) {
    auto hasher = [](FaceHandle f) -> int {
        return f.idx();
    };
    std::unordered_set<FaceHandle, decltype(hasher)> visited(0, hasher);
    std::queue<FaceHandle> queue;
    FaceHandle ff = m_mesh->face_handle(root_heh);
    queue.push(m_mesh->face_handle(root_heh));
    
    boundary.clear();
    while (!queue.empty()) {
        FaceHandle current_f = queue.front();
        queue.pop();
        visited.insert(current_f);
        
        for (const HalfedgeHandle cfh : m_mesh->fh_range(current_f)) {
            EdgeHandle eh = m_mesh->edge_handle(cfh);
            FaceHandle neighbor_f = m_mesh->opposite_face_handle(cfh);
            if (seams.find(eh) != seams.end()) {
                boundary.insert(cfh);
            } else if (neighbor_f.is_valid()
                       && visited.find(neighbor_f) == visited.end()) {
                visited.insert(neighbor_f);
                queue.push(neighbor_f);
            }
        }
    }
    
    return boundary.empty();
}

bool PlushPatternGenerator::get_joint_boundary(std::set<HalfedgeHandle> &boundary1,
                                               std::set<HalfedgeHandle> &boundary2,
                                               std::vector<HalfedgeHandle> &joint_seam_segment,
                                               std::set<HalfedgeHandle> &joint_boundary,
                                               std::set<EdgeHandle> &seams) {
    // Eliminate joint edges temporarily
    for (HalfedgeHandle heh : joint_seam_segment) {
        auto found_it = seams.find(m_mesh->edge_handle(heh));
        assert(found_it != seams.end());
        seams.erase(found_it);
    }
    
    bool result = get_adjacent_boundary(*joint_seam_segment.begin(), seams, joint_boundary);
    
    // Then restore them back after calculation
    for (HalfedgeHandle heh : joint_seam_segment) {
        seams.insert(m_mesh->edge_handle(heh));
    }
    return result;
}

/**
 @brief Return rings of boundary halfedges. May contains multiple rings
 @param getInteriorHalfedge If true, return the opposite halfedge of boundary halfedge. These halfedge are not boundary halfedge.
 */
std::vector< std::vector<HalfedgeHandle> >  PlushPatternGenerator::getBoundaryOfOpenedMesh(const FilteredTriMesh &mesh,
                                                    bool getInteriorHalfedge) {
    std::vector< std::vector<HalfedgeHandle> > boundaries;
    
    auto hasher = [](HalfedgeHandle heh) -> int {
        return heh.idx();
    };
    std::unordered_set<HalfedgeHandle, decltype(hasher)> visited(0, hasher);
    for (EdgeHandle eh : mesh.boundary_edges()) {
        for (int i = 0; i < 2; i++) {
            HalfedgeHandle heh = mesh.halfedge_handle(eh, i);
            if (!mesh.is_boundary(heh)
            ||  visited.find(heh) != visited.end()) {
                continue;
            }
            
            std::vector<HalfedgeHandle> boundary;
            
            HalfedgeHandle start_heh = heh;
            HalfedgeHandle current_heh = start_heh;
            do {
                visited.insert(current_heh);
                boundary.push_back(current_heh);
                current_heh = mesh.next_boundary_halfedge_handle(current_heh);
            } while (current_heh != start_heh);
            
            assert(boundary.size() > 0);
            assert(mesh.from_vertex_handle(start_heh) == mesh.to_vertex_handle(boundary[boundary.size()-1]));
            
            // The loop of opposite halfedges is in reversed direction
            if (getInteriorHalfedge) {
                std::reverse(boundary.begin(), boundary.end());
                for (size_t i = 0; i < boundary.size(); i++) {
                    boundary[i] = mesh.opposite_halfedge_handle(boundary[i]);
                    assert(!mesh.is_boundary(boundary[i]) || mesh.is_dual_boundary(boundary[i]));
                }
            }
            boundaries.push_back(boundary);
        }
    }
    return boundaries;
}

/**
 * @brief This function returns a set of rings of halfedges. These rings are separated by given seams.
 * @param closed_seams Result loops will be stored here.
 * @param seams given as separatror.
 */
void PlushPatternGenerator::get_closed_boundaries_of_seams(std::vector< std::vector<HalfedgeHandle> > &closed_seams,
                                                           std::set<EdgeHandle> &seams) {
    std::set<HalfedgeHandle> isVisited;
    for (EdgeHandle eh : seams) {
        HalfedgeHandle heh = m_mesh->halfedge_handle(eh, 0);
        VertexHandle startV = m_mesh->from_vertex_handle(heh);
        
        for (TriMesh::VertexOHalfedgeIter voh_it = m_mesh->voh_iter(startV); voh_it; voh_it++) {
            // Skip visited edges and non-spanning-tree edges
            if (isVisited.find(*voh_it) != isVisited.end()
            ||  seams.find(m_mesh->edge_handle(*voh_it)) == seams.end()
            ||  m_mesh->is_boundary(*voh_it)) {
                continue;
            }
            
            // We got a non-visited edge, lets start searching from this edge
            HalfedgeHandle start_heh = *voh_it;
            HalfedgeHandle current_heh = start_heh;
            isVisited.insert(start_heh);
            
            std::vector<HalfedgeHandle> boundary;
            boundary.push_back(current_heh);
            do {
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
                } while (std::find(seams.begin(), seams.end(), m_mesh->edge_handle(next_heh)) == seams.end());
                
                // Flip it to point to the next vertex
                current_heh = m_mesh->opposite_halfedge_handle(next_heh);
                isVisited.insert(current_heh);
                
                boundary.push_back(current_heh);
            } while (current_heh != start_heh);
            
            // Delete last duplicated halfedge
            boundary.pop_back();
            
            closed_seams.push_back(boundary);
        }
    }
}

FilteredTriMesh PlushPatternGenerator::get_subMesh_with_boundary(FaceHandle root_face, std::set<HalfedgeHandle> &seams)
{
    std::set<EdgeHandle> eh_seams;
    for (HalfedgeHandle heh : seams) {
        eh_seams.insert(m_mesh->edge_handle(heh));
    }
    return get_subMesh_with_boundary(root_face, eh_seams);
}

/**
 Create a new sub mesh from m_mesh. The new sub mesh is a connected component containing root_face and bounded by given boundary.
 This function modifies the tagged() property of m_mesh.
 @param new_mesh Resulting sub meshes will be stored here.
 @param root_face The given face handle used for identify which sub mesh should be extract.
 @param seams The given boundary of original mesh.
 */
FilteredTriMesh PlushPatternGenerator::get_subMesh_with_boundary(FaceHandle root_face, std::set<EdgeHandle> &seams)
{
    assert(!seams.empty() && "seams should not be empty.");
    
    auto hasherF = [](FaceHandle fh) -> int {
        return fh.idx();
    };
    std::unordered_set<FaceHandle, decltype(hasherF)> visitedF(0, hasherF);
    
    std::queue<FaceHandle> queue;
    queue.push(root_face);
    visitedF.insert(root_face);
    
    // Explore
    std::set<FaceHandle> faces;
    std::set<EdgeHandle> boundary_edges;
    while (!queue.empty()) {
        FaceHandle f = queue.front();
        queue.pop();
        
        faces.insert(f);
        
        for (const HalfedgeHandle cfh : m_mesh->fh_range(f)) {
            // Explore other faces through non-boundary halfedges
            FaceHandle neighborFi = m_mesh->opposite_face_handle(cfh);
            if (seams.find(m_mesh->edge_handle(cfh)) == seams.end()) {
                if (visitedF.find(neighborFi) == visitedF.end()) {
                    queue.push(neighborFi);
                    visitedF.insert(neighborFi);
                }
            } else {
                boundary_edges.insert(m_mesh->edge_handle(cfh));
            }
        }
    }
    
    FilteredTriMesh filtered_mesh(m_mesh, faces, boundary_edges);
    return filtered_mesh;
}

/**
 Split existing m_mesh into sub meshes with boundary as seams.
 @param subMeshes
    Resulting sub meshes will be stored here.
 @param subMesh_graph
    This variable store the relationship between each sub-mesh as a graph.
    The vertices in this graph have a property of pointer pointing to the returned sub-mesh.
    If the returned variable is destroyed, the pointers in this graph also become invalid.
 */
PlushPatternGenerator::SubMesh_graph PlushPatternGenerator::get_subMeshes_with_boundary(std::set<EdgeHandle> &seams)
{
    SubMesh_graph subMesh_graph;
    
    // Init faces to invalid sub-mesh id.
    for (FaceHandle f : m_mesh->faces()) {
        m_mesh->property(face_to_patch_idx_handle, f) = -1;
    }
    
    // Use seam segments to generate sub-meshes
    size_t edge_counter = 0;
    std::vector< std::vector<HalfedgeHandle> > segments = get_halfedge_segments_from_seams(seams);
    for (std::vector<HalfedgeHandle> &segment : segments) {
        // Choose a non-clustered face as the root to extract sub-mesh
        FaceHandle root_f[2];
        root_f[0] = m_mesh->face_handle(segment.front());
        root_f[1] = m_mesh->opposite_face_handle(segment.front());
        
        SubMesh_graph::vertex_descriptor subMesh_idx[2];
        for (int i = 0; i < 2; i++) {
            assert(root_f[i].is_valid());
            if (m_mesh->property(face_to_patch_idx_handle, root_f[i]) < 0) {
                // Add this subMesh to SubMesh_graph as a vertex
                subMesh_idx[i] = boost::add_vertex(get_subMesh_with_boundary(root_f[i], seams), subMesh_graph);
                FilteredTriMesh &subMesh = boost::get(boost::vertex_owner, subMesh_graph, subMesh_idx[i]);
                
                subMesh.merged_subMesh_idx.resize(segments.size()*2);
                subMesh.merged_subMesh_idx[subMesh_idx[i]] = 1;
                
                // Assign sub-mesh id to all faces of this sub-mesh
                for (FaceHandle f : subMesh.faces()) {
                    m_mesh->property(face_to_patch_idx_handle, f) = subMesh_idx[i];
                }
            } else {
                subMesh_idx[i] = m_mesh->property(face_to_patch_idx_handle, root_f[i]);
            }
        }
        if (subMesh_idx[0] != subMesh_idx[1]) {
            auto result = boost::add_edge(subMesh_idx[0], subMesh_idx[1], subMesh_graph);
            bool success = result.second;
            assert(success);
            
            SubMesh_graph::edge_descriptor e = result.first;
            boost::put(boost::edge_index, subMesh_graph, e, edge_counter++);
            boost::put(boost::edge_owner, subMesh_graph, e, segment);
        }
    }
    
    // Re-adjust the length of ids
    int n_vertices = boost::num_vertices(subMesh_graph);
    int n_edges = boost::num_edges(subMesh_graph);
    SubMesh_graph::vertex_iterator v_it, v_ite;
    for (boost::tie(v_it, v_ite) = boost::vertices(subMesh_graph); v_it != v_ite; v_it++) {
        FilteredTriMesh &filtered_mesh = boost::get(boost::vertex_owner, subMesh_graph, *v_it);
        filtered_mesh.merged_subMesh_idx.resize(n_vertices);
        filtered_mesh.merged_seam_idx.resize(n_edges);
    }
    
    return subMesh_graph;
}

// I write this function converting Patch_boundary to FilteredTriMesh
bool PlushPatternGenerator::convert_patch_boundary_to_FilteredTriMesh(FilteredTriMesh &subMesh, Patch_boundary &patch_boundary) {
    std::set<FaceHandle> faces;
    SubMesh_graph &subMesh_graph = m_mesh->property(subMeshes_handle);
    for (SubMesh_graph::vertex_descriptor v1 = patch_boundary.merged_subMesh_idx.find_first();
         v1 != patch_boundary.merged_subMesh_idx.npos;
         v1 = patch_boundary.merged_subMesh_idx.find_next(v1)) {
        FilteredTriMesh &mesh = boost::get(boost::vertex_owner, subMesh_graph, v1);
        for (FaceHandle fh : mesh.faces()) {
            faces.insert(fh);
        }
    }
    
    std::set<EdgeHandle> boundary_edges;
    SubMesh_graph::edge_iterator e_it, e_ite;
    for (boost::tie(e_it, e_ite) = boost::edges(subMesh_graph); e_it != e_ite; e_it++) {
        auto idx = boost::get(boost::edge_index, subMesh_graph, *e_it);
        if (patch_boundary.boundary_seam_idx[idx]) {
            auto seam_segment = boost::get(boost::edge_owner, subMesh_graph, *e_it);
            for (HalfedgeHandle heh : seam_segment) {
                boundary_edges.insert(m_mesh->edge_handle(heh));
            }
        }
    }
    subMesh = FilteredTriMesh(m_mesh, faces, boundary_edges);
    
    return true;
}

/*
 @brief Create a TriMesh corresponds to the input FilteredTriMesh respect to its boundary.
 Due to the mechanism of OpenFlipper, one should only use the TriMesh create by OpenFlipper.
 That is, we can not "new" our TriMesh by ourself. It should be create first then passed as parameter
 @param result_triMesh The result triMesh will be stored here. One should pass the TriMesh pointer created by OpenFlipper.
 @param subMesh The input FilteredTriMesh containing sub-mesh & boundary information.
 @param use_flattened_position (Debug) Use original position or flattened position calculated by parameterization.
 */
bool PlushPatternGenerator::get_triMesh_from_subMesh(TriMesh *result_triMesh, FilteredTriMesh &subMesh, bool use_flattened_position) {
    OpenMesh::VPropHandleT<VertexHandle> inverse_mapping = getInverseMappingHandle(result_triMesh);

    // We use the idx from subMesh to map old vertice to new vertice
    std::map<int, VertexHandle> vertex_mapping;
    
    // Insert new vertices
    for (EdgeHandle eh : subMesh.boundary_edges()) {
        for (int i = 0; i < 2; i++) {
            HalfedgeHandle heh = subMesh.halfedge_handle(eh, i);
            if (subMesh.is_boundary(heh)) {
                VertexHandle vh = subMesh.from_vertex_handle(heh);
                
                TriMesh::Point p = subMesh.flattened_boundary_point(heh);
                if (!use_flattened_position) {
                    p = subMesh.point(vh);
                }
                VertexHandle new_vh = result_triMesh->add_vertex(p);
                
                result_triMesh->property(inverse_mapping, new_vh) = vh;
                int v_id = subMesh.get_boundary_vertex_idx(heh);
                vertex_mapping.emplace(v_id, new_vh);
            }
        }
    }
    for (VertexHandle vh : subMesh.vertices()) {
        if (!subMesh.is_boundary(vh)) {
            TriMesh::Point p = subMesh.flattened_non_boundary_point(vh);
            if (!use_flattened_position) {
                p = subMesh.point(vh);
            }
            VertexHandle new_vh = result_triMesh->add_vertex(p);
            result_triMesh->property(inverse_mapping, new_vh) = vh;
            int v_id = subMesh.get_non_boundary_vertex_idx(vh);
            vertex_mapping.emplace(v_id, new_vh);
        }
    }

    // Now new vertices are ready, we can use them to insert new faces
    for (FaceHandle fh : subMesh.faces()) {
        std::vector<VertexHandle> face_vhs;
        // Search for corresponding new vertex
        for (const HalfedgeHandle cfh : subMesh.fh_range(fh)) {
            VertexHandle vi = subMesh.from_vertex_handle(cfh);
            int v_id = -1;
            if (subMesh.is_boundary(vi)) {
                HalfedgeHandle heh = cfh;
                v_id = subMesh.get_boundary_vertex_idx(heh);
            } else {
                v_id = subMesh.get_non_boundary_vertex_idx(vi);
            }
            assert(vertex_mapping[v_id].is_valid());
            face_vhs.push_back(vertex_mapping[v_id]);
        }
        assert(face_vhs.size() == 3);
        
        FaceHandle newF = result_triMesh->add_face(face_vhs);
        result_triMesh->set_color(newF, TriMesh::Color(1, 1, 1, 1));
    }

    // Assigning edge color for boundary
    result_triMesh->request_edge_colors();
    ACG::ColorCoder color_coder(0, 255, false);
    for (EdgeHandle eh : result_triMesh->edges()) {
        if (result_triMesh->is_boundary(eh)) {
            EdgeHandle original_eh = get_original_handle(result_triMesh, eh);
            int segment_no = m_mesh->property(segment_no_handle, original_eh);
            
            TriMesh::Color color = color_coder.color_float4(segment_no);
            
            result_triMesh->set_color(eh, color);
        } else {
            result_triMesh->set_color(eh, TriMesh::Color(1, 1, 1, 0));
        }
    }
    
    result_triMesh->update_normals();

    // Update seams
//    std::vector<int> subSeamsId;
//    
//    std::set<EdgeHandle> &seams = subMesh->property(seamsHandle);
//    std::set<EdgeHandle> &oldSeams = oldSubMesh.property(oldSeamsHandle);
//    for (auto e_it = oldSeams.begin(); e_it != oldSeams.end(); e_it++) {
//        HalfedgeHandle old_heh = oldSubMesh.halfedge_handle(*e_it, 0);
//        VertexHandle oldV1 = oldSubMesh.from_vertex_handle(old_heh);
//        VertexHandle oldV2 = oldSubMesh.to_vertex_handle(old_heh);
//        EdgeHandle he;
//        bool edgeExist = PlushPatternGenerator::getEdge(subMesh, he, oldToNewMapping[oldV1], oldToNewMapping[oldV2]);
//        assert(edgeExist);
//        
//        // Update inverse map
//        seams.insert(he);
//        subSeamsId.push_back(he.idx());
//        
//        VertexHandle oldOriginalV1 = oldSubMesh.property(oldInverseMappingHandle, oldV1);
//        VertexHandle oldOriginalV2 = oldSubMesh.property(oldInverseMappingHandle, oldV2);
//        EdgeHandle original_he;
//        edgeExist = PlushPatternGenerator::getEdge(m_triMeshObj->mesh(), original_he, oldOriginalV1, oldOriginalV2);
//        assert(edgeExist);
//        seamsId.push_back(original_he.idx());
//    }

    return true;
}

/*
 * v4 -------- v1
 *  \         / \
 *   \       /   \
 *    \     /     \
 *     \   /       \
 *      \ /         \
 *      v2 -------- v3
 */
Eigen::SparseMatrix<double> PlushPatternGenerator::calc_parameterization_weight_matrix(FilteredTriMesh &mesh,
                                                                                       InteriorParameterizationMethod method) {
    
    int n = mesh.n_duplicated_vertices();
    Eigen::SparseMatrix<double> M(n,n);
    
    // Prepare for matrix
    std::vector< Eigen::Triplet<double> > tripletList;
    double *rowSum = new double[n];
    memset(rowSum, 0, sizeof(double)*n);
    
    for (EdgeHandle eh : mesh.edges()) {
        if (mesh.is_boundary(eh)) {
            continue;
        }
        
        HalfedgeHandle heh12 = mesh.halfedge_handle(eh, 0);
        HalfedgeHandle heh21 = mesh.halfedge_handle(eh, 1);
        
        VertexHandle v1 = mesh.from_vertex_handle(heh12);
        VertexHandle v2 = mesh.to_vertex_handle(heh12);
        
        HalfedgeHandle heh31 = mesh.prev_halfedge_handle(heh12);
        HalfedgeHandle heh42 = mesh.prev_halfedge_handle(heh21);
        
        double weight12 = 0, weight21 = 0;
        if (method == Barycentric) {
            weight12 = 1;
            weight21 = 1;
        } else if (method == MeanValue) {
            double angle312 = mesh.calc_sector_angle(heh31);
            double angle214 = mesh.calc_sector_angle(heh21);
            double angle123 = mesh.calc_sector_angle(heh12);
            double angle421 = mesh.calc_sector_angle(heh42);
            weight12 = (tan(angle312/2) + tan(angle214/2)) / mesh.calc_edge_length(heh12);
            weight21 = (tan(angle123/2) + tan(angle421/2)) / mesh.calc_edge_length(heh21);
        } else if (method == Conformal) {
            double angle231 = mesh.calc_sector_angle(mesh.next_halfedge_handle(heh12));
            double angle142 = mesh.calc_sector_angle(mesh.next_halfedge_handle(heh21));
            double cotR231 = tan(M_PI_2 - angle231);
            double cotR142 = tan(M_PI_2 - angle142);
            weight12 = cotR231 + cotR142;
            weight21 = weight12;
        } else {
            assert("Invalid InteriorParameterizationMethod");
        }
        
        int idx1, idx2;
        if (!mesh.is_boundary(v1)) {
            idx1 = mesh.get_non_boundary_vertex_idx(v1);
        } else {
            idx1 = mesh.get_boundary_vertex_idx(heh12);
        }
        if (!mesh.is_boundary(v2)) {
            idx2 = mesh.get_non_boundary_vertex_idx(v2);
        } else {
            idx2 = mesh.get_boundary_vertex_idx(heh21);
        }
        
        assert(idx1 >= 0 && idx1 < n
            && idx2 >= 0 && idx2 < n);
        
        if (!mesh.is_boundary(v1)) {
            tripletList.push_back(Eigen::Triplet<double>(idx1,
                                                         idx2,
                                                         weight12));
            rowSum[idx1] += weight12;
        }
        if (!mesh.is_boundary(v2)) {
            tripletList.push_back(Eigen::Triplet<double>(idx2,
                                                         idx1,
                                                         weight21));
            rowSum[idx2] += weight21;
        }
    }
    
    for (VertexHandle vh : mesh.vertices()) {
        if (!mesh.is_boundary(vh)) {
            int idx = mesh.get_non_boundary_vertex_idx(vh);
            tripletList.push_back(Eigen::Triplet<double>(idx, idx, -rowSum[idx]));
        }
    }
    for (HalfedgeHandle heh : mesh.halfedges()) {
        if (mesh.is_boundary(heh)) {
            int idx = mesh.get_boundary_vertex_idx(heh);
            tripletList.push_back(Eigen::Triplet<double>(idx, idx, 1));
        }
    }
    // Filling matrix
    M.setFromTriplets(tripletList.begin(), tripletList.end());
    M.makeCompressed();
    
    return M;
}

void PlushPatternGenerator::calc_parameterization_weight_matrix( TriMesh *mesh,
                                                                 Eigen::SparseMatrix<double> &M,
                                                                 InteriorParameterizationMethod method)  {
    int n = mesh->n_vertices();
    
    // Prepare for matrix
    std::vector< Eigen::Triplet<double> > tripletList;
    double *rowSum = new double[n];
    memset(rowSum, 0, sizeof(double)*n);
    
    for (EdgeHandle eh : mesh->edges()) {
        HalfedgeHandle heh12 = mesh->halfedge_handle(eh, 0);
        HalfedgeHandle heh21 = mesh->halfedge_handle(eh, 1);
        
        VertexHandle v1 = mesh->from_vertex_handle(heh12);
        VertexHandle v2 = mesh->to_vertex_handle(heh12);
        if (mesh->is_boundary(v1) && mesh->is_boundary(v2)) {
            continue;
        }
        
        HalfedgeHandle heh31 = mesh->prev_halfedge_handle(heh12);
        HalfedgeHandle heh42 = mesh->prev_halfedge_handle(heh21);
        
        double weight12 = 0, weight21 = 0;
        if (method == Barycentric) {
            weight12 = 1;
            weight21 = 1;
        } else if (method == MeanValue) {
            double angle312 = mesh->calc_sector_angle(heh31);
            double angle214 = mesh->calc_sector_angle(heh21);
            double angle123 = mesh->calc_sector_angle(heh12);
            double angle421 = mesh->calc_sector_angle(heh42);
            weight12 = (tan(angle312/2) + tan(angle214/2)) / mesh->calc_edge_length(heh12);
            weight21 = (tan(angle123/2) + tan(angle421/2)) / mesh->calc_edge_length(heh21);
        } else if (method == Conformal) {
            double angle231 = mesh->calc_sector_angle(mesh->next_halfedge_handle(heh12));
            double angle142 = mesh->calc_sector_angle(mesh->next_halfedge_handle(heh21));
            double cotR231 = tan(M_PI_2 - angle231);
            double cotR142 = tan(M_PI_2 - angle142);
            weight12 = cotR231 + cotR142;
            weight21 = weight12;
        } else {
            assert("Invalid InteriorParameterizationMethod");
            return;
        }
        
        if (!mesh->is_boundary(v1)) {
            tripletList.push_back(Eigen::Triplet<double>(v1.idx(), v2.idx(), weight12));
            rowSum[v1.idx()] += weight12;
        }
        if (!mesh->is_boundary(v2)) {
            tripletList.push_back(Eigen::Triplet<double>(v2.idx(), v1.idx(), weight21));
            rowSum[v2.idx()] += weight21;
        }
    }
    
    for (VertexHandle v : mesh->vertices()) {
        if (mesh->is_boundary(v)) {
            tripletList.push_back(Eigen::Triplet<double>(v.idx(), v.idx(), 1));
        } else {
            tripletList.push_back(Eigen::Triplet<double>(v.idx(), v.idx(), -rowSum[v.idx()]));
        }
    }
    // Filling matrix
    M.setFromTriplets(tripletList.begin(), tripletList.end());
    M.makeCompressed();
}
