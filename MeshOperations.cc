#include "PlushPatternGenerator.hh"

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/copy.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>

#include <queue>
#include <unordered_set>

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

double PlushPatternGenerator::calcArea(OpenMesh::Vec3d p1, OpenMesh::Vec3d p2, OpenMesh::Vec3d p3) {
    OpenMesh::Vec3d v12 = p2-p1;
    OpenMesh::Vec3d v13 = p3-p1;
    return (v12 % v13).norm() / 2;
}

bool PlushPatternGenerator::is_different_texture(TriMesh *mesh, HalfedgeHandle heh) {
    return is_different_texture(mesh, mesh->edge_handle(heh));
}

bool PlushPatternGenerator::is_different_texture(TriMesh *mesh, EdgeHandle eh) {
    return false;
    
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


//bool PlushPatternGenerator::is_loop(std::set<EdgeHandle> &boundary,
//                                    std::set<EdgeHandle> &seams) {
//    HalfedgeHandle start_heh = m_mesh->halfedge_handle(*boundary.begin(), 0);
//    VertexHandle start_v = m_mesh->from_vertex_handle(start_heh);
//
//    std::set<EdgeHandle> visited;
//    std::queue<HalfedgeHandle> queue;
//    queue.push(start_heh);
//
//    do {
//        HalfedgeHandle heh = queue.front();
//        queue.pop();
//
//        visited.insert(m_mesh->edge_handle(heh));
//        VertexHandle v = m_mesh->to_vertex_handle(heh);
//        // Return true if it meet the starting vertex
//        if (v == start_v) {
//            return true;
//        }
//
//        for (HalfedgeHandle neighbor_heh : m_mesh->voh_range(v)) {
//            EdgeHandle neighbor_eh = m_mesh->edge_handle(neighbor_heh);
//            if (seams.find(neighbor_eh) != seams.end()
//            &&  visited.find(neighbor_eh) == visited.end()) {
//                visited.insert(neighbor_eh);
//                queue.push(neighbor_heh);
//            }
//        }
//    } while (!queue.empty());
//
//    return false;
//}

void PlushPatternGenerator::get_intersection_points(std::set<EdgeHandle> *seams,
                                                    std::set<VertexHandle> &intersection_points) {
    // Find intersection points with valance != 2
    for (EdgeHandle eh : *seams) {
        HalfedgeHandle heh = m_mesh->halfedge_handle(eh, 0);
        VertexHandle v[2];
        v[0] = m_mesh->from_vertex_handle(heh);
        v[1] = m_mesh->to_vertex_handle(heh);
        
        for (int i = 0; i < 2; i++) {
            int valance = 0;
            for (HalfedgeHandle heh : m_mesh->voh_range(v[i])) {
                if (seams->find(m_mesh->edge_handle(heh)) != seams->end()) {
                    valance++;
                }
            }
            if (valance != 2) {
                intersection_points.insert(v[i]);
            }
        }
    }
}

void PlushPatternGenerator::get_segments_from_seams(std::vector< std::vector<EdgeHandle> > &segments, std::set<EdgeHandle> *seams) {
    std::vector< std::vector< HalfedgeHandle> > heh_segments;
    get_segments_from_seams(heh_segments, seams);
    for (auto heh_segment : heh_segments) {
        std::vector<EdgeHandle> segment;
        for (HalfedgeHandle heh : heh_segment) {
            segment.push_back(m_mesh->edge_handle(heh));
        }
        segments.push_back(segment);
    }
}

void PlushPatternGenerator::get_segments_from_seams(std::vector< std::vector<HalfedgeHandle> > &segments, std::set<EdgeHandle> *seams) {
    std::set<VertexHandle> intersection_points;
    get_intersection_points(seams, intersection_points);
    get_segments_from_seams(segments, intersection_points, seams);
}
/**
 Extract oriented path segments from seams. segments are divided by seam points with valance != 2.
 The unused (isolated) intersection_points will be removed.
 @param <#parameter#>
 @return <#retval#>
 @retval <#meaning#>
 */
void PlushPatternGenerator::get_segments_from_seams(std::vector< std::vector<HalfedgeHandle> > &segments,
                                                    std::set<VertexHandle> intersection_points,
                                                    std::set<EdgeHandle> *seams) {
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
                if (seams->find(cve) != seams->end()) {
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
                        if (seams->find(cve) != seams->end()
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
    for (EdgeHandle eh : *seams) {
        std::vector<HalfedgeHandle> segment;
        
        HalfedgeHandle heh = m_mesh->halfedge_handle(eh, 0);
        VertexHandle v = m_mesh->from_vertex_handle(heh);
        
        bool is_all_visited = true;
        do {
            is_all_visited = true;
            for (const HalfedgeHandle cvoh : m_mesh->voh_range(v)) {
                EdgeHandle cve = m_mesh->edge_handle(cvoh);
                VertexHandle neighbor_v = m_mesh->to_vertex_handle(cvoh);
                if (seams->find(cve) != seams->end()
                    &&  visited.find(cve) == visited.end()) {
                    visited.insert(cve);
                    
                    segment.push_back(cvoh);
                    
                    v = neighbor_v;
                    is_all_visited = false;
                    break;
                }
            }
        } while (!is_all_visited);
    }
    return;
    
    // Check if result is valid
    assert(visited.size() == seams->size());
    
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
            ||  seams->find(m_mesh->edge_handle(*voh_it)) == seams->end()
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

bool PlushPatternGenerator::extract_mesh_with_boundary(TriMesh *new_mesh, FaceHandle root_face, std::set<HalfedgeHandle> *seams)
{
    std::set<EdgeHandle> eh_seams;
    for (HalfedgeHandle heh : *seams) {
        eh_seams.insert(m_mesh->edge_handle(heh));
    }
    return extract_mesh_with_boundary(new_mesh, root_face, &eh_seams);
}
/**
 Create a new sub mesh from m_mesh. The new sub mesh is a connected component containing root_face and bounded by given boundary.
 This function modifies the tagged() property of m_mesh.
 @param new_mesh Resulting sub meshes will be stored here.
 @param root_face The given face handle used for identify which sub mesh should be extract.
 @param seams The given boundary of original mesh.
 */
bool PlushPatternGenerator::extract_mesh_with_boundary(TriMesh *new_mesh, FaceHandle root_face, std::set<EdgeHandle> *seams)
{
    assert(!seams->empty() && "seams should not be empty.");
    
    OpenMesh::VPropHandleT<VertexHandle> inverse_mapping = getInverseMappingHandle(new_mesh);
    
    // For non-boundary vertices, we construct a one-to-one mapping from old to new
    std::map<VertexHandle, VertexHandle> vertex_to_vertex_mapping;
    // But for boundary vertices, one old vertex may maps to two new vertices. This is the problem.
    // So we use halfedge-vertex mapping because two old halfedges can map to two new vertices.
    std::map<HalfedgeHandle, VertexHandle> halfedge_to_vertex_mapping;
    
    auto hasherV = [](VertexHandle vh) -> int {
        return vh.idx();
    };
    std::unordered_set<VertexHandle, decltype(hasherV)> seams_vertices(0, hasherV);
    for (EdgeHandle eh : *seams) {
        HalfedgeHandle heh = m_mesh->halfedge_handle(eh, 0);
        seams_vertices.insert(m_mesh->from_vertex_handle(heh));
        seams_vertices.insert(m_mesh->to_vertex_handle(heh));
    }
    
    auto hasherF = [](FaceHandle fh) -> int {
        return fh.idx();
    };
    std::unordered_set<FaceHandle, decltype(hasherF)> visitedF(0, hasherF);
    
    std::queue<FaceHandle> queue;
    queue.push(root_face);
    visitedF.insert(root_face);
    
    // Explore & insert new vertices
    while (!queue.empty()) {
        FaceHandle f = queue.front();
        queue.pop();
        
        for (const HalfedgeHandle cfh : m_mesh->fh_range(f)) {
            VertexHandle vi = m_mesh->from_vertex_handle(cfh);
            // Insert vertices into new mesh
            // If this is a boundary vertex...
            if (seams_vertices.find(vi) != seams_vertices.end()) {
                // And only add new vertex if this is a boundary halfedge
                if (seams->find(m_mesh->edge_handle(cfh)) != seams->end()) {
                    VertexHandle new_v = new_mesh->add_vertex(m_mesh->point(vi));
                    new_mesh->property(inverse_mapping, new_v) = vi;
                    halfedge_to_vertex_mapping.emplace(cfh, new_v);
                }
            }
            // If not a boundary vertex, we have to check if this vertex is already added
            else if (vertex_to_vertex_mapping.find(vi) == vertex_to_vertex_mapping.end()) {
                VertexHandle new_v = new_mesh->add_vertex(m_mesh->point(vi));
                new_mesh->property(inverse_mapping, new_v) = vi;
                vertex_to_vertex_mapping.emplace(vi, new_v);
            }
            
            // Explore other faces through non-boundary halfedges
            FaceHandle neighborFi = m_mesh->opposite_face_handle(cfh);
            if (seams->find(m_mesh->edge_handle(cfh)) == seams->end()
                &&  visitedF.find(neighborFi) == visitedF.end()) {
                queue.push(neighborFi);
                visitedF.insert(neighborFi);
            }
        }
    }
    
    // Insert new faces
    new_mesh->request_face_colors();
    for (FaceHandle f : visitedF) {
        std::vector<VertexHandle> face_vhs;
        // Search for corresponding new vertex
        for (const HalfedgeHandle cfh : m_mesh->fh_range(f)) {
            VertexHandle vi = m_mesh->from_vertex_handle(cfh);
            // Non-boundary vertex
            if (seams_vertices.find(vi) == seams_vertices.end()) {
                face_vhs.push_back(vertex_to_vertex_mapping[vi]);
            }
            // boundary vertex on boundary halfedge
            else if (seams->find(m_mesh->edge_handle(cfh)) != seams->end()) {
                face_vhs.push_back(halfedge_to_vertex_mapping[cfh]);
            }
            // boundary vertex on non-boundary halfedge, we need to search for neighbor boundary halfedge
            else {
                // loop until we find a boundary halfedge
                HalfedgeHandle current_heh = m_mesh->next_halfedge_handle(m_mesh->opposite_halfedge_handle(cfh));
                for (unsigned int count = 0; seams->find(m_mesh->edge_handle(current_heh)) == seams->end(); count++) {
                    current_heh = m_mesh->next_halfedge_handle(m_mesh->opposite_halfedge_handle(current_heh));
                    assert(count < m_mesh->valence(vi));
                }
                face_vhs.push_back(halfedge_to_vertex_mapping[current_heh]);
            }
        }
        assert(!face_vhs.empty());

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
