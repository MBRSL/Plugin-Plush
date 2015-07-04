#include "OpenMesh_Boost_Wrapper.hh"
#include "PlushPatternGenerator.hh"
#include "WeightFunctor.hh"

#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/depth_first_search.hpp>

#include <queue>

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

void PlushPatternGenerator::get_segments_from_seams(std::vector< std::vector<EdgeHandle> > &segments) {
    std::vector< std::vector< HalfedgeHandle> > heh_segments;
    get_segments_from_seams(heh_segments);
    for (auto heh_segment : heh_segments) {
        std::vector<EdgeHandle> segment;
        for (HalfedgeHandle heh : heh_segment) {
            segment.push_back(m_mesh->edge_handle(heh));
        }
        segments.push_back(segment);
    }
}

void PlushPatternGenerator::get_segments_from_seams(std::vector< std::vector<HalfedgeHandle> > &segments) {
    std::set<EdgeHandle> *seams = getSeams();
    std::set<VertexHandle> intersection_points;
    get_intersection_points(seams, intersection_points);
    get_segments_from_seams(segments, intersection_points);
}
/**
 Extract oriented path segments from seams. segments are divided by seam points with valance != 2.
 @param <#parameter#>
 @return <#retval#>
 @retval <#meaning#>
 */
void PlushPatternGenerator::get_segments_from_seams(std::vector< std::vector<HalfedgeHandle> > &segments,
                                                    std::set<VertexHandle> intersection_points) {
    std::set<EdgeHandle> *seams = getSeams();
    
    // Use DFS to find segments
    std::set<EdgeHandle> visited;

    for (VertexHandle start_v : intersection_points) {
        // Loop until no new segment for this intersection point
        std::vector<HalfedgeHandle> segment;
        do {
            segment.clear();
            VertexHandle v = start_v;
            do {
                for (const HalfedgeHandle cvoh : m_mesh->voh_range(v)) {
                    EdgeHandle cve = m_mesh->edge_handle(cvoh);
                    VertexHandle neighbor_v = m_mesh->to_vertex_handle(cvoh);
                    if (seams->find(cve) != seams->end()
                    &&  visited.find(cve) == visited.end()) {
                        visited.insert(cve);
                        
                        segment.push_back(cvoh);

                        v = neighbor_v;
                        break;
                    }
                }
            } while (intersection_points.find(v) == intersection_points.end());

            // When reaching intersection point, create a new segment
            if (!segment.empty()) {
                segments.push_back(segment);
            }
        } while (!segment.empty());
    }
    
    // For the rest of the seams, they are indivisual loops.
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
    
    std::set<EdgeHandle> duplication_test;
    for (auto segment : segments) {
        for (HalfedgeHandle heh : segment) {
            EdgeHandle eh = m_mesh->edge_handle(heh);
            if (duplication_test.find(eh) == visited.end()) {
                duplication_test.insert(eh);
            } else {
                assert("Duplicated edge handle");
            }
        }
    }
}

/// Find the two adjacent circular boundary along segment
bool PlushPatternGenerator::get_adjacent_boundary(HalfedgeHandle root_heh,
                                                    std::set<EdgeHandle> &seams,
                                                    std::vector<HalfedgeHandle> &boundaries) {
    std::queue<HalfedgeHandle> queue;
    std::set<HalfedgeHandle> visited;
    
    queue.push(root_heh);
    visited.insert(root_heh);
    
    while (!queue.empty()) {
        HalfedgeHandle heh = queue.front();
        queue.pop();
        
        VertexHandle v = m_mesh->to_vertex_handle(heh);
        HalfedgeHandle current_heh = heh;
        do {
            current_heh = m_mesh->opposite_halfedge_handle(
                                                           m_mesh->next_halfedge_handle(current_heh));
        } while (current_heh != heh && seams.find(m_mesh->edge_handle(current_heh)) == seams.end());
        
        // It seems that this is the end of an one-way seam. Don't worry, we'll continue by reversing direction
        if (current_heh == heh) { }
        
        current_heh = m_mesh->opposite_halfedge_handle(current_heh);
        if (visited.find(current_heh) == visited.end()) {
            queue.push(current_heh);
            visited.insert(current_heh);
        }
    }
    
    if (visited.empty()) {
        return false;
    }
//    for (HalfedgeHandle heh : visited) {
//        VertexHandle v = m_mesh->from_vertex_handle(heh);
//        int valance = 0;
//        for (EdgeHandle cve : m_mesh->ve_range(v)) {
//            HalfedgeHandle cve_heh1 = m_mesh->halfedge_handle(cve, 0);
//            HalfedgeHandle cve_heh2 = m_mesh->halfedge_handle(cve, 1);
//            if (visited.find(cve_heh1) != visited.end()
//                ||  visited.find(cve_heh2) != visited.end()) {
//                valance++;
//            }
//        }
//        if (valance != 2) {
//            return false;
//        }
//    }
    
    // Forming a oriented vector of joint boundary
    HalfedgeHandle starting_heh = *visited.begin();
    HalfedgeHandle current_heh = starting_heh;
    do {
        boundaries.push_back(current_heh);
        HalfedgeHandle next_heh = current_heh;
        do {
            next_heh = m_mesh->opposite_halfedge_handle(
                                                        m_mesh->next_halfedge_handle(next_heh));
        } while (next_heh != current_heh
                 && visited.find(m_mesh->opposite_halfedge_handle(next_heh)) == visited.end());
        
        EdgeHandle eh = m_mesh->edge_handle(next_heh);
        current_heh = m_mesh->opposite_halfedge_handle(next_heh);
    } while (current_heh != starting_heh);
    
    return true;
}

bool PlushPatternGenerator::get_joint_boundary(std::vector<HalfedgeHandle> &boundaries1,
                                               std::vector<HalfedgeHandle> &boundaries2,
                                               std::set<EdgeHandle> &seams,
                                               std::vector<HalfedgeHandle> &joint_boundary) {
    std::set<HalfedgeHandle> visited;
    
    visited.insert(boundaries1.begin(), boundaries1.end());
    visited.insert(boundaries2.begin(), boundaries2.end());
    
    // Eliminate adjacent edges
    bool changed = true;
    while (changed) {
        changed = false;
        for (HalfedgeHandle heh : visited) {
            if (visited.find(m_mesh->opposite_halfedge_handle(heh)) != visited.end()) {
                visited.erase(heh);
                visited.erase(m_mesh->opposite_halfedge_handle(heh));
                changed = true;
                break;
            }
        }
    }
    
    // Check if joint boundary forms a ring
    if (visited.empty()) {
        return false;
    }
    for (HalfedgeHandle heh : visited) {
        VertexHandle v = m_mesh->from_vertex_handle(heh);
        int valance = 0;
        for (EdgeHandle cve : m_mesh->ve_range(v)) {
            HalfedgeHandle cve_heh1 = m_mesh->halfedge_handle(cve, 0);
            HalfedgeHandle cve_heh2 = m_mesh->halfedge_handle(cve, 1);
            if (visited.find(cve_heh1) != visited.end()
            ||  visited.find(cve_heh2) != visited.end()) {
                valance++;
            }
        }
        if (valance != 2) {
            return false;
        }
    }
    
    // Forming a oriented vector of joint boundary
    HalfedgeHandle starting_heh = *visited.begin();
    HalfedgeHandle current_heh = starting_heh;
    do {
        joint_boundary.push_back(current_heh);
        HalfedgeHandle next_heh = current_heh;
        do {
            next_heh = m_mesh->opposite_halfedge_handle(
                       m_mesh->next_halfedge_handle(next_heh));
        } while (next_heh != current_heh
              && visited.find(m_mesh->opposite_halfedge_handle(next_heh)) == visited.end());
        
        EdgeHandle eh = m_mesh->edge_handle(next_heh);
        current_heh = m_mesh->opposite_halfedge_handle(next_heh);
    } while (current_heh != starting_heh);
    
    return true;
}

void PlushPatternGenerator::optimize_patches(double threshold, bool step) {
    isJobCanceled = false;
    
    std::map< std::vector<HalfedgeHandle>, double > &boundary_distortion = m_mesh->property(joint_boundary_distortion_handle);
    std::map< std::vector<HalfedgeHandle>, double > &boundary_area = m_mesh->property(joint_boundary_area_handle);

    VertexHandle dummy;
    WeightFunctor weightFunctor(m_mesh, dummy, NULL);

    std::set<EdgeHandle> *seams = getSeams();
//    std::set<VertexHandle> intersection_points;
//    get_intersection_points(seams, intersection_points);
//    MeshSelection::clearEdgeSelection(m_mesh);

    bool cant_be_merged_anymore = false;
    while (!cant_be_merged_anymore) {
        cant_be_merged_anymore = true;
        
        std::vector< std::vector<HalfedgeHandle> > heh_segments;
//        get_segments_from_seams(heh_segments, intersection_points);
        get_segments_from_seams(heh_segments);
        
        // Sort from high to low
        typedef std::pair<std::vector<HalfedgeHandle>, double> Segment_pair;
        auto comparator = [](Segment_pair p1, Segment_pair p2) ->bool {
            return p1.second < p2.second;
        };
        std::priority_queue<Segment_pair, std::vector<Segment_pair>, decltype(comparator)> sorted_segments(comparator);
        for (auto heh_segment : heh_segments) {
            bool different_texture = false;
            for (HalfedgeHandle heh : heh_segment) {
                
                // If this segment contains different texture, it should not be merge
                if (is_different_texture(m_mesh, heh)) {
                    different_texture = true;
                    break;
                }
            }
            if (!different_texture) {
                double weight = weightFunctor(heh_segment);
                sorted_segments.push(std::make_pair(heh_segment, weight));
            }
        }
        
        int iterations = 0;
        int num_original_segments = sorted_segments.size();
        int count2 = 0;
        while (!sorted_segments.empty()) {
            if (isJobCanceled) {
                emit log(LOGINFO, "Merging canceled.");
                return;
            }

            Segment_pair pair = sorted_segments.top();
            std::vector<HalfedgeHandle> segment = pair.first;
            sorted_segments.pop();

            bool can_be_merge = true;
            double A = 0, L2_M = 0;
            
            std::vector<HalfedgeHandle> boundaries[3];
            // First two boundaries are the left & right boundary.
            get_adjacent_boundary(*segment.begin(), *seams, boundaries[0]);
            get_adjacent_boundary(m_mesh->opposite_halfedge_handle(*segment.begin()), *seams, boundaries[1]);
            
            // The third one is the joint boundary of the first two.
            bool is_valid = get_joint_boundary(boundaries[0], boundaries[1], *seams, boundaries[2]);

            // It's possible that after merging, some segment becomes "isolated" and will not form a valid joint boundary
            // In such case, just discard this segment
            double max_area_diff[3];
            double sum_area[3];
            if (!is_valid) {
                can_be_merge = true;
            } else {
                for (int i = 0; i < 3; i++) {
                    if (boundary_distortion.find(boundaries[i]) == boundary_distortion.end()) {
                        TriMesh adj_submesh;
                        extract_mesh_with_boundary(&adj_submesh, m_mesh->face_handle(*boundaries[i].begin()), &boundaries[i]);
                        std::map<VertexHandle, OpenMesh::Vec3d> boundaryPosition;
                        calcLPFB(&adj_submesh, &boundaryPosition);
                        calcInteriorPoints(&adj_submesh, &boundaryPosition);
                        calcDistortion(adj_submesh);
                        
                        max_area_diff[i] = std::numeric_limits<double>::min();
                        sum_area[i] = 0;
                        for (FaceHandle f : adj_submesh.faces()) {
                            sum_area[i] += adj_submesh.calc_sector_area(adj_submesh.halfedge_handle(f));

                            FaceHandle original_f = get_original_handle(&adj_submesh, f);
                            max_area_diff[i] = max(max_area_diff[i], abs(m_mesh->property(distortionFHandle, original_f)));
                        }
                        boundary_distortion[boundaries[i]] = max_area_diff[i];
                        boundary_area[boundaries[i]] = sum_area[i];
                    } else {
                        max_area_diff[i] = boundary_distortion[boundaries[i]];
                        sum_area[i] = boundary_area[boundaries[i]];
                    }
                }
//                    double sum_area_3D = 0;
//                    double sum_area_2D = 0;
//                    double sum_L2_Ti_mul_A_prime_Ti = 0;
//                    auto inverse_mapping = getInverseMappingHandle(&merged_mesh);
//                    for (FaceHandle f : merged_mesh.faces()) {
//                        FaceHandle original_f = get_original_handle(&merged_mesh, f);
//                        double area_3D = m_mesh->calc_sector_area(m_mesh->halfedge_handle(original_f));
//                        sum_area_3D += area_3D;
//                        double area_2D = merged_mesh.calc_sector_area(merged_mesh.halfedge_handle(f));
//                        sum_area_2D += area_2D;
//                    
//                        // Texture stretch from
//                        // P.V. Sander, J. Snyder, S.J. Gortler, and H. Hoppe,
//                        // “Texture Mapping Progressive Meshes,” Proc. ACM SIGGRAPH ’01, pp. 409-416, 2001
//                        TriMesh::ConstFaceVertexIter cfv_it = merged_mesh.cfv_iter(f);
//                        TriMesh::Point p1 = merged_mesh.point(*cfv_it);
//                        TriMesh::Point q1 = m_mesh->point(merged_mesh.property(inverse_mapping, *cfv_it++));
//                        TriMesh::Point p2 = merged_mesh.point(*cfv_it);
//                        TriMesh::Point q2 = m_mesh->point(merged_mesh.property(inverse_mapping, *cfv_it++));
//                        TriMesh::Point p3 = merged_mesh.point(*cfv_it);
//                        TriMesh::Point q3 = m_mesh->point(merged_mesh.property(inverse_mapping, *cfv_it++));
//                        
//                        OpenMesh::Vec3d Ss = (q1 * (p2[1] - p3[1]) + q2 * (p3[1] - p1[1]) + q3 * (p1[1] - p2[1]))/ (2 * area_2D);
//                        OpenMesh::Vec3d St = (q1 * (p3[0] - p2[0]) + q2 * (p1[0] - p3[0]) + q3 * (p2[0] - p1[0]))/ (2 * area_2D);
//                        double a = Ss | Ss;
//                        double c = St | St;
//                        double L2_Ti = sqrt((a + c)/2);
//                        sum_L2_Ti_mul_A_prime_Ti += L2_Ti * L2_Ti * area_3D;
//                    }
//                    L2_M = sqrt(sum_L2_Ti_mul_A_prime_Ti/sum_area_3D);
//                    A = abs(sum_area_2D - sum_area_3D)/sum_area_3D;
                
                // Now we should have the distortion of these boundary
                // The joint boundary is not developable if it satisfied all of the followings:
                // 1. max_area_diff greater than threshold
                // 2. The distortion doesn't increase too much (This usually happens when merging large patch with tiny patch)
                if (max_area_diff[2] > threshold) {
                    if (min(sum_area[0], sum_area[1])/max(sum_area[0], sum_area[1]) < 0.01) {
                        can_be_merge = true;
                    } else {
                        can_be_merge = false;
                    }
                }
            }
            
            iterations++;
            emit setJobState((double)iterations/num_original_segments * 100);

            // visualization
            if (can_be_merge && is_valid) {
                for (HalfedgeHandle heh : prevBoundary) {
                    EdgeHandle eh = m_mesh->edge_handle(heh);
                    m_mesh->set_color(eh, TriMesh::Color(1,1,1,0));
                    m_mesh->status(eh).set_selected(true);
                }
                for (HalfedgeHandle heh : prevSegment) {
                    EdgeHandle eh = m_mesh->edge_handle(heh);
                    m_mesh->set_color(eh, TriMesh::Color(1,1,1,0));
                }
                for (HalfedgeHandle heh : boundaries[2]) {
                    EdgeHandle eh = m_mesh->edge_handle(heh);
                    m_mesh->set_color(eh, TriMesh::Color(0,1,0,1));
                    m_mesh->status(eh).set_selected(false);
                }
                prevBoundary = boundaries[2];
                prevSegment = segment;
                
                for (HalfedgeHandle heh : segment) {
                    EdgeHandle eh = m_mesh->edge_handle(heh);
                    seams->erase(eh);
                    m_mesh->set_color(eh, TriMesh::Color(0,0,1,1));
                    m_mesh->status(eh).set_selected(false);
                }
                
                cant_be_merged_anymore = false;
                emit updateView();
                emit log(LOGINFO, QString::number(m_mesh->property(merge_iterations_handle)));
                m_mesh->property(merge_iterations_handle)++;

                if (step && m_mesh->property(merge_iterations_handle) % 10 == 0) {
                    return;
                } else {
                    break;
                }
            }
        }
    }
    
    // Assign color for segment
    std::vector< std::vector<HalfedgeHandle> > heh_segments;
    get_segments_from_seams(heh_segments);

    int count = 0;
    for (auto heh_segment : heh_segments) {
        for (HalfedgeHandle heh : heh_segment) {
            EdgeHandle eh = m_mesh->edge_handle(heh);
            m_mesh->property(segment_no_handle, eh) = count;
        }
        count = (count+20) % 715;
    }
    
//    // Remove isolated segments after merging
//    std::vector< std::vector<EdgeHandle> > segments;
//    get_segments_from_seams(segments);
//    bool changed = false;
//    do {
//        changed = false;
//        segments.clear();
//        get_segments_from_seams(segments);
//        for (std::vector<EdgeHandle> segment : segments) {
//            std::vector<HalfedgeHandle> joint_boundary;
//
//            // Find the joint boundary set
//            std::queue<HalfedgeHandle> queue;
//            std::set<HalfedgeHandle> visited;
//            HalfedgeHandle heh1 = m_mesh->halfedge_handle(*segment.begin(), 0);
//            HalfedgeHandle heh2 = m_mesh->halfedge_handle(*segment.begin(), 1);
//            queue.push(heh1);
//            queue.push(heh2);
//            while (!queue.empty()) {
//                HalfedgeHandle heh = queue.front();
//                queue.pop();
//                
//                VertexHandle v = m_mesh->to_vertex_handle(heh);
//                HalfedgeHandle current_heh = heh;
//                do {
//                    current_heh = m_mesh->opposite_halfedge_handle(
//                                                                   m_mesh->next_halfedge_handle(current_heh));
//                } while (current_heh != heh && seams->find(m_mesh->edge_handle(current_heh)) == seams->end());
//                
//                // It seems that this is the end of an one-way seam. Don't worry, we'll continue by reversing direction
//                if (current_heh == heh) { }
//                
//                current_heh = m_mesh->opposite_halfedge_handle(current_heh);
//                if (visited.find(current_heh) == visited.end()) {
//                    queue.push(current_heh);
//                    visited.insert(current_heh);
//                }
//            }
//            // Elliminate segment from candidate
//            for (EdgeHandle eh : segment) {
//                visited.erase(m_mesh->halfedge_handle(eh, 0));
//                visited.erase(m_mesh->halfedge_handle(eh, 1));
//            }
//            
//            // It's possible that after merging, some segment becomes "isolated" and will not form a valid joint boundary
//            // In such case, just discard this segment
//            if (visited.empty()) {
//                for (EdgeHandle eh : segment) {
//                    seams->erase(eh);
//                    m_mesh->status(eh).set_selected(false);
//                }
//                changed = true;
//                break;
//            }
//        }
//    } while (changed);
}