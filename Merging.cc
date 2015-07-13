#include "OpenMesh_Boost_Wrapper.hh"
#include "PlushPatternGenerator.hh"
#include "WeightFunctor.hh"
#include "FilteredTriMesh.hh"

#include <queue>
#include <unordered_set>
/*
 @brief Merge two adjacent patches with seam_segment as boundary
 @param seam_segment_idx
    The idx of seam segment which is used to merged these two patch. This is the edge# of subMesh_graph
    If you don't want to use it, just pass -1.
 @param seam_segment 
    Indicate which boundary edge should be merged (mark as non-boundary) in the merged patch.
 */
FilteredTriMesh PlushPatternGenerator::merge_patch(FilteredTriMesh &patch1,
                                                   FilteredTriMesh &patch2,
                                                   int seam_segment_idx,
                                                   std::set<EdgeHandle> &seam_segment) {
    std::set<FaceHandle> faces;
    std::set<EdgeHandle> boundary_edges;

    for (int i = 0; i < 2; i++) {
        FilteredTriMesh &patch = (i == 0) ? patch1 : patch2;
        for (EdgeHandle eh : patch.boundary_edges()) {
            // Ignore edges to be merged
            if (seam_segment.find(eh) == seam_segment.end()
            &&  !patch1.merged_edge_idx[eh.idx()]
            &&  !patch2.merged_edge_idx[eh.idx()]) {
                boundary_edges.insert(eh);
            }
        }
        for (FaceHandle fh : patch.faces()) {
            faces.insert(fh);
        }
    }
    
    FilteredTriMesh merged_patch(m_mesh, faces, boundary_edges);
    
    // Update idx information by union the set
    merged_patch.merged_subMesh_idx = patch1.merged_subMesh_idx | patch2.merged_subMesh_idx;
    merged_patch.merged_seam_idx = patch1.merged_seam_idx | patch2.merged_seam_idx;
    merged_patch.merged_edge_idx = patch1.merged_edge_idx | patch2.merged_edge_idx;
    
    if (seam_segment_idx >= 0) {
        merged_patch.merged_seam_idx[seam_segment_idx] = true;
    }

    merged_patch.n_merged_seams = merged_patch.merged_seam_idx.count();
    
    for (EdgeHandle eh : seam_segment) {
        merged_patch.merged_edge_idx[eh.idx()] = true;
    }
    return merged_patch;
}

/**
 @brief Construct sets of patches from sub-meshes separated by seams.
 Here we defined the term "sub-mesh" as the smallest mesh separated by seams and
 "patch" is generated by merging one adjacent sub-meshes.
 The process is illustrated as following:
 1. Treat each sub-mesh as a patch consist only one sub-mesh. These are the level 0 patches.
 2. For each level 0 patch, find its neigboring sub-mesh
    and merge them if the distortion is smaller than threshold.
    The merged patches are the level 1 patches.
 3. Repeat step 2 for next level until no patches in this level. i.e., no patches can be merged.
 @param threshold The distortion (area change of faces) threshold used for testing whether two patches can be merged
 */
void PlushPatternGenerator::construct_subsets(double threshold) {
    isJobCanceled = false;
    
    std::set<EdgeHandle> &seams = m_mesh->property(seams_handle);
    
    SubMesh_graph subMesh_graph = get_subMeshes_with_boundary(seams);
    // Record the best layout given Patch_idx
    std::map<FilteredTriMesh::Patch_idx, FilteredTriMesh> optimal_patches;

    // This stores which patch has been calculated
    std::set<FilteredTriMesh::Seam_idx> is_patch_calculated;
    
    // The merged result of this iteration, will be used as base patch in the next iteration
    std::map<FilteredTriMesh::Patch_idx, FilteredTriMesh> base_patches;
    
    // Construct first level of patches from sub-meshes
    // Compute the distortion of each patch
    int progress_counter = 0;
    SubMesh_graph::vertex_iterator v_it, v_ite;
    for (boost::tie(v_it, v_ite) = boost::vertices(subMesh_graph); v_it != v_ite; v_it++) {
        FilteredTriMesh &subMesh = boost::get(boost::vertex_owner, subMesh_graph, *v_it);
        optimal_patches.emplace(subMesh.merged_subMesh_idx, subMesh);
        base_patches.emplace(subMesh.merged_subMesh_idx, subMesh);

        std::map<HalfedgeHandle, OpenMesh::Vec3d> boundaryPosition;
        calcLPFB(subMesh, boundaryPosition);
        calcInteriorPoints(subMesh, boundaryPosition);
        calcDistortion(subMesh);
        subMesh.seam_score = 0;
        is_patch_calculated.insert(subMesh.merged_seam_idx);
        
        emit setJobState((double)progress_counter++/boost::num_vertices(subMesh_graph) * 100);
    }
    
    VertexHandle dummy;
    WeightFunctor weightFunctor(m_mesh, dummy, NULL,
                                m_distanceCoefficient,
                                m_textureCoefficient,
                                m_curvatureCoefficient,
                                m_skeletonCoefficient,
                                m_pathCoefficient
                                );

    // Construct larger patches from previous level
//    for (int level = 1; hierarchical_patches[level-1].size() > 0; level++) {
        //    for (int level = 1; level < 30; level++) {
    int level;
    for (level = 0; !base_patches.empty(); level++) {
        std::map<FilteredTriMesh::Patch_idx, FilteredTriMesh> next_base_patches;

        // Extend container
//        hierarchical_patches.resize(hierarchical_patches.size()+1);

        // For every patch, test if it can be merged with its neighbors
        progress_counter = 0;
        for (auto &pair: base_patches) {
            FilteredTriMesh &patch = pair.second;
            if (isJobCanceled) {
                emit log(LOGINFO, "Subset calculation canceled.");
                return;
            }
            
            // For each sub-mesh in this patch, find its neigboring sub-mesh for merging
            for (SubMesh_graph::vertex_descriptor v1 = patch.merged_subMesh_idx.find_first();
                 v1 != boost::dynamic_bitset<>::npos;
                 v1 = patch.merged_subMesh_idx.find_next(v1)) {
                // Use subMesh_grpah to find adjacent sub-meshes
                SubMesh_graph::out_edge_iterator e_it, e_ite;
                for (boost::tie(e_it, e_ite) = boost::out_edges(v1, subMesh_graph); e_it != e_ite; e_it++) {
                    int seam_segment_idx = boost::get(boost::edge_index, subMesh_graph, *e_it);
                    
                    // Only deal with edges that are not merged yet
                    if (!patch.merged_seam_idx[seam_segment_idx]) {
                        SubMesh_graph::vertex_descriptor v2 = boost::target(*e_it, subMesh_graph);
                        FilteredTriMesh &neigboring_subMesh = boost::get(boost::vertex_owner, subMesh_graph, v2);
                    
                        // Create a merged sub-mesh
                        boost::dynamic_bitset<> merged_seam_idx = patch.merged_seam_idx;
                        merged_seam_idx[seam_segment_idx] = true;
                        
                        if (is_patch_calculated.find(merged_seam_idx) == is_patch_calculated.end()) {
                            is_patch_calculated.insert(merged_seam_idx);

                            std::vector<HalfedgeHandle> seam_segment = boost::get(boost::edge_owner, subMesh_graph, *e_it);
                            std::set<EdgeHandle> seam_segment_set;
                            for (HalfedgeHandle heh : seam_segment) {
                                seam_segment_set.insert(m_mesh->edge_handle(heh));
                            }
                            FilteredTriMesh merged_subMesh = merge_patch(patch, neigboring_subMesh, seam_segment_idx, seam_segment_set);
                            
                            // check the guassian curvature of the segment to be merged
                            // if all curvature < threshold, then this segment can be merged
                            // note that boundary vertices do not need to have 0 curvature
                            bool non_developable = false;
                            for (HalfedgeHandle heh : seam_segment) {
                                VertexHandle v1 = merged_subMesh.from_vertex_handle(heh);
                                VertexHandle v2 = merged_subMesh.from_vertex_handle(heh);
                                
                                double gaussian_curvature1 = m_mesh->property(gaussianCurvatureHandle, v1);
                                double gaussian_curvature2 = m_mesh->property(gaussianCurvatureHandle, v2);
                                if ((!merged_subMesh.is_boundary(v1) && gaussian_curvature1 > 0.05) ||
                                    (!merged_subMesh.is_boundary(v2) && gaussian_curvature2 > 0.05)) {
                                    non_developable = true;
                                    break;
                                }
                            }
                            
                            if (!non_developable) {
                                // Distortion not calculated yet. Calculate now.
                                std::map<HalfedgeHandle, OpenMesh::Vec3d> boundaryPosition;
                                calcLPFB(merged_subMesh, boundaryPosition);
                                calcInteriorPoints(merged_subMesh, boundaryPosition);
                                calcDistortion(merged_subMesh);
                                
                                if (merged_subMesh.max_distortion < threshold) {
                                    merged_subMesh.seam_score = weightFunctor(seam_segment);
                                    if (patch.seam_score >= 0) {
                                        merged_subMesh.seam_score += patch.seam_score;
                                    }
                                    
                                    auto found_prev_patch = optimal_patches.find(merged_subMesh.merged_subMesh_idx);
                                    if (found_prev_patch != optimal_patches.end()) {
                                        double prev_seam_score = found_prev_patch->second.seam_score;
                                        if (prev_seam_score == -1 || merged_subMesh.seam_score < prev_seam_score) {
                                            found_prev_patch->second = merged_subMesh;

                                            // Replace record in next_base_patches to prevent from calculating more candidate patches
                                            auto found_next_base_patch = next_base_patches.find(merged_subMesh.merged_subMesh_idx);
                                            if (found_next_base_patch != next_base_patches.end()) {
                                                found_next_base_patch->second = merged_subMesh;
                                            } else {
                                                next_base_patches.emplace(merged_subMesh.merged_subMesh_idx, merged_subMesh);
                                            }
                                        }
                                    } else {
                                        optimal_patches.emplace(merged_subMesh.merged_subMesh_idx, merged_subMesh);
                                        next_base_patches.emplace(merged_subMesh.merged_subMesh_idx, merged_subMesh);
                                    }
                                }
                            }
                        }
                    }
                }
            }
            progress_counter++;
            emit setJobState((double)progress_counter/base_patches.size() * 100);
        }
        printf("Level %d: %lu patches\n", level, base_patches.size());
        
        base_patches = std::move(next_base_patches);
        next_base_patches.clear();
    }
    // Add result to property
    for (auto &pair : optimal_patches) {
        FilteredTriMesh &patch = pair.second;
        if (patch.n_merged_seams == level-1) {
            m_mesh->property(merged_patches_handle).push_back(patch);
        }
    }
}

void PlushPatternGenerator::optimize_patches(double threshold, bool step) {
    isJobCanceled = false;
    
    std::map< std::set<HalfedgeHandle>, double > &boundary_distortion = m_mesh->property(joint_boundary_distortion_handle);
    std::map< std::set<HalfedgeHandle>, double > &boundary_area = m_mesh->property(joint_boundary_area_handle);
    
    VertexHandle dummy;
    WeightFunctor weightFunctor(m_mesh, dummy, NULL,
                                m_distanceCoefficient,
                                m_textureCoefficient,
                                m_curvatureCoefficient,
                                m_skeletonCoefficient,
                                m_pathCoefficient
                                );
    
    std::set<EdgeHandle> &seams = m_mesh->property(seams_handle);
    
    std::set<VertexHandle> &intersection_points = m_mesh->property(intersection_points_handle);
    if (intersection_points.empty()) {
        intersection_points = get_intersection_points(seams);
    }
    
    bool cant_be_merged_anymore = false;
    while (!cant_be_merged_anymore) {
        cant_be_merged_anymore = true;
        
        std::vector< std::vector<HalfedgeHandle> > heh_segments;
//        heh_segments = get_halfedge_segments_from_seams(seams);
        heh_segments = get_halfedge_segments_from_seams(intersection_points, seams);
        
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
        while (!sorted_segments.empty()) {
            if (isJobCanceled) {
                emit log(LOGINFO, "Merging canceled.");
                return;
            }
            
            Segment_pair pair = sorted_segments.top();
            std::vector<HalfedgeHandle> segment = pair.first;
            sorted_segments.pop();
            
            bool can_be_merge = true;
//            double A = 0, L2_M = 0;
            double max_area_diff[3];
            double sum_area[3];
            std::set<HalfedgeHandle> boundaries[3];
            get_adjacent_boundary(*segment.begin(), seams, boundaries[0]);
            get_adjacent_boundary(m_mesh->opposite_halfedge_handle(*segment.begin()), seams, boundaries[1]);
            
            // The third one is the joint boundary of the first two.
            get_joint_boundary(boundaries[0], boundaries[1], segment, boundaries[2], seams);

            for (int i = 0; i < 3; i++) {
                if (boundary_distortion.find(boundaries[i]) == boundary_distortion.end()) {
                    FilteredTriMesh submesh = get_subMesh_with_boundary(m_mesh->face_handle(*boundaries[i].begin()), boundaries[i]);

                    std::map<HalfedgeHandle, OpenMesh::Vec3d> boundaryPosition;
                    calcLPFB(submesh, boundaryPosition);
                    calcInteriorPoints(submesh, boundaryPosition);
                    calcDistortion(submesh);
                
                    max_area_diff[i] = std::numeric_limits<double>::min();
                    sum_area[i] = 0;
                    for (FaceHandle f : submesh.faces()) {
                        sum_area[i] += submesh.calc_sector_area(submesh.halfedge_handle(f));
                        
                        double distortion = m_mesh->property(distortionFHandle, f);
                        max_area_diff[i] = max(max_area_diff[i], abs(distortion));
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
            
            iterations++;
            emit setJobState((double)iterations/num_original_segments * 100);
            
            if (!can_be_merge) {
                // Restore seams if not developable
                for (HalfedgeHandle heh : segment) {
                    seams.insert(m_mesh->edge_handle(heh));
                }
            } else {
                // visualization
                for (HalfedgeHandle heh : prevBoundary) {
                    EdgeHandle eh = m_mesh->edge_handle(heh);
                    m_mesh->set_color(eh, TriMesh::Color(1,1,1,0));
                    m_mesh->status(eh).set_selected(true);
//                    m_mesh->status(heh).set_selected(false);
                }
                for (HalfedgeHandle heh : prevSegment) {
                    EdgeHandle eh = m_mesh->edge_handle(heh);
                    m_mesh->set_color(eh, TriMesh::Color(1,1,1,0));
                    m_mesh->status(eh).set_selected(false);
                }

                prevBoundary.clear();
                if (boundaries[2].empty()) {
                    emit log(LOGERR, "No boundary is found in submesh[2]");
                    return;
                }
                for (HalfedgeHandle heh : boundaries[2]) {
                    EdgeHandle eh = m_mesh->edge_handle(heh);
                    m_mesh->set_color(eh, TriMesh::Color(0,1,0,1));
                    m_mesh->status(eh).set_selected(false);
//                    m_mesh->status(heh).set_selected(true);
                }
                prevBoundary = boundaries[2];
                
                for (HalfedgeHandle heh : segment) {
                    EdgeHandle eh = m_mesh->edge_handle(heh);
                    seams.erase(eh);
                    m_mesh->set_color(eh, TriMesh::Color(0,0,1,1));
                    m_mesh->status(eh).set_selected(false);
                }
                prevSegment = segment;
                
                cant_be_merged_anymore = false;
                
                emit updateView();
                
                if (step) {
                    return;
                } else {
                    break;
                }
            }
        }
    }
    
    // Assign color for segment
    std::vector< std::vector<HalfedgeHandle> > heh_segments = get_halfedge_segments_from_seams(seams);
    
    int count = 0;
    for (auto heh_segment : heh_segments) {
        for (HalfedgeHandle heh : heh_segment) {
            EdgeHandle eh = m_mesh->edge_handle(heh);
            m_mesh->property(segment_no_handle, eh) = count;
        }
        count = (count+20) % 255;
    }
}
