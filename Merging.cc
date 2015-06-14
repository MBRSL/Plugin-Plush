#include "OpenMesh_Boost_Wrapper.hh"
#include "PlushPatternGenerator.hh"
#include "WeightFunctor.hh"

#include <queue>
#include <unordered_set>
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
    
    std::set<EdgeHandle> *seams = getSeams();
    
    std::set<VertexHandle> &intersection_points = m_mesh->property(intersection_points_handle);
    if (intersection_points.empty()) {
        get_intersection_points(seams, intersection_points);
    }
    
    bool cant_be_merged_anymore = false;
    while (!cant_be_merged_anymore) {
        cant_be_merged_anymore = true;
        
        std::vector< std::vector<HalfedgeHandle> > heh_segments;
//        get_segments_from_seams(heh_segments, seams);
        get_segments_from_seams(heh_segments, intersection_points, seams);
        
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
            get_adjacent_boundary(*segment.begin(), *seams, boundaries[0]);
            get_adjacent_boundary(m_mesh->opposite_halfedge_handle(*segment.begin()), *seams, boundaries[1]);
            
            // The third one is the joint boundary of the first two.
            get_joint_boundary(boundaries[0], boundaries[1], segment, boundaries[2], *seams);

            for (int i = 0; i < 3; i++) {
                if (boundary_distortion.find(boundaries[i]) == boundary_distortion.end()) {
                    TriMesh submesh;
                    extract_mesh_with_boundary(&submesh, m_mesh->face_handle(*boundaries[i].begin()), &boundaries[i]);
                    std::map<VertexHandle, OpenMesh::Vec3d> boundaryPosition;
                    calcLPFB(&submesh, &boundaryPosition);
                    calcInteriorPoints(&submesh, &boundaryPosition);
                    calcDistortion(submesh);
                
                    max_area_diff[i] = std::numeric_limits<double>::min();
                    sum_area[i] = 0;
                    for (FaceHandle f : submesh.faces()) {
                        sum_area[i] += submesh.calc_sector_area(submesh.halfedge_handle(f));
                        
                        FaceHandle original_f = get_original_handle(&submesh, f);
                        double distortion = m_mesh->property(distortionFHandle, original_f);
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
                    seams->insert(m_mesh->edge_handle(heh));
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
                    seams->erase(eh);
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
    std::vector< std::vector<HalfedgeHandle> > heh_segments;
    get_segments_from_seams(heh_segments, seams);
    
    int count = 0;
    for (auto heh_segment : heh_segments) {
        for (HalfedgeHandle heh : heh_segment) {
            EdgeHandle eh = m_mesh->edge_handle(heh);
            m_mesh->property(segment_no_handle, eh) = count;
        }
        count = (count+20) % 255;
    }
}