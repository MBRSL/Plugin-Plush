//
//  Geodesic.cc
//  OpenFlipper
//
//  Created by 饒亢 on 2015/1/13.
//
//
#include "OpenMesh_Boost_Wrapper.hh"
#include "PlushPatternGenerator.hh"

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/boykov_kolmogorov_max_flow.hpp>
#include <boost/graph/copy.hpp>

#include <MeshTools/MeshSelectionT.hh>
#include <queue>

#include <QTextStream>
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

bool PlushPatternGenerator::calcSeams(std::vector<VertexHandle> selectedVertices,
                                      double developable_threshold,
                                      int limitNum,
                                      bool elimination,
                                      bool allPath) {
    std::map< std::vector<HalfedgeHandle>, double > &joint_boundary_distortion = m_mesh->property(joint_boundary_distortion_handle);
    joint_boundary_distortion.clear();

    std::vector<TriMesh> subMeshes;

    std::set<EdgeHandle> *seams = getSeams();
    if (seams) {
        seams->clear();
    }
    
//    calcCircularSeams(m_mesh, developable_threshold, false);
    calcStructralSeams(m_mesh, selectedVertices, limitNum, elimination, allPath);
//    calcLocalSeams(m_mesh, developable_threshold);
    return true;
}

/**
 @brief Generate local seams for non-developable sub-meshes

 The result is append to seamsHandle property of m_mesh.
 @param <#parameter#>
 @return <#retval#>
 @retval <#meaning#>
 */
bool PlushPatternGenerator::calcLocalSeams(TriMesh *sub_mesh, double developable_threshold) {
    unsigned int num_prev_seam = 0;
    std::set<EdgeHandle> *seams = getSeams();
    std::vector<TriMesh> sub_meshes;
    std::set<std::vector<HalfedgeHandle> > calculated;
    
    do {
        num_prev_seam = seams->size();
        
        std::vector< std::vector<HalfedgeHandle> > closed_seams;
        get_closed_boundaries_of_seams(&closed_seams, seams);
        
        for (std::vector<HalfedgeHandle> closed_seam : closed_seams) {
            if (calculated.find(closed_seam) != calculated.end()) {
                continue;
            }
            calculated.insert(closed_seam);
            
            TriMesh sub_mesh;
            extract_mesh_with_boundary(&sub_mesh, m_mesh->face_handle(*closed_seam.begin()), seams);
            
            std::map<VertexHandle, OpenMesh::Vec3d> boundaryPosition;
            calcLPFB(&sub_mesh, &boundaryPosition);
            calcInteriorPoints(&sub_mesh, &boundaryPosition);
            calcDistortion(sub_mesh);
            
            calcCircularSeams(&sub_mesh, developable_threshold, true);
        }
        seams = getSeams();
    } while (seams->size() != num_prev_seam);
    
    return true;
}

/**
 @brief Generate circular seams that are not possible to get with calcSeams.
 The step of this function:
 1. Find vertices with high curvature (curvature threshold = mean + std).
 2. Categorize these vertices using connected component.
 3. Discard component with only one boundary.
 4. Calculate a max flow problem using faces correspond to vertices component.
 5. The resulting cut corresponds to the seam edge.

 The result is append to seamsHandle property of m_mesh.
 @param <#parameter#>
 @return <#retval#>
 @retval <#meaning#>
 */
bool PlushPatternGenerator::calcCircularSeams(TriMesh *mesh, double threshold, bool is_local_seam) {
    OpenMesh::VPropHandleT<VertexHandle> inverseMapping = getInverseMappingHandle(mesh);
    // Calculate mean & std
    double sum_curvature = 0;
    double sum_curvature_sqr = 0;
    int count = 0;
    for (VertexHandle v : mesh->vertices()) {
        VertexHandle original_v = get_original_handle(mesh, v);
        double curvature = m_mesh->property(maxCurvatureHandle, original_v);
        if (curvature < 0) {
            sum_curvature += curvature;
            sum_curvature_sqr += curvature * curvature;
            count++;
        }
    }
    double mean = sum_curvature/count;
    double std = sqrt(sum_curvature_sqr/mesh->n_vertices() - mean*mean);
    if (!is_local_seam) {
        threshold = mean/2;
    }
    
    double avg_edge_length = 0;
    for (EdgeHandle eh : mesh->edges()) {
        avg_edge_length += mesh->calc_edge_length(eh);
    }
    avg_edge_length /= mesh->n_edges();

    // Select vertices with curvature greater than threshold
    auto curvatureFilterFunctor = [&](VertexHandle v) -> bool {
        VertexHandle original_v = get_original_handle(mesh, v);
        double curvature = m_mesh->property(maxCurvatureHandle, original_v);
        return curvature < threshold;
    };

    auto distortionFilterFunctor = [&](VertexHandle v) -> bool {
        double distortion = 0;
        for (const FaceHandle cvf : mesh->vf_range(v)) {
            FaceHandle original_f = get_original_handle(mesh, cvf);
            distortion += m_mesh->property(distortionFHandle, original_f);
        }
        return abs(distortion/mesh->valence(v)) > threshold;
    };

    struct Filter {
    private:
        TriMesh *m_mesh;
        std::function<bool(VertexHandle v)> m_functor;
    public:
        Filter() {}
        Filter(TriMesh *mesh, std::function<bool(VertexHandle v)> functor) : m_mesh(mesh), m_functor(functor){}
        bool operator()(const VertexHandle& v) const {
            return m_functor(v);
        }
    };

    Filter filter;
    if (is_local_seam) {
        filter = Filter(mesh, distortionFilterFunctor);
    } else {
        filter = Filter(mesh, curvatureFilterFunctor);
    }
    typedef boost::filtered_graph<TriMesh, boost::keep_all, Filter> Filtered_graph;
    Filtered_graph curvature_filtered_mesh(*mesh, boost::keep_all(), filter);
    
    typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
    VertexHandle, HalfedgeHandle > Graph;
    Graph curvature_graph;

    struct Graph_copier {
    public:
        Graph *m_g_dst;
        Graph_copier(Graph *g_dst) : m_g_dst(g_dst) {}
        void operator()(const Filtered_graph::vertex_descriptor v_src, Graph::vertex_descriptor &v_dst) {
            (*m_g_dst)[v_dst] = v_src;
        }
        void operator()(const Filtered_graph::edge_descriptor e_src, Graph::edge_descriptor &e_dst) {
            (*m_g_dst)[e_dst] = e_src;
        }
    };
    Graph_copier graph_copier(&curvature_graph);
    
    // Split vertices into connected component
    // There seems to be some problems when applying connected components algorithm on filtered graph.
    // So we copy first it to a regular graph then apply the algorithm.
    std::map<int, int> componentMap;
    boost::copy_graph(curvature_filtered_mesh, curvature_graph,
                      boost::vertex_copy(graph_copier)
                      .edge_copy(graph_copier));
    int nComponents = boost::connected_components(curvature_graph, boost::make_assoc_property_map(componentMap));
    
    std::vector< std::set<VertexHandle> > components;
    components.resize(nComponents);
    for (size_t i = 0; i < num_vertices(curvature_graph); i++) {
        VertexHandle v = curvature_graph[i];
        components[componentMap[i]].insert(v);
    }
    
    // For each component
    for (auto component_it = components.begin(); component_it != components.end(); component_it++) {
        std::set<VertexHandle> component = *component_it;

        std::vector< std::set<VertexHandle> > boundariesVertices(2);

        // Use expand-shrink to form a much complete selection without holes
        const int nRings = 5;
        expandVertice(mesh, component, nRings);
        for (VertexHandle v : component) {
            for (const VertexHandle cvv : mesh->vv_range(v)) {
                // Any vertex connecting to non-component vertex is treat as boundary vertex
                if (component.find(cvv) == component.end()) {
                    boundariesVertices[0].insert(v);
                }
            }
        }
        std::set<VertexHandle> component_shrinked = component;
        shrinkVertice(mesh, component_shrinked, nRings);
        boundariesVertices[1].insert(component_shrinked.begin(), component_shrinked.end());
        
        // Visualize
//        if (is_local_seam) {
//            for (VertexHandle v : component) {
//                VertexHandle original_v = get_original_handle(mesh, v);
//                m_mesh->status(original_v).set_selected(true);
//            }
//            for (VertexHandle v : boundariesVertices[1]) {
//                VertexHandle original_v = get_original_handle(mesh, v);
//                m_mesh->status(original_v).set_selected(true);
//            }
//        }

//        std::vector< std::set<VertexHandle> > boundariesVertices;
//        boundariesVertices.resize(nBoundaries);
//        for (size_t i = 0; i < num_vertices(boundary_graph); i++) {
//            VertexHandle v = boundary_graph[i];
//            boundariesVertices[boundaryMap[i]].insert(v);
//        }
        
        // Find faces corresponds to boundariesVertices, also categorize all faces into first two boundary
        std::vector<FaceHandle> facesCandidates;
        std::vector< std::set<FaceHandle> > boundariesFaces(2);
        for (const FaceHandle cf : mesh->faces()) {
            TriMesh::ConstFaceVertexIter cfv_it = mesh->cfv_iter(cf);
            VertexHandle v1 = *cfv_it++;
            VertexHandle v2 = *cfv_it++;
            VertexHandle v3 = *cfv_it++;

            if (component.find(v1) != component.end()
            &&  component.find(v2) != component.end()
            &&  component.find(v3) != component.end()) {
                facesCandidates.push_back(cf);
                
                for (int i = 0; i < 2; i++) {
                    if (boundariesVertices[i].find(v1) != boundariesVertices[i].end()
                    ||  boundariesVertices[i].find(v2) != boundariesVertices[i].end()
                    ||  boundariesVertices[i].find(v3) != boundariesVertices[i].end()) {
                    // This face contains boundary vertices, add it to source/drain
                        boundariesFaces[i].insert(cf);
                    }
                }
            }
        }
                
        // Create a graph for min cut
        typedef unsigned int Id;
        typedef double EdgeWeight;
        typedef boost::adjacency_list_traits<boost::vecS, boost::vecS, boost::directedS> Traits;
        typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS,
            boost::property < boost::vertex_index_t, Id,
        // Stores inverse mapping here
            boost::property < boost::vertex_owner_t, FaceHandle,
            boost::property < boost::vertex_color_t, boost::default_color_type,
            boost::property < boost::vertex_distance_t, EdgeWeight,
            boost::property < boost::vertex_predecessor_t, Traits::edge_descriptor > > > > >,
        
            boost::property<boost::edge_index_t, Id,
            boost::property<boost::edge_capacity_t, EdgeWeight,
            boost::property<boost::edge_residual_capacity_t, EdgeWeight,
            boost::property<boost::edge_reverse_t, Traits::edge_descriptor> > > > > FlowGraph;

        FlowGraph flow_graph;
        auto capacity_map = get(boost::edge_capacity, flow_graph);
        auto reverse_edge_map = get(boost::edge_reverse, flow_graph);

        auto addBidirectionalEdge = [&](Id source, Id target, EdgeWeight capacity) {
            int nEdges = boost::num_edges(flow_graph);
            FlowGraph::edge_descriptor edge = boost::add_edge(source, target, nEdges+1, flow_graph).first;
            FlowGraph::edge_descriptor reverseEdge = boost::add_edge(target, source, nEdges+2, flow_graph).first;
            
            reverse_edge_map[edge] = reverseEdge;
            reverse_edge_map[reverseEdge] = edge;

            capacity_map[edge] = capacity;
            capacity_map[reverseEdge] = capacity;
        };
        
        std::map<FaceHandle, Id> idMapping;
        for (auto f_it = facesCandidates.begin(); f_it != facesCandidates.end(); f_it++) {
            Id vId = boost::add_vertex(flow_graph);
            put(boost::vertex_owner, flow_graph, vId, *f_it);
            idMapping.emplace(*f_it, vId);
        }
        std::set<FaceHandle> visitedF;
        for (FaceHandle f : facesCandidates) {
            visitedF.insert(f);
            for (const HalfedgeHandle cfh : mesh->fh_range(f)) {
                if (!mesh->is_boundary(mesh->opposite_halfedge_handle(cfh))) {
                    FaceHandle neighborF = mesh->opposite_face_handle(cfh);
                    if (visitedF.find(neighborF) == visitedF.end()
                    &&  idMapping.find(neighborF) != idMapping.end()) {
                        Id source = idMapping[f];
                        Id target = idMapping[neighborF];
                        
                        double weight = 1e9;
                        
                        VertexHandle v1 = mesh->from_vertex_handle(cfh);
                        VertexHandle v2 = mesh->to_vertex_handle(cfh);
                        if (is_local_seam) {
                            double distortion1 = 0;
                            for (const VertexHandle cvv : mesh->vv_range(v1)) {
                                double distortion_cvv = 0;
                                for (const FaceHandle cvf : mesh->vf_range(cvv)) {
                                    distortion_cvv += m_mesh->property(distortionFHandle, get_original_handle(mesh, cvf));
                                }
                                distortion1 += distortion_cvv/mesh->valence(cvv);
                            }
                            distortion1 /= m_mesh->valence(get_original_handle(mesh, v1));
                            distortion1 += m_mesh->property(distortionVHandle, get_original_handle(mesh, v1));
                            double distortion2 = 0;
                            for (const VertexHandle cvv : mesh->vv_range(v2)) {
                                double distortion_cvv = 0;
                                for (const FaceHandle cvf : mesh->vf_range(cvv)) {
                                    FaceHandle original_f = get_original_handle(mesh, cvf);
                                    distortion_cvv += m_mesh->property(distortionFHandle, get_original_handle(mesh, cvf));
                                }
                                distortion2 += distortion_cvv/mesh->valence(cvv);
                            }
                            distortion2 /= m_mesh->valence(get_original_handle(mesh, v2));
                            distortion2 += m_mesh->property(distortionVHandle, get_original_handle(mesh, v2));
                            weight = 1/(1 + abs(distortion1) + abs(distortion2));
                        } else {
                            const double alpha = 0.8;
//                        double angWij = pow((1-(1-abs(cos(theta)))*convexityFac),2);
//                        double angWij = pow(cos(theta),2);

                            double curvature1 = m_mesh->property(maxCurvatureHandle, get_original_handle(mesh, v1));
                            double curvature2 = m_mesh->property(maxCurvatureHandle, get_original_handle(mesh, v2));
                            if (curvature1 < 0 && curvature2 < 0) {
                                double curvatureWeight = abs((2 * mean)/(curvature1 + curvature2));
                                
                                double theta = mesh->calc_dihedral_angle(cfh);
                                double angleWeight = abs(cos(theta));
                                
    //                            double edge_length_weight = mesh->calc_edge_length(*cfh_it)/avg_edge_length;
                                weight = alpha * curvatureWeight + (1-alpha) * angleWeight;
                            }
                        }
                        addBidirectionalEdge(source, target, weight);
                    }
                }
            }
        }
        
        // Insert source & drain node
        for (size_t i = 0; i < boundariesFaces.size(); i++) {
            // Id of Source & drain are facesCandidates.size(), facesCandidates.size()+1 respectively.
            Id source = facesCandidates.size()+i;
            for (auto f_it = boundariesFaces[i].begin(); f_it != boundariesFaces[i].end(); f_it++) {
                Id target = idMapping[*f_it];
                addBidirectionalEdge(source, target, 1e9);
            }
        }

        // Because floating point error will cause assertion error in max flow, we just temporarily disable it.
#ifdef NDEBUG
#undef NDEBUG
#define _UN_NDEBUG_
#endif
        boost::boykov_kolmogorov_max_flow(flow_graph, facesCandidates.size(), facesCandidates.size()+1);
#ifdef _UN_NDEBUG_
#undef _UN_NDEBUG_
#define NDEBUG
#endif
        // Find target cut(edges)
        std::vector<EdgeHandle> final_edges;
        FlowGraph::edge_iterator e_it, e_ite;
        for (tie(e_it, e_ite) = boost::edges(flow_graph); e_it != e_ite; e_it++) {
            Id v1 = boost::source(*e_it, flow_graph);
            Id v2 = boost::target(*e_it, flow_graph);
            
            // Ignore virtual source/drain
            if (v1 < facesCandidates.size() &&  v2 < facesCandidates.size()) {
                boost::default_color_type color1 = get(boost::vertex_color, flow_graph, v1);
                boost::default_color_type color2 = get(boost::vertex_color, flow_graph, v2);
            
                if (color1 != color2
                && (color1 == boost::black_color || color2 == boost::black_color)) {
                    FaceHandle f1 = get(boost::vertex_owner, flow_graph, v1);
                    FaceHandle f2 = get(boost::vertex_owner, flow_graph, v2);
                    for (const HalfedgeHandle cfh : mesh->fh_range(f1)) {
                        if (!mesh->is_boundary(mesh->opposite_halfedge_handle(cfh))
                       &&   mesh->opposite_face_handle(cfh) == f2) {
                            HalfedgeHandle original_heh = get_original_handle(mesh, cfh);
                            final_edges.push_back(m_mesh->edge_handle(original_heh));
                        }
                    }
                }
            }
        }
        OpenMesh::MPropHandleT< std::set<EdgeHandle> > seamsHandle = getSeamsHandle(m_mesh);
        std::set<EdgeHandle> &seams = m_mesh->property(seamsHandle);
        seams.insert(final_edges.begin(), final_edges.end());
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
bool PlushPatternGenerator::calcStructralSeams(TriMesh *mesh,
                                               std::vector<VertexHandle> selectedVertices,
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
    OpenMesh::MPropHandleT< std::set<EdgeHandle> > seamsHandle = getSeamsHandle(m_mesh);
    std::set<EdgeHandle> &seams = m_mesh->property(seamsHandle);
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
            bool edgeExist = getEdge(m_mesh, eh, path[i-1], path[i]);
            assert(edgeExist);
            seams.insert(eh);
        }
        
        QString msg = QString("Weight of path #%1 from %2 to %3: %4").arg(count+1).arg((path.begin())->idx()).arg((path.end()-1)->idx()).arg(it->first);
        emit log(LOGINFO, msg);
    }
    
    // Add texture seams that are not possible to be generated in the above proccess
//    for (EdgeHandle eh : m_mesh->edges()) {
//        if (is_different_texture(m_mesh, eh)) {
//            seams.insert(eh);
//        }
//    }
    
    return true;
}

bool PlushPatternGenerator::load_seams() {
    QFile file(m_meshName + "_seams.txt");
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
        emit log(LOGERR, QString("Unable to read file %1").arg(m_meshName + "_seams.txt"));
        return 0;
    }
    
    QTextStream in(&file);

    OpenMesh::MPropHandleT< std::set<EdgeHandle> > seamsHandle = PlushPatternGenerator::getSeamsHandle(m_mesh);
    std::set<EdgeHandle> &seams = m_mesh->property(seamsHandle);

    QString eh_str;
    while(!in.atEnd()) {
        eh_str = in.readLine();
        bool ok;
        int eh_id = eh_str.toInt(&ok);
        if (ok) {
            EdgeHandle eh = m_mesh->edge_handle(eh_id);
            m_mesh->status(eh).set_selected(true);
            seams.insert(eh);
        }
    }
    
    emit updateView();
    return true;
}

bool PlushPatternGenerator::save_seams() {
    OpenMesh::MPropHandleT< std::set<EdgeHandle> > seamsHandle = PlushPatternGenerator::getSeamsHandle(m_mesh);
    std::set<EdgeHandle> &seams = m_mesh->property(seamsHandle);
    
    // Prepare file for saving data
    QFile file(m_meshName + "_seams.txt");
    file.open(QIODevice::WriteOnly | QIODevice::Text);
    QTextStream out(&file);
    
    for (EdgeHandle eh : seams) {
        out << eh.idx() << "\n";
    }
    file.close();
    
    return true;
}
