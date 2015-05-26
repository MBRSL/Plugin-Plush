#include "OpenMesh_Boost_Wrapper.hh"
#include "LPFB.hh"
#include "PlushPatternGenerator.hh"

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/copy.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>

#include <cassert>
#include <iostream>

#include <queue>

using namespace Ipopt;

// constructor
LPFB_NLP::LPFB_NLP(TriMesh *mesh, std::map<VertexHandle, OpenMesh::Vec3d> *boundaryPosition) : m_mesh(mesh), m_boundaryPosition(boundaryPosition)
{
    std::vector< std::vector<HalfedgeHandle> > boundaries;
    PlushPatternGenerator::getBoundaryOfOpenedMesh(boundaries, m_mesh, true);
    
    assert(!boundaries.empty() && "No cut is found on given sub mesh.");
    
    // Because of unknown reason, the multiple source version of Dijkastra doesn't work.
    // We have to constrcut a new graph with an additional starting point to achieve it.
    typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
    boost::property<boost::vertex_owner_t, VertexHandle>,

    boost::property<boost::edge_weight_t, double,
    boost::property<boost::edge_owner_t, HalfedgeHandle> >
    > Graph;
    Graph multi_source_graph;
    
    struct Graph_copier {
    public:
        TriMesh *m_mesh;
        Graph *m_g_dst;
        Graph_copier(TriMesh *mesh, Graph *g_dst) : m_mesh(mesh), m_g_dst(g_dst) {}
        void operator()(const VertexHandle v_src, Graph::vertex_descriptor &v_dst) {
            put(boost::vertex_owner, *m_g_dst, v_dst, v_src);
        }
        void operator()(const HalfedgeHandle e_src, Graph::edge_descriptor &e_dst) {
            put(boost::edge_owner, *m_g_dst, e_dst, e_src);
            put(boost::edge_weight, *m_g_dst, e_dst, m_mesh->calc_edge_length(e_src));
        }
    };
    Graph_copier graph_copier(m_mesh, &multi_source_graph);

    // It seems that we can't directly appy copy_graph on TriMesh, but wrap it with filtered_graph does the trick
    typedef boost::filtered_graph<TriMesh, boost::keep_all> Filtered_graph;
    Filtered_graph filtered_mesh(*m_mesh, boost::keep_all());
    boost::copy_graph(filtered_mesh, multi_source_graph,
                      boost::vertex_copy(graph_copier)
                      .edge_copy(graph_copier));

    // We use the boundary with most vertices number as base boundary. Other boundaries are treated as inner boundaries.
    // The target set. Virtual cut searching will stop when reaching target set.
    std::vector<HalfedgeHandle> *baseBoundary = &boundaries[0];
    for (size_t j = 1; j < boundaries.size(); j++) {
        if (boundaries[j].size() > baseBoundary->size()) {
            baseBoundary = &boundaries[j];
        }
    }

    m_boundary3D.insert(m_boundary3D.begin(), baseBoundary->begin(), baseBoundary->end());

    // For inner boundaries, calculate virtual cut and insert them into appropriate position of m_boundary3D
    for (size_t i = 0; i < boundaries.size(); i++) {
        std::vector<HalfedgeHandle> &sourceBoundary = boundaries[i];
        if (sourceBoundary == *baseBoundary) {
            continue;
        }

        std::set<VertexHandle> sourceBoundaryVertices;
        for (auto e_it = sourceBoundary.begin(); e_it != sourceBoundary.end(); e_it++) {
            sourceBoundaryVertices.insert(m_mesh->from_vertex_handle(*e_it));
        }

        // we solve multi-source shoretest path problem by manually construct a new graph
        // with new additional vertex as starting point and inserting 0-weight edges
        // from starting point to other point on source boundary
        Graph::vertex_descriptor startingV = boost::add_vertex(multi_source_graph);
        Graph::vertex_iterator v_it, v_ite;
        for (tie(v_it, v_ite) = vertices(multi_source_graph); v_it != v_ite; v_it++) {
            VertexHandle v = get(boost::vertex_owner, multi_source_graph, *v_it);
            if (sourceBoundaryVertices.find(v) != sourceBoundaryVertices.end()) {
                Graph::edge_descriptor e = boost::add_edge(startingV, *v_it, multi_source_graph).first;
                put(boost::edge_weight, multi_source_graph, e, 0);
            }
        }

        std::vector<Graph::vertex_descriptor> predecessor_map(num_vertices(multi_source_graph));
        auto predecessor_pmap = boost::make_iterator_property_map(predecessor_map.begin(), get(boost::vertex_index, multi_source_graph));
        
        std::vector<double> distance_map(num_vertices(multi_source_graph));
        auto distance_pmap = boost::make_iterator_property_map(distance_map.begin(), get(boost::vertex_index, multi_source_graph));
        
        boost::dijkstra_shortest_paths(multi_source_graph, startingV,
                                       boost::predecessor_map(predecessor_pmap)
                                       .distance_map(distance_pmap)
                                      );
        
        // Don't forget to remove starting point so that multi_source_graph can be reuse for next boundary
        boost::clear_vertex(startingV, multi_source_graph);
        boost::remove_vertex(startingV, multi_source_graph);
        
        // Find the shortest distance on target boundary
        double minDistance = std::numeric_limits<double>::max();
        Graph::vertex_descriptor minDstV;
        for (size_t j = 0; j < baseBoundary->size(); j++) {
            VertexHandle v = m_mesh->to_vertex_handle(baseBoundary->at(j));
            if (distance_map[v.idx()] < minDistance) {
                minDistance = distance_map[v.idx()];
                minDstV = v.idx();
            }
        }
        assert(minDistance > 0);
        
        // Back tracing
        std::vector<HalfedgeHandle> virtualCut;
        std::vector<HalfedgeHandle> virtualCutReverse;
        Graph::vertex_descriptor currentV = minDstV, nextV = predecessor_map[minDstV];
        while (nextV != startingV) {
            assert(std::find(sourceBoundaryVertices.begin(), sourceBoundaryVertices.end(), m_mesh->vertex_handle(currentV)) == sourceBoundaryVertices.end());
            
            VertexHandle v_from = get(boost::vertex_owner, multi_source_graph, currentV);
            VertexHandle v_to = get(boost::vertex_owner, multi_source_graph, nextV);
            HalfedgeHandle heh;
            bool edgeExist = PlushPatternGenerator::getHalfedge(m_mesh, heh, v_from, v_to);
            assert(edgeExist);
            
            virtualCut.push_back(heh);
            virtualCutReverse.push_back(m_mesh->opposite_halfedge_handle(heh));
            
            currentV = nextV;
            nextV = predecessor_map[nextV];
        }
        std::reverse(virtualCutReverse.begin(), virtualCutReverse.end());
        
        // Find the insertion position (the index of starting point in source boundary)
        int insertionPosition = -1;
        for (size_t j = 0; j < sourceBoundary.size(); j++) {
            if ((unsigned long)m_mesh->from_vertex_handle(sourceBoundary[j]).idx() == currentV) {
                insertionPosition = j;
                break;
            }
        }
        assert(insertionPosition >= 0 && "The end point of virtual cut is not found in source boundary.");
        
        // Rotate sourceBoundary so that the 0th element matches insertion position
        std::rotate(sourceBoundary.begin(), sourceBoundary.begin() + insertionPosition, sourceBoundary.end());
        
        // Insert this boundary and virtual cut into overall boundary
        for (auto he_it = m_boundary3D.begin(); he_it != m_boundary3D.end(); he_it++) {
            HalfedgeHandle heh = *he_it;
            // Find the insertion position
            if (m_mesh->from_vertex_handle(heh) == m_mesh->from_vertex_handle(virtualCut[0])) {
                he_it = m_boundary3D.insert(he_it, virtualCut.begin(), virtualCut.end());
                he_it = m_boundary3D.insert(he_it+virtualCut.size(), sourceBoundary.begin(), sourceBoundary.end());
                he_it = m_boundary3D.insert(he_it+sourceBoundary.size(), virtualCutReverse.begin(), virtualCutReverse.end());
                
                // Check if this boundary is well connected together
                for (size_t edgeNo = 0; edgeNo < m_boundary3D.size(); edgeNo++) {
                    assert(m_mesh->to_vertex_handle(m_boundary3D[edgeNo])
                        == m_mesh->from_vertex_handle(m_boundary3D[(edgeNo+1)%m_boundary3D.size()]));
                }
                break;
            }
        }
    }
    
    // Calculating coincide pairs
    // At least the base boundary should be connected
    m_coincidentPair.push_back(std::make_pair(0, m_boundary3D.size()));
    // Other inner boundaries
    std::map<VertexHandle, int> vertexPosition;
    bool cycle = false;
    for (size_t position = 0; position < m_boundary3D.size(); position++) {
        VertexHandle v = m_mesh->from_vertex_handle(m_boundary3D[position]);
        
        // If this vertex is visited before, there is a loop. And we need to add it to coincide pairs
        if (vertexPosition.find(v) != vertexPosition.end()) {
            int prevPosition = vertexPosition[v];
            m_coincidentPair.push_back(std::make_pair(prevPosition, position));
            if (!cycle) {
                m_boundaryCoincidentPair.push_back(std::make_pair(prevPosition, position));
            }
            // We need to delete it first before update it.
            vertexPosition.erase(v);
            cycle = true;
        } else {
            cycle = false;
        }
        vertexPosition.emplace(v, position);
    }
    // There is at least one coincident pair on boundary
    m_boundaryCoincidentPair.push_back(std::make_pair(0, m_boundary3D.size()));
    
    // Calculate inner angles & edges length for later use.
    int n = m_boundary3D.size();
    for (int j = 0; j < n; j++) {
        double sumInnerAngle = PlushPatternGenerator::getSumInnerAngle(m_mesh, m_boundary3D[(j-1+n) % n], m_boundary3D[j]);
        m_innerAngle3D.push_back(sumInnerAngle);
        
        m_edgeLengths.push_back(m_mesh->calc_edge_length(m_boundary3D[j]));
    }
}

//destructor
LPFB_NLP::~LPFB_NLP()
{}

// returns the size of the problem
bool LPFB_NLP::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                            Index& nnz_h_lag, TNLP::IndexStyleEnum& index_style)
{
    // variables
    n = m_boundary3D.size();
    
    // Constraints
    m = 0;
    // Angle constraints
    m += m_boundaryCoincidentPair.size();
    // Position constraints
    m += m_coincidentPair.size() * 2;
    
    // Jacobian
    nnz_jac_g = 0;
    // Angle constraints
    for (size_t boundary_no = 0; boundary_no < m_boundaryCoincidentPair.size(); boundary_no++) {
        int startingIdx = m_boundaryCoincidentPair[boundary_no].first;
        int endingIdx = m_boundaryCoincidentPair[boundary_no].second;

        nnz_jac_g += endingIdx - startingIdx;
    }
    // Position constraints
    for (size_t pairNo = 0; pairNo < m_coincidentPair.size(); pairNo++) {
        int endingIdx = m_coincidentPair[pairNo].second;
        
        nnz_jac_g += 2 * endingIdx;
    }
    
    // the hessian is dense
//    nnz_h_lag = n*n;
    
    // use the C style indexing (0-based)
    index_style = TNLP::C_STYLE;
    
    return true;
}

// returns the variable bounds
bool LPFB_NLP::get_bounds_info(Index n, Number* x_l, Number* x_u,
                                Index m, Number* g_l, Number* g_u)
{
    // lower/upper bounds
    for (Index i = 0; i < n; i++) {
        x_l[i] = 0;
        x_u[i] = 2 * M_PI;
    }
    
    // All constraints should be 0
    for (Index i = 0; i < m; i++) {
        g_l[i] = g_u[i] = 0;
    }
    
    return true;
}

// returns the initial point for the problem
bool LPFB_NLP::get_starting_point(Index n, bool init_x, Number* x,
                                   bool init_z, Number* z_L, Number* z_U,
                                   Index m, bool init_lambda,
                                   Number* lambda)
{
    // initialize to the given starting point
    for (int i = 0; i < n; i++) {
        x[i] = m_innerAngle3D[i];
    }
    
    return true;
}

// returns the value of the objective function
bool LPFB_NLP::eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
{
    double sumAngleDiff = 0;
    for (int i = 0; i < n; i++) {
        double angleDiff = x[i] - m_innerAngle3D[i];
        sumAngleDiff += 0.5 * angleDiff * angleDiff;
    }
    
    obj_value = sumAngleDiff;
    return true;
}

// return the gradient of the objective function grad_{x} f(x)
bool LPFB_NLP::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
{
    for (int i = 0; i < n; i++) {
        grad_f[i] = x[i] - m_innerAngle3D[i];
    }
    return true;
}

// return the value of the constraints: g(x)
bool LPFB_NLP::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
{
    double sumAngle = (n - 2) * M_PI;
    
    double *posX = new double[n+1];
    double *posY = new double[n+1];
    
    double phi_k = 0;
    posX[0] = posY[0] = 0;
    for (Index k = 0; k < n; k++) {
        sumAngle -= x[k];
        phi_k += M_PI - x[k];
        posX[k+1] = posX[k] + m_edgeLengths[k] * cos(phi_k);
        posY[k+1] = posY[k] + m_edgeLengths[k] * sin(phi_k);
    }
    Index constraintIdx = 0;
    
    // Angle constraints, inner turning angle sums up to 2PI
    for (size_t boundary_no = 0; boundary_no < m_boundaryCoincidentPair.size(); boundary_no++) {
        int startingIdx = m_boundaryCoincidentPair[boundary_no].first;
        int endingIdx = m_boundaryCoincidentPair[boundary_no].second;
        
        double sumInnerAngle = 0;
        int nInnerAngles = 0;
        for (int i = startingIdx; i < endingIdx; i++, nInnerAngles++) {
            sumInnerAngle += x[i];
        }
        
        g[constraintIdx++] = (nInnerAngles - 2) * M_PI - sumInnerAngle;
    }

    // Position constraints
    for (size_t pairNo = 0; pairNo < m_coincidentPair.size(); pairNo++) {
        int startingIdx = m_coincidentPair[pairNo].first;
        int endingIdx = m_coincidentPair[pairNo].second;
        
        g[constraintIdx++] = posX[endingIdx] - posX[startingIdx];
        g[constraintIdx++] = posY[endingIdx] - posY[startingIdx];
    }
    delete[] posX;
    delete[] posY;
    
    return true;
}

// return the structure or values of the jacobian
bool LPFB_NLP::eval_jac_g(Index n, const Number* x, bool new_x,
                           Index m, Index nele_jac, Index* iRow, Index *jCol,
                           Number* values)
{
    if (values == NULL) {
        // return the structure of the jacobian
        // this particular jacobian is dense
        
        // Angle constraints
        Index idx = 0;
        Index constraintIdx = 0;

        for (size_t boundary_no = 0; boundary_no < m_boundaryCoincidentPair.size(); boundary_no++, constraintIdx++) {
            int startingIdx = m_boundaryCoincidentPair[boundary_no].first;
            int endingIdx = m_boundaryCoincidentPair[boundary_no].second;

            for (Index j = startingIdx; j < endingIdx; j++, idx++) {
                iRow[idx] = constraintIdx;
                jCol[idx] = j;
            }
        }
        
        // Position constraints
        for (size_t pairNo = 0; pairNo < m_coincidentPair.size(); pairNo++) {
            Index endingIdx = m_coincidentPair[pairNo].second;
            for (Index j = 0; j < endingIdx; j++, idx++) {
                // X pos
                iRow[idx] = constraintIdx;
                jCol[idx] = j;
                // Y pos
                iRow[idx+endingIdx] = constraintIdx+1;
                jCol[idx+endingIdx] = j;
            }
            constraintIdx += 2;
            idx += endingIdx;
        }
        
        assert(idx == nele_jac);
    }
    else {
        // return the values of the jacobian of the constraints
        double *dPosX = new double[n+1];
        double *dPosY = new double[n+1];

        double phi_k = 0;
        dPosX[0] = dPosY[0] = 0;
        for (Index k = 0; k < n; k++) {
            phi_k += M_PI - x[k];
            dPosX[k+1] = dPosX[k] + m_edgeLengths[k] * sin(phi_k);
            dPosY[k+1] = dPosY[k] + m_edgeLengths[k] * -cos(phi_k);
        }

        // Angle constraints
        Index idx = 0;
        
        for (size_t boundary_no = 0; boundary_no < m_boundaryCoincidentPair.size(); boundary_no++) {
            int startingIdx = m_boundaryCoincidentPair[boundary_no].first;
            int endingIdx = m_boundaryCoincidentPair[boundary_no].second;
            
            for (Index j = startingIdx; j < endingIdx; j++) {
                values[idx++] = -1;
            }
        }

        // Position constraints
        for (size_t pairNo = 0; pairNo < m_coincidentPair.size(); pairNo++) {
            Index startingIdx = m_coincidentPair[pairNo].first;
            Index endingIdx = m_coincidentPair[pairNo].second;

            for (Index k = 0; k < endingIdx; k++, idx++) {
                if (k < startingIdx) {
                    values[idx] = dPosX[endingIdx] - dPosX[startingIdx];
                    values[idx+endingIdx] = dPosY[endingIdx] - dPosY[startingIdx];
                } else {
                    values[idx] = dPosX[endingIdx] - dPosX[k];
                    values[idx+endingIdx] = dPosY[endingIdx] - dPosY[k];
                }
            }
            idx += endingIdx;
        }
        delete[] dPosX;
        delete[] dPosY;
    }
    
    return true;
}

//return the structure or values of the hessian
/*
 There are some problem for this function. I do calculate the answer correctly, but the values seems to be modified after this function and became incorrect. I now use built-in hessian_approximation to replace this function.
 This function is NOT updated for inner virtual cut.
 */
//bool LPFB_NLP::eval_h(Index n, const Number* x, bool new_x,
//                       Number obj_factor, Index m, const Number* lambda,
//                       bool new_lambda, Index nele_hess, Index* iRow,
//                       Index* jCol, Number* values)
//{
//    if (values == NULL) {
//        // the hessian for this problem is dense
//        Index idx=0;
//        for (Index i = 0; i < n; i++) {
//            for (Index j = 0; j < n; j++) {
//                iRow[idx] = i;
//                jCol[idx] = j;
//                idx++;
//            }
//        }
//        
//        assert(idx == nele_hess);
//    }
//    else {
//        // return the values.
//        // fill the objective portion which is identical matrix
//        int idx = 0;
//        for (Index i = 0; i < n; i++) {
//            for (Index j = 0; j < n; j++) {
//                if (i == j) {
//                    values[idx] = obj_factor;
//                } else {
//                    values[idx] = 0;
//                }
//                idx++;
//            }
//        }
//        
//        // 1st constraint is 0 after differentiation
//        // add the portion for the 2nd, 3rd constraints
//        double *dPosX = new double[n+1];
//        double *dPosY = new double[n+1];
//        
//        double phi_k = 0;
//        dPosX[0] = dPosY[0] = 0;
//        for (Index k = 0; k < n; k++) {
//            phi_k += M_PI - x[k];
//            dPosX[k+1] = dPosX[k] + m_edgeLengths[k] * -cos(phi_k);
//            dPosY[k+1] = dPosY[k] + m_edgeLengths[k] * -sin(phi_k);
//        }
//        
//        idx = 0;
//        for (Index i = 0; i < n; i++) {
//            for (Index j = 0; j < n; j++) {
//                // First constaints equals 0 after differentition
//
//                // 2nd & 3rd constraints, closed loop position
//                values[idx] += lambda[1] * (dPosX[n] - dPosX[max(i,j)]);
//                values[idx] += lambda[2] * (dPosY[n] - dPosY[max(i,j)]);
//                idx++;
//            }
//        }
//        delete[] dPosX;
//        delete[] dPosY;
//    }
//    
//    return true;
//}

void LPFB_NLP::finalize_solution(SolverReturn status,
                                  Index n, const Number* x, const Number* z_L, const Number* z_U,
                                  Index m, const Number* g, const Number* lambda,
                                  Number obj_value,
                                  const IpoptData* ip_data,
                                  IpoptCalculatedQuantities* ip_cq)
{
    std::cout << std::endl << std::endl << "Objective value" << std::endl;
    std::cout << "f(x*) = " << obj_value << std::endl;
    
    std::cout << std::endl << "Final value of the constraints:" << std::endl;
    for (Index i=0; i<m ;i++) {
        std::cout << "g(" << i << ") = " << g[i] << std::endl;
    }
    
    double *phi_k = new double[n];
    double *posX = new double[n+1];
    double *posY = new double[n+1];
    posX[0] = posY[0] = 0;
    for (Index k = 0; k < n; k++) {
        if (k == 0) {
            phi_k[k] = M_PI - x[k];
        } else {
            phi_k[k] = phi_k[k-1] + M_PI - x[k];
        }
        posX[k+1] = posX[k] + m_edgeLengths[k] * cos(phi_k[k]);
        posY[k+1] = posY[k] + m_edgeLengths[k] * sin(phi_k[k]);
        
//        printf("%d : %d | %.3lf | %.3lf\n", k, m_mesh->from_vertex_handle(m_boundary3D[k]).idx(), x[k], x[k]-m_innerAngle3D[k]);
    }
    
    OpenMesh::MPropHandleT< std::set<EdgeHandle> > seamsHandle = PlushPatternGenerator::getSeamsHandle(m_mesh);
    std::set<EdgeHandle> &seams = m_mesh->property(seamsHandle);
    for (Index k = 0; k < n; k++) {
        // If k == 0, it will be placed at (0,0,0)
        seams.insert(m_mesh->edge_handle(m_boundary3D[k]));
        VertexHandle v = m_mesh->from_vertex_handle(m_boundary3D[k]);
        if (m_mesh->is_boundary(v)) {
            m_boundaryPosition->emplace(v, OpenMesh::Vec3d(posX[k], posY[k], 0));
        }
    }
    
    delete[] phi_k;
    delete[] posX;
    delete[] posY;
}
