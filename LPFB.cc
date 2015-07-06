#include "OpenMesh_Boost_Wrapper.hh"
#include "LPFB.hh"
#include "PlushPatternGenerator.hh"

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>

#include <cassert>
#include <iostream>

#include <queue>

// constructor
LPFB_NLP::LPFB_NLP(FilteredTriMesh &mesh, std::map<HalfedgeHandle, OpenMesh::Vec3d> &boundaryPosition) : m_mesh(mesh), m_boundaryPosition(boundaryPosition)
{
    std::vector< std::vector<HalfedgeHandle> > boundaries = PlushPatternGenerator::getBoundaryOfOpenedMesh(m_mesh, true);
    
    assert(!boundaries.empty() && "No cut is found on given sub mesh.");
    
    // We use the boundary with most vertices number as base boundary. Other boundaries are treated as inner boundaries.
    // The target set. Virtual cut searching will stop when reaching target set.
    std::vector<HalfedgeHandle> *baseBoundary = &boundaries[0];
    for (size_t j = 1; j < boundaries.size(); j++) {
        if (boundaries[j].size() > baseBoundary->size()) {
            baseBoundary = &boundaries[j];
        }
    }

    m_boundary3D.insert(m_boundary3D.begin(), baseBoundary->begin(), baseBoundary->end());

    // If there are more than one boundaries, we need to generate virtual cuts between these boundaries
    // and link these boundaries together so that LPFB still works
    std::set<VertexHandle> virtual_cut_v;
    if (boundaries.size() > 1) {
        // Because of unknown reason, the multiple source version of Dijkastra doesn't work.
        // We have to constrcut a new graph with an additional starting point to achieve it.
        typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
        boost::property<boost::vertex_owner_t, VertexHandle>,
        boost::property<boost::edge_weight_t, double>
        > Graph;
        Graph multi_source_graph;
        
        std::map<VertexHandle, unsigned long> triMesh_to_graph_vertex_mapping;
        for (const VertexHandle vh : m_mesh.vertices()) {
            Graph::vertex_descriptor vid = boost::add_vertex(vh, multi_source_graph);
            triMesh_to_graph_vertex_mapping.emplace(vh, vid);
        }
        for (const EdgeHandle eh : m_mesh.edges()) {
            HalfedgeHandle heh = m_mesh.halfedge_handle(eh, 0);
            VertexHandle v1 = m_mesh.from_vertex_handle(heh);
            VertexHandle v2 = m_mesh.to_vertex_handle(heh);
            
            Graph::vertex_descriptor vid1 = triMesh_to_graph_vertex_mapping[v1];
            Graph::vertex_descriptor vid2 = triMesh_to_graph_vertex_mapping[v2];
            boost::add_edge(vid1, vid2, m_mesh.calc_edge_length(eh), multi_source_graph);
        }

        // For inner boundaries, calculate virtual cut and insert them into appropriate position of m_boundary3D
        for (size_t i = 0; i < boundaries.size(); i++) {
            std::vector<HalfedgeHandle> &sourceBoundary = boundaries[i];
            if (sourceBoundary == *baseBoundary) {
                continue;
            }

            std::set<VertexHandle> sourceBoundaryVertices;
            for (HalfedgeHandle heh : sourceBoundary) {
                sourceBoundaryVertices.insert(m_mesh.from_vertex_handle(heh));
            }

            // we solve multi-source shoretest path problem by manually construct a new graph
            // with new additional vertex as starting point and inserting 0-weight edges
            // from starting point to other point on source boundary
            Graph::vertex_descriptor startingV = boost::add_vertex(multi_source_graph);
            for (VertexHandle vh : sourceBoundaryVertices) {
                assert(triMesh_to_graph_vertex_mapping.find(vh) != triMesh_to_graph_vertex_mapping.end());
                Graph::vertex_descriptor vid = triMesh_to_graph_vertex_mapping[vh];
                boost::add_edge(startingV, vid, 0, multi_source_graph);
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
                VertexHandle v = m_mesh.to_vertex_handle(baseBoundary->at(j));
                Graph::vertex_descriptor vid = triMesh_to_graph_vertex_mapping[v];
                if (distance_map[vid] < minDistance) {
                    minDistance = distance_map[vid];
                    minDstV = vid;
                }
            }
            assert(minDistance > 0 && "It's pretty likely that your boundary is wrong.");
            
            // Back tracing
            std::vector<HalfedgeHandle> virtualCut;
            std::vector<HalfedgeHandle> virtualCutReverse;
            Graph::vertex_descriptor currentV = minDstV, nextV = predecessor_map[minDstV];
            while (nextV != startingV) {
                VertexHandle v_from = get(boost::vertex_owner, multi_source_graph, currentV);
                VertexHandle v_to = get(boost::vertex_owner, multi_source_graph, nextV);
                
                // Make sure multi-source dijkstra is correct.
                // i.e., the path doesn't pass through any vertex on source boundary except source vertex
                assert(std::find(sourceBoundaryVertices.begin(), sourceBoundaryVertices.end(), v_from) == sourceBoundaryVertices.end());

                HalfedgeHandle heh;
                bool edgeExist = PlushPatternGenerator::getHalfedge(m_mesh, heh, v_from, v_to);
                assert(edgeExist);
                
                virtualCut.push_back(heh);
                virtualCutReverse.push_back(m_mesh.opposite_halfedge_handle(heh));
                virtual_cut_v.insert(v_from);
                virtual_cut_v.insert(v_to);
                
                currentV = nextV;
                nextV = predecessor_map[nextV];
            }
            std::reverse(virtualCutReverse.begin(), virtualCutReverse.end());
            
            // Find the insertion position (the index of starting point in source boundary)
            int insertionPosition = -1;
            for (size_t j = 0; j < sourceBoundary.size(); j++) {
                if (triMesh_to_graph_vertex_mapping[m_mesh.from_vertex_handle(sourceBoundary[j])] == currentV) {
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
                if (m_mesh.from_vertex_handle(heh) == m_mesh.from_vertex_handle(virtualCut[0])) {
                    he_it = m_boundary3D.insert(he_it, virtualCut.begin(), virtualCut.end());
                    he_it = m_boundary3D.insert(he_it+virtualCut.size(), sourceBoundary.begin(), sourceBoundary.end());
                    he_it = m_boundary3D.insert(he_it+sourceBoundary.size(), virtualCutReverse.begin(), virtualCutReverse.end());
                    
                    // Check if this boundary is well connected together
                    for (size_t edgeNo = 0; edgeNo < m_boundary3D.size(); edgeNo++) {
                        assert(m_mesh.to_vertex_handle(m_boundary3D[edgeNo])
                            == m_mesh.from_vertex_handle(m_boundary3D[(edgeNo+1)%m_boundary3D.size()]));
                    }
                    break;
                }
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
        VertexHandle v = m_mesh.from_vertex_handle(m_boundary3D[position]);
        
        // If this vertex is visited before, there is a loop. And we need to add it to coincident pairs
        // Only vertex on virtual cut may be coincident
        if (virtual_cut_v.find(v) != virtual_cut_v.end()
        &&  vertexPosition.find(v) != vertexPosition.end()) {
            int prevPosition = vertexPosition[v];
            m_coincidentPair.push_back(std::make_pair(prevPosition, position));
//            if (!cycle) {
//                m_boundaryCoincidentPair.push_back(std::make_pair(prevPosition, position));
//            }
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
        double sumInnerAngle = m_mesh.getSumInnerAngle(m_boundary3D[(j-1+n) % n], m_boundary3D[j]);
        m_innerAngle3D.push_back(sumInnerAngle);
        
        m_edgeLengths.push_back(m_mesh.calc_edge_length(m_boundary3D[j]));
    }
}

void LPFB_NLP::newtons_method() {
    int n, m, nnz_jac_g, nnz_h_lag;
    get_nlp_info(n, m, nnz_jac_g);

    double *cos_phi = new double[n];
    double *sin_phi = new double[n];
    Eigen::VectorXd X(n);
    Eigen::VectorXd g(m);
    Eigen::VectorXd grad_f(n);
    Eigen::MatrixXd grad_g(m, n);
    
    get_starting_point(n, X);
    get_phi(n, X, cos_phi, sin_phi);
    eval_grad_f(n, X, grad_f);
    eval_g(n, X, g, cos_phi, sin_phi);
    eval_jac_g(n, X, m, grad_g, cos_phi, sin_phi);

    Eigen::VectorXd lambda(m);
    lambda.setOnes();

    double derivative = 1e9;
    Eigen::VectorXd delta_theta(n);
    Eigen::VectorXd delta_lambda(m);
    delta_theta.setOnes();
    delta_lambda.setOnes();
    for (int iterations = 0; iterations < 20 && derivative > 1e-2 && delta_theta.norm() > 1e-2 && delta_lambda.norm() > 1e-2; iterations++) {
        get_delta(n, m,
                  X,
                  g,
                  grad_f,
                  grad_g,
                  lambda, delta_theta, delta_lambda);
        
//        printf("X:");
//        for (int i = 0; i < n; i++) {
//            printf("%.3lf, ", X[i]);
//        }
//        printf("\n");
//        printf("delta_theta:");
//        for (int i = 0; i < n; i++) {
//            printf("%.3lf, ", delta_theta[i]);
//        }
//        printf("\n");
//        printf("delta_lambda:");
//        for (int i = 0; i < m; i++) {
//            printf("%.3lf, ", delta_lambda[i]);
//        }
//        printf("\n\n");

        X += delta_theta;
        lambda += delta_lambda;
        for (int i = 0; i < n; i++) {
            X[i] = min(max(X[i], 0.0), 2*M_PI);
        }
        
        get_phi(n, X, cos_phi, sin_phi);
        eval_grad_f(n, X, grad_f);
        eval_g(n, X, g, cos_phi, sin_phi);
        eval_jac_g(n, X, m, grad_g, cos_phi, sin_phi);
        
        double component_sqr = 0;
        for (int p = 0; p < m; p++) {
            grad_g.row(p) *= lambda[p];
        }
        for (int i = 0; i < n; i++) {
            component_sqr += pow(grad_f[i] + grad_g.col(i).sum(), 2);
        }
        // derivative for lambda_theta, lambda_0x, lambda_0y, ...
        for (int p = 0; p < m; p++) {
            component_sqr += pow(g[p], 2);
        }
        
        derivative = sqrt(component_sqr);
    }
    
    finalize_solution(n, X, m, lambda, cos_phi, sin_phi);
    
    delete[] cos_phi;
    delete[] sin_phi;
}

bool LPFB_NLP::get_delta(int n, int m,
                        Eigen::VectorXd &X,
                        Eigen::VectorXd &g,
                        Eigen::VectorXd &grad_f,
                        Eigen::MatrixXd &grad_g,
                        Eigen::VectorXd &lambda,
                        Eigen::VectorXd &delta_theta,
                        Eigen::VectorXd &delta_lambda) {
    // B_lambda
    Eigen::VectorXd B_lambda = g * -1;

    // B_theta
    Eigen::VectorXd B_theta(n);
    for (int i = 0; i < n; i++) {
        B_theta[i] = grad_f[i];
        for (int p = 0; p < m; p++) {
            B_theta[i] += grad_g(p, i) * lambda[p];
        }
    }
    B_theta *= -1;
    
    Eigen::MatrixXd big_lambda = grad_g;
    Eigen::MatrixXd big_lambda_x_big_lambda = big_lambda * big_lambda.transpose();
    Eigen::VectorXd lambda_rhs = big_lambda * B_theta - B_lambda;
    delta_lambda = big_lambda_x_big_lambda.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(lambda_rhs);
    delta_theta = B_theta - big_lambda.transpose() * delta_lambda;

    return true;
}

bool LPFB_NLP::get_phi(int n, Eigen::VectorXd &X, double *cos_phi, double *sin_phi) {
    double phi_k = 0;
    for (int k = 0; k < n; k++) {
        phi_k += M_PI - X[k];
        cos_phi[k] = cos(phi_k);
        sin_phi[k] = sin(phi_k);
    }
    
    return true;
}

// returns the size of the problem
bool LPFB_NLP::get_nlp_info(int& n, int& m, int& nnz_jac_g)
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
    return true;
}

// returns the variable bounds
//bool LPFB_NLP::get_bounds_info(int n, Number* x_l, Number* x_u,
//                                int m, Number* g_l, Number* g_u)
//{
//    // lower/upper bounds
//    for (int i = 0; i < n; i++) {
//        x_l[i] = 0;
//        x_u[i] = 2 * M_PI;
//    }
//    
//    // All constraints should be 0
//    for (int i = 0; i < m; i++) {
//        g_l[i] = g_u[i] = 0;
//    }
//    
//    return true;
//}

// returns the initial point for the problem
bool LPFB_NLP::get_starting_point(int n, Eigen::VectorXd &X)
{
    // initialize to the given starting point
    for (int i = 0; i < n; i++) {
        X[i] = m_innerAngle3D[i];
    }
    
    return true;
}

// returns the value of the objective function
bool LPFB_NLP::eval_f(int n, Eigen::VectorXd &X, double &obj_value)
{
    double sumAngleDiff = 0;
    for (int i = 0; i < n; i++) {
        double angleDiff = X[i] - m_innerAngle3D[i];
        sumAngleDiff += 0.5 * angleDiff * angleDiff;
    }
    
    obj_value = sumAngleDiff;
    return true;
}

// return the gradient of the objective function grad_{x} f(x)
bool LPFB_NLP::eval_grad_f(int n, Eigen::VectorXd &X, Eigen::VectorXd &grad_f)
{
    for (int i = 0; i < n; i++) {
        grad_f[i] = X[i] - m_innerAngle3D[i];
    }
    return true;
}

// return the value of the constraints: g(x)
bool LPFB_NLP::eval_g(int n, Eigen::VectorXd &X, Eigen::VectorXd &g, double *cos_phi, double *sin_phi)
{
    double *posX = new double[n+1];
    double *posY = new double[n+1];
    posX[0] = posY[0] = 0;
    for (int k = 0; k < n; k++) {
        posX[k+1] = posX[k] + m_edgeLengths[k] * cos_phi[k];
        posY[k+1] = posY[k] + m_edgeLengths[k] * sin_phi[k];
    }
    int constraintIdx = 0;
    
    // Angle constraints, inner turning angle sums up to 2PI
    for (size_t boundary_no = 0; boundary_no < m_boundaryCoincidentPair.size(); boundary_no++) {
        int startingIdx = m_boundaryCoincidentPair[boundary_no].first;
        int endingIdx = m_boundaryCoincidentPair[boundary_no].second;
        
        double sumInnerAngle = 0;
        int nInnerAngles = 0;
        for (int i = startingIdx; i < endingIdx; i++, nInnerAngles++) {
            sumInnerAngle += X[i];
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
bool LPFB_NLP::eval_jac_g(int n, Eigen::VectorXd &X,
                          int m, Eigen::MatrixXd &values,
                          double *cos_phi, double *sin_phi)
{
    values.setZero();

    // return the values of the jacobian of the constraints
    double *dPosX = new double[n+1];
    double *dPosY = new double[n+1];

    dPosX[0] = dPosY[0] = 0;
    for (int k = 0; k < n; k++) {
        dPosX[k+1] = dPosX[k] + m_edgeLengths[k] * sin_phi[k];
        dPosY[k+1] = dPosY[k] + m_edgeLengths[k] * -cos_phi[k];
    }

    // Angle constraints
    int constraintIdx = 0;
    for (size_t boundary_no = 0; boundary_no < m_boundaryCoincidentPair.size(); boundary_no++) {
        int startingIdx = m_boundaryCoincidentPair[boundary_no].first;
        int endingIdx = m_boundaryCoincidentPair[boundary_no].second;
        
        for (int j = startingIdx; j < endingIdx; j++) {
            values(constraintIdx, j) = -1;
        }
        constraintIdx++;
    }

    // Position constraints
    for (size_t pairNo = 0; pairNo < m_coincidentPair.size(); pairNo++) {
        int startingIdx = m_coincidentPair[pairNo].first;
        int endingIdx = m_coincidentPair[pairNo].second;

        for (int k = 0; k < endingIdx; k++) {
            if (k < startingIdx) {
                values(constraintIdx, k) = dPosX[endingIdx] - dPosX[startingIdx];
                values(constraintIdx+1, k) = dPosY[endingIdx] - dPosY[startingIdx];
            } else {
                values(constraintIdx, k) = dPosX[endingIdx] - dPosX[k];
                values(constraintIdx+1, k) = dPosY[endingIdx] - dPosY[k];
            }
        }
        constraintIdx += 2;
    }
    delete[] dPosX;
    delete[] dPosY;
    
    return true;
}

//return the structure or values of the hessian
/*
 There are some problem for this function. I do calculate the answer correctly, but the values seems to be modified after this function and became incorrect. I now use built-in hessian_approximation to replace this function.
 This function is NOT updated for inner virtual cut.
 */
//bool LPFB_NLP::eval_h(int n, const Number* x, bool new_x,
//                       Number obj_factor, int m, const Number* lambda,
//                       bool new_lambda, int nele_hess, int* iRow,
//                       int* jCol, Number* values)
//{
//    if (values == NULL) {
//        // the hessian for this problem is dense
//        int idx=0;
//        for (int i = 0; i < n; i++) {
//            for (int j = 0; j < n; j++) {
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
//        for (int i = 0; i < n; i++) {
//            for (int j = 0; j < n; j++) {
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
//        for (int k = 0; k < n; k++) {
//            phi_k += M_PI - x[k];
//            dPosX[k+1] = dPosX[k] + m_edgeLengths[k] * -cos(phi_k);
//            dPosY[k+1] = dPosY[k] + m_edgeLengths[k] * -sin(phi_k);
//        }
//        
//        idx = 0;
//        for (int i = 0; i < n; i++) {
//            for (int j = 0; j < n; j++) {
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

void LPFB_NLP::finalize_solution(int n, Eigen::VectorXd &X,
                                 int m, Eigen::VectorXd &lambda,
                                 double *cos_phi, double *sin_phi) {
    double obj_value;
    eval_f(n, X, obj_value);
    std::cout << std::endl << std::endl << "Objective value" << std::endl;
    std::cout << "f(x*) = " << obj_value << std::endl;

    Eigen::VectorXd g(m);
    eval_g(n, X, g, cos_phi, sin_phi);
    std::cout << std::endl << "Final value of the constraints:" << std::endl;
    for (int i=0; i<m ;i++) {
        std::cout << "g(" << i << ") = " << g[i] << std::endl;
    }
    
    double *posX = new double[n+1];
    double *posY = new double[n+1];
    posX[0] = posY[0] = 0;
    for (int k = 0; k < n; k++) {
        posX[k+1] = posX[k] + m_edgeLengths[k] * cos_phi[k];
        posY[k+1] = posY[k] + m_edgeLengths[k] * sin_phi[k];
        
//        printf("%d : %d | %.3lf | %.3lf\n", k, m_mesh.from_vertex_handle(m_boundary3D[k]).idx(), x[k], x[k]-m_innerAngle3D[k]);
    }
    
    for (int k = 0; k < n; k++) {
        // If k == 0, it will be placed at (0,0,0)
        // Ignore virtual cuts, they are not real boundary
        // Also, since we use the opposite of boundary halfedge which is not boundary, we need to use edge handle to test both side.
        if (m_mesh.is_boundary(m_mesh.edge_handle(m_boundary3D[k]))) {
            m_boundaryPosition.emplace(m_boundary3D[k], OpenMesh::Vec3d(posX[k], posY[k], 0));
        }
    }
    delete[] posX;
    delete[] posY;
}
