#include "LPFB.hh"
#include "PlushPatternGenerator.hh"

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

    // We use the 1st boundary as base boundary. Other boundaries are treated as inner boundaries.
    // The target set. Virtual cut searching will stop when reaching target set.
    std::vector<HalfedgeHandle> &baseBoundary = boundaries[0];
    std::set<VertexHandle> targetBoudarySet;
    for (size_t j = 0; j < baseBoundary.size(); j++) {
        VertexHandle v = m_mesh->to_vertex_handle(baseBoundary[j]);
        targetBoudarySet.insert(v);
    }

    m_boundary3D.insert(m_boundary3D.begin(), baseBoundary.begin(), baseBoundary.end());

    // For inner boundaries, calculate virtual cut and insert them into appropriate position of m_boundary3D
    for (size_t i = 1; i < boundaries.size(); i++) {
        std::vector<HalfedgeHandle> &sourceBoundary = boundaries[i];

        // Calculate the virtual cut using BFS
        // Randomly choose one vertex as starting point
        std::set<VertexHandle> sourceBoundarySet;
        for (size_t j = 0; j < sourceBoundary.size(); j++) {
            VertexHandle v = m_mesh->to_vertex_handle(sourceBoundary[j]);
            sourceBoundarySet.insert(v);
        }
        
        std::queue<VertexHandle> queue;
        std::set<VertexHandle> visited;
        std::map<VertexHandle, VertexHandle> predecessor;
        
        VertexHandle currentV = m_mesh->from_vertex_handle(sourceBoundary[0]);
        queue.push(currentV);
        visited.insert(currentV);
        
        while (!queue.empty() && targetBoudarySet.find(currentV) == targetBoudarySet.end()) {
            for (TriMesh::VertexVertexIter vv_it = m_mesh->vv_iter(currentV); vv_it; vv_it++) {
                if (visited.find(*vv_it) == visited.end()
                    // Don't go through sourceBoundary itself
                &&  sourceBoundarySet.find(*vv_it) == sourceBoundarySet.end()) {
                    queue.push(*vv_it);
                    visited.insert(*vv_it);
                    
                    predecessor.emplace(*vv_it, currentV);
                }
            }
            currentV = queue.front();
            queue.pop();
        }
        visited.clear();

        // We got to one of the targets
        std::vector<HalfedgeHandle> virtualCut;
        std::vector<HalfedgeHandle> virtualCutReverse;
        while (predecessor.find(currentV) != predecessor.end()) {
            VertexHandle nextV = predecessor[currentV];
            HalfedgeHandle heh;
            PlushPatternGenerator::getHalfedge(m_mesh, heh, currentV, nextV);
            virtualCut.push_back(heh);
            virtualCutReverse.push_back(m_mesh->opposite_halfedge_handle(heh));
            currentV = nextV;
        }
        std::reverse(virtualCutReverse.begin(), virtualCutReverse.end());
        
        // Insert this boundary and virtual cut into overall boundary
        for (auto he_it = m_boundary3D.begin(); he_it != m_boundary3D.end(); he_it++) {
            HalfedgeHandle heh = *he_it;
            // Find the insertion position
            if (m_mesh->from_vertex_handle(heh) == m_mesh->from_vertex_handle(virtualCut[0])) {
                he_it = m_boundary3D.insert(he_it, virtualCut.begin(), virtualCut.end());
                he_it = m_boundary3D.insert(he_it+virtualCut.size(), sourceBoundary.begin(), sourceBoundary.end());
                he_it = m_boundary3D.insert(he_it+sourceBoundary.size(), virtualCutReverse.begin(), virtualCutReverse.end());
                break;
            }
        }
    }
    
    // Calculating coincide pairs
    // At least the base boundary should be connected
    m_coincidePair.push_back(std::make_pair(0, m_boundary3D.size()));
    // Other inner boundaries
    std::map<VertexHandle, int> vertexPosition;
    for (size_t position = 0; position < m_boundary3D.size(); position++) {
        VertexHandle v = m_mesh->from_vertex_handle(m_boundary3D[position]);
        
        // If this vertex is visited before, there is a loop. And we need to add it to coincide pairs
        if (vertexPosition.find(v) != vertexPosition.end()) {
            int prevPosition = vertexPosition[v];
            m_coincidePair.push_back(std::make_pair(prevPosition, position));
            
            // We need to delete it first before update it.
            vertexPosition.erase(v);
        }
        vertexPosition.emplace(v, position);
    }
    
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
    
    // constraints
    m = 1 + m_coincidePair.size() * 2;
    
    // jacobian
    nnz_jac_g = n;
    for (size_t pairNo = 0; pairNo < m_coincidePair.size(); pairNo++) {
        int endingIdx = m_coincidePair[pairNo].second;

        nnz_jac_g += 2*endingIdx;
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
    // We set the range little wider than [0, 2*M_PI]
    for (Index i = 0; i < n; i++) {
        x_l[i] = -1 * M_PI;
        x_u[i] = 3 * M_PI;
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
    g[0] = sumAngle;

    for (size_t pairNo = 0; pairNo < m_coincidePair.size(); pairNo++) {
        int startingIdx = m_coincidePair[pairNo].first;
        int endingIdx = m_coincidePair[pairNo].second;
        
        double v1 = posX[endingIdx];
        double v2 = posX[startingIdx];
        g[1 + 2*pairNo] = posX[endingIdx] - posX[startingIdx];
        g[2 + 2*pairNo] = posY[endingIdx] - posY[startingIdx];
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
        
        // First constraints
        // N variables
        Index idx = 0;
        for (Index j = 0; j < n; j++) {
            iRow[idx] = 0;
            jCol[idx] = j;
            idx++;
        }
        
        for (size_t pairNo = 0; pairNo < m_coincidePair.size(); pairNo++) {
            Index endingIdx = m_coincidePair[pairNo].second;
            for (Index j = 0; j < endingIdx; j++) {
                iRow[idx] = 1 + 2*pairNo;
                jCol[idx] = j;
                idx++;
            }
            for (Index j = 0; j < endingIdx; j++) {
                iRow[idx] = 2 + 2*pairNo;
                jCol[idx] = j;
                idx++;
            }
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

        // First constaints, inner turning angle sums up to 2PI
        for (Index k = 0; k < n; k++) {
            values[k] = -1;
        }
        Index idx = n;

        // 2nd & 3rd constraints, closed loop position
        for (size_t pairNo = 0; pairNo < m_coincidePair.size(); pairNo++) {
            Index startingIdx = m_coincidePair[pairNo].first;
            Index endingIdx = m_coincidePair[pairNo].second;

            for (Index k = 0; k < endingIdx; k++, idx++) {
                if (k < startingIdx) {
                    values[idx] = dPosX[endingIdx] - dPosX[startingIdx];
                } else {
                    values[idx] = dPosX[endingIdx] - dPosX[k];
                }
            }
            for (Index k = 0; k < endingIdx; k++, idx++) {
                if (k < startingIdx) {
                    values[idx] = dPosY[endingIdx] - dPosY[startingIdx];
                } else {
                    values[idx] = dPosY[endingIdx] - dPosY[k];
                }
            }
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
        
        printf("%d : %d | %.3lf\n", k, m_mesh->from_vertex_handle(m_boundary3D[k]).idx(), x[k]);
    }
    
    for (Index k = 0; k < n; k++) {
        // If k == 0, it will be placed at (0,0,0)
        VertexHandle v = m_mesh->from_vertex_handle(m_boundary3D[k]);
        if (m_mesh->is_boundary(v)) {
            m_boundaryPosition->emplace(v, OpenMesh::Vec3d(posX[k], posY[k], 0));
        }
    }
}