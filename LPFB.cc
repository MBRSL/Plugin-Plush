#include "LPFB.hh"
#include "PlushPatternGenerator.hh"

#include <cassert>
#include <iostream>

using namespace Ipopt;

// constructor
LPFB_NLP::LPFB_NLP(TriMesh *mesh) : m_mesh(mesh)
{
    std::vector< std::vector<HalfedgeHandle> > boundaries;
    PlushPatternGenerator::getBoundaryOfOpenedMesh(boundaries, m_mesh, true);
    m_boundary3D = boundaries[0];
    int n = m_boundary3D.size();

    // Calculate inner angles & edges length for later use.
    for (int j = 0; j < n; j++) {
        double sumInnerAngle = PlushPatternGenerator::getSumInnerAngle(m_mesh, m_boundary3D[(j-1+n) % n], m_boundary3D[j]);
        m_innerAngle3D.push_back(sumInnerAngle);
        
        TriMesh::Point v1 = m_mesh->point(m_mesh->from_vertex_handle(m_boundary3D[j]));
        TriMesh::Point v2 = m_mesh->point(m_mesh->to_vertex_handle(m_boundary3D[j]));
        double length = (v1-v2).norm();
        m_edgeLengths.push_back(length);
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
    m = 3;
    
    // jacobian is dense
    nnz_jac_g = 3 * n;
    
    // the hessian is dense
    nnz_h_lag = n*n;
    
    // use the C style indexing (0-based)
    index_style = TNLP::C_STYLE;
    
    return true;
}

// returns the variable bounds
bool LPFB_NLP::get_bounds_info(Index n, Number* x_l, Number* x_u,
                                Index m, Number* g_l, Number* g_u)
{
    assert(m == 3);
    
    // lower/upper bounds
    // We set the range little wider than [0, 2*M_PI]
    for (Index i = 0; i < n; i++) {
        x_l[i] = -1 * M_PI;
        x_u[i] = 3 * M_PI;
    }
    
    // All constraints should be 0
    g_l[0] = g_u[0] = g_l[1] = g_u[1] = g_l[2] = g_u[2] = 0;
    
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
        x[i] = min(2 * M_PI, m_innerAngle3D[i]);
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
    assert(m == 3);
    
    double loopAngleConstrain = (n - 2) * M_PI;
    
    double loopPositionConstrainsX = 0;
    double loopPositionConstrainsY = 0;
    double phi_i = 0;
    for (int i = 0; i < n; i++) {
        double angleDiff = x[i] - m_innerAngle3D[i];
        
        loopAngleConstrain -= x[i];
        
        phi_i += M_PI - x[i];
        loopPositionConstrainsX += m_edgeLengths[i] * cos(phi_i);
        loopPositionConstrainsY += m_edgeLengths[i] * sin(phi_i);
    }
    
    g[0] = loopAngleConstrain;
    g[1] = loopPositionConstrainsX;
    g[2] = loopPositionConstrainsY;
    
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
        
        // M constraints
        Index idx = 0;
        for (Index i = 0; i < m; i++) {
            // N variables
            for (Index j = 0; j < n; j++) {
                iRow[idx] = i;
                jCol[idx] = j;
                idx++;
            }
        }
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

        for (Index j = 0; j < n; j++) {
            // First constaints, inner turning angle sums up to 2PI
            values[j] = -1;

            // 2nd & 3rd constraints, closed loop position
            values[n + j] = dPosX[n] - dPosX[j];
            values[2*n + j] = dPosY[n] - dPosY[j];
        }
        delete[] dPosX;
        delete[] dPosY;
    }
    
    return true;
}

//return the structure or values of the hessian
/*
 There are some problem for this function. I do calculate the answer correctly, but the values seems to be modified after this function and became incorrect. I now use built-in hessian_approximation to replace this function.
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
//    std::cout << std::endl << std::endl << "Solution of the primal variables, x" << std::endl;
//    for (Index i=0; i<n; i++) {
//        std::cout << "x[" << i << "] = " << x[i] << std::endl;
//    }
//    
//    std::cout << std::endl << std::endl << "Solution of the bound multipliers, z_L and z_U" << std::endl;
//    for (Index i=0; i<n; i++) {
//        std::cout << "z_L[" << i << "] = " << z_L[i] << std::endl;
//    }
//    for (Index i=0; i<n; i++) {
//        std::cout << "z_U[" << i << "] = " << z_U[i] << std::endl;
//    }
//    
    std::cout << std::endl << std::endl << "Objective value" << std::endl;
    std::cout << "f(x*) = " << obj_value << std::endl;
    
    std::cout << std::endl << "Final value of the constraints:" << std::endl;
    for (Index i=0; i<m ;i++) {
        std::cout << "g(" << i << ") = " << g[i] << std::endl;
    }
    
    double phi_k = 0;
    double posX = 0, posY = 0;
    for (int k = 0; k <= n; k++) {
        if (k > 0) {
            phi_k += M_PI - x[k-1];
            
            double lk = m_edgeLengths[k-1];
            
            posX += lk * cos(phi_k);
            posY += lk * sin(phi_k);
        }
        // If k == 0, it will be placed at (0,0,0)
        VertexHandle v = m_mesh->from_vertex_handle(m_boundary3D[k % n]);
        TriMesh::Point &p = m_mesh->point(v);
        p[0] = posX;
        p[1] = posY;
        p[2] = 0;
    }
}
