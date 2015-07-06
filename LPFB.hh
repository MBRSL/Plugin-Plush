#ifndef LPFB_HH
#define LPFB_HH

#include "Common.hh"
#include "FilteredTriMesh.hh"

#include <Eigen/Dense>
#include <Eigen/SVD>

class LPFB_NLP
{
public:
    LPFB_NLP(FilteredTriMesh &mesh, std::map<HalfedgeHandle, OpenMesh::Vec3d> &boundaryPosition);
    
    void newtons_method();
    
    /** Method to return some info about the nlp */
    bool get_nlp_info(int& n, int& m, int& nnz_jac_g);

    /** Method to return the starting point for the algorithm */
    bool get_starting_point(int n, Eigen::VectorXd &X);
    
    /** Method to return the objective value */
    bool eval_f(int n, Eigen::VectorXd &X, double &obj_value);
    
    /** Method to return the gradient of the objective */
    bool eval_grad_f(int n, Eigen::VectorXd &X, Eigen::VectorXd &grad_f);
    
    /** Method to return the constraint residuals */
    bool eval_g(int n, Eigen::VectorXd &X, Eigen::VectorXd &g,
                double *cos_phi, double *sin_phi);
    
    /** Method to return:
     *   1) The structure of the jacobian (if "values" is NULL)
     *   2) The values of the jacobian (if "values" is not NULL)
     */
    bool eval_jac_g(int n, Eigen::VectorXd &X,
                    int m, Eigen::MatrixXd &values,
                    double *cos_phi, double *sin_phi);
    
    bool get_delta(int n, int m,
                   Eigen::VectorXd &X,
                   Eigen::VectorXd &g,
                   Eigen::VectorXd &grad_f,
                   Eigen::MatrixXd &grad_g,
                   Eigen::VectorXd &lambda,
                   Eigen::VectorXd &delta_theta,
                   Eigen::VectorXd &delta_lambda);
    
    bool get_phi(int n, Eigen::VectorXd &X, double *cos_phi, double *sin_phi);
    //@}
    
    /** @name Solution Methods */
    //@{
    /** This method is called when the algorithm is complete so the TNLP can store/write the solution */
    void finalize_solution(int n, Eigen::VectorXd &X,
                           int m, Eigen::VectorXd &lambda,
                           double *cos_phi, double *sin_phi);
    //@}
    
private:
    FilteredTriMesh &m_mesh;
    std::map<HalfedgeHandle, OpenMesh::Vec3d> &m_boundaryPosition;
    std::vector<double> m_innerAngle3D;
    std::vector<double> m_edgeLengths;
    std::vector<HalfedgeHandle> m_boundary3D;
    /// Same as m_coincidentPair except that it only contains pairs of boundary points. No cuts points.
    /// This is used for calculating angle constraints
    std::vector< std::pair<int, int> > m_boundaryCoincidentPair;
    /// Records pairs of coincident points. Including points on cut.
    std::vector< std::pair<int, int> > m_coincidentPair;
};


#endif
