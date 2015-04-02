#ifndef LPFB_HH
#define LPFB_HH

#include "IpTNLP.hpp"
#include "Common.hh"

class LPFB_NLP : public Ipopt::TNLP
{
public:
    LPFB_NLP(TriMesh *mesh);
    virtual ~LPFB_NLP();
    
    /**@name Overloaded from TNLP */
    //@{
    /** Method to return some info about the nlp */
    virtual bool get_nlp_info(Ipopt::Index& n, Ipopt::Index& m, Ipopt::Index& nnz_jac_g,
                              Ipopt::Index& nnz_h_lag, IndexStyleEnum& index_style);
    
    /** Method to return the bounds for my problem */
    virtual bool get_bounds_info(Ipopt::Index n, Ipopt::Number* x_l, Ipopt::Number* x_u,
                                 Ipopt::Index m, Ipopt::Number* g_l, Ipopt::Number* g_u);
    
    /** Method to return the starting point for the algorithm */
    virtual bool get_starting_point(Ipopt::Index n, bool init_x, Ipopt::Number* x,
                                    bool init_z, Ipopt::Number* z_L, Ipopt::Number* z_U,
                                    Ipopt::Index m, bool init_lambda,
                                    Ipopt::Number* lambda);
    
    /** Method to return the objective value */
    virtual bool eval_f(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Number& obj_value);
    
    /** Method to return the gradient of the objective */
    virtual bool eval_grad_f(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Number* grad_f);
    
    /** Method to return the constraint residuals */
    virtual bool eval_g(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Index m, Ipopt::Number* g);
    
    /** Method to return:
     *   1) The structure of the jacobian (if "values" is NULL)
     *   2) The values of the jacobian (if "values" is not NULL)
     */
    virtual bool eval_jac_g(Ipopt::Index n, const Ipopt::Number* x, bool new_x,
                            Ipopt::Index m, Ipopt::Index nele_jac, Ipopt::Index* iRow, Ipopt::Index *jCol,
                            Ipopt::Number* values);
    
    /** Method to return:
     *   1) The structure of the hessian of the lagrangian (if "values" is NULL)
     *   2) The values of the hessian of the lagrangian (if "values" is not NULL)
     */
//    virtual bool eval_h(Ipopt::Index n, const Ipopt::Number* x, bool new_x,
//                        Ipopt::Number obj_factor, Ipopt::Index m, const Ipopt::Number* lambda,
//                        bool new_lambda, Ipopt::Index nele_hess, Ipopt::Index* iRow,
//                        Ipopt::Index* jCol, Ipopt::Number* values);
    
    //@}
    
    /** @name Solution Methods */
    //@{
    /** This method is called when the algorithm is complete so the TNLP can store/write the solution */
    virtual void finalize_solution(Ipopt::SolverReturn status,
                                   Ipopt::Index n, const Ipopt::Number* x, const Ipopt::Number* z_L, const Ipopt::Number* z_U,
                                   Ipopt::Index m, const Ipopt::Number* g, const Ipopt::Number* lambda,
                                   Ipopt::Number obj_value,
                                   const Ipopt::IpoptData* ip_data,
                                   Ipopt::IpoptCalculatedQuantities* ip_cq);
    //@}
    
private:
    TriMesh *m_mesh;
    std::vector<double> m_innerAngle3D;
    std::vector<double> m_edgeLengths;
    std::vector<HalfedgeHandle> m_boundary3D;
    std::vector< std::pair<int, int> > m_coincidePair;
};


#endif
