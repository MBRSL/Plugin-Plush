#include "PlushPatternGenerator.hh"
#include <unsupported/Eigen/NonLinearOptimization>

// Generic functor
template<typename _Scalar, int NX=Eigen::Dynamic, int NY=Eigen::Dynamic>
struct Functor
{
    typedef _Scalar Scalar;
    enum {
        InputsAtCompileTime = NX,
        ValuesAtCompileTime = NY
    };
    typedef Eigen::Matrix<Scalar,InputsAtCompileTime,1> InputType;
    typedef Eigen::Matrix<Scalar,ValuesAtCompileTime,1> ValueType;
    typedef Eigen::Matrix<Scalar,ValuesAtCompileTime,InputsAtCompileTime> JacobianType;
    
    const int m_inputs, m_values;
    
    Functor() : m_inputs(InputsAtCompileTime), m_values(ValuesAtCompileTime) {}
    Functor(int inputs, int values) : m_inputs(inputs), m_values(values) {}
    
    int inputs() const { return m_inputs; }
    int values() const { return m_values; }
    
    // you should define that in the subclass :
    //  void operator() (const InputType& x, ValueType* v, JacobianType* _j=0) const;
};

class lmder_functor : public Functor<double>
{
private:
    /// Number of vertices. It's the same number of edges
    int n;

    std::vector<double> edgeLengths;
    
    /**
     *              inner edges
     *                \  |  /
     *     loop[i-1]   \ | /    loop[i]
     *  -------------->  v -------------->
     *
     * For loop[i], it records the sum of inner angles between loop[i-1] and loop[i]
     */
    std::vector<double> innerAngle3D;
    
    double lambdaTheta = 1;
    double lambdaTheta_0x = 1;
    double lambdaTheta_0y = 1;
    
public:
    /**
     Constructor of cost function object for a given loop.
     @param n Number of vertices(edges) in loop.
     @param edgeLengths Length of edge (v1, v2), (v2, v3),...,(vn-1, vn).
     @param innerAngle3D Inner angles for v1, v2,...,vn.
     @param lambdaTheta The value of Lagrange multiplier lambda theta.
     @param lambdaTheta_0x The value of Lagrange multiplier lambda theta 0x.
     @param lambdaTheta_0y The value of Lagrange multiplier lambda theta 0y.
     */
    lmder_functor(const int n, const std::vector<double> &edgeLengths, const std::vector<double> &innerAngle3D)
        : Functor<double>(n, n),
    n(n), edgeLengths(edgeLengths), innerAngle3D(innerAngle3D) {
//        lambdaTheta(lambdaTheta), lambdaTheta_0x(lambdaTheta_0x), lambdaTheta_0y(lambdaTheta_0y) {
            lambdaTheta = 1;
            lambdaTheta_0x = 1;
            lambdaTheta_0y = 1;
        }
    
    /// Cost function
    int operator()(const Eigen::VectorXd &x, Eigen::VectorXd &fvec) const
    {
        // Our cost function consist of only 4 terms. Set others to 0.
        fvec.setZero();
        
        double sumAngleDiff = 0;

        double loopAngleConstrain = (n - 2) * M_PI;
        
        double loopPositionConstrainsX = 0;
        double loopPositionConstrainsY = 0;
        double phi_i = 0;
        for (int i = 0; i < n; i++) {
            double theta_i = x[i];
            
            double angleDiff = theta_i - innerAngle3D[i];
            sumAngleDiff += 0.5 * angleDiff * angleDiff;
            
            loopAngleConstrain -= theta_i;
            
            phi_i += M_PI - theta_i;
            loopPositionConstrainsX += edgeLengths[i] * cos(phi_i);
            loopPositionConstrainsY += edgeLengths[i] * sin(phi_i);
        }
        loopAngleConstrain *= lambdaTheta;
        loopPositionConstrainsX *= lambdaTheta_0x;
        loopPositionConstrainsY *= lambdaTheta_0y;
        
        fvec[0] = sumAngleDiff;
        fvec[1] = loopAngleConstrain;
        fvec[2] = loopPositionConstrainsX;
        fvec[3] = loopPositionConstrainsY;
        
//        printf("%.3lf, %.3lf, %.3lf, %.3lf\n", fvec[0], fvec[1], fvec[2], fvec[3]);
        return 0;
    }

    int df(const Eigen::VectorXd &x, Eigen::MatrixXd &fjac) const
    {
        fjac.setZero();
        double *phi_k = new double[n];
        double *sumPosX = new double[n];
        double *sumPosY = new double[n];
        
        phi_k[0] = M_PI - x[0];
        for (int k = 1; k < n; k++) {
            phi_k[k] = phi_k[k-1] + M_PI - x[k];
        }
        double sumPosXTmp = 0, sumPosYTmp = 0;
        for (int k = n-1; k >= 0; k--) {
            sumPosXTmp += edgeLengths[k] * sin(phi_k[k]);
            sumPosYTmp += edgeLengths[k] * -cos(phi_k[k]);
            sumPosX[k] = sumPosXTmp;
            sumPosY[k] = sumPosYTmp;
        }
        
        for (int i = 0; i < n; i++) {
            fjac(0, i) = x[i] - innerAngle3D[i];
            fjac(1, i) = -lambdaTheta;
            fjac(2, i) = lambdaTheta_0x * sumPosX[i];
            fjac(3, i) = lambdaTheta_0y * sumPosY[i];
        }
        delete[] phi_k;
        return 0;
    }
};

bool PlushPatternGenerator::calcFlattenedGraph(std::vector< std::vector<HalfedgeHandle> > &loops)
{
    isJobCanceled = false;
    m_flattenedGraph.clear();
    for (size_t i = 0; i < loops.size(); i++) {
        if (isJobCanceled) {
            emit log(LOGINFO, "Flattening calculation canceled.");
            return false;
        }

        std::vector<HalfedgeHandle> loop = loops[i];
        int n = loop.size();
        
        // Calculate inner angles for later use.
        Eigen::VectorXd theta(n);
        std::vector<double> innerAngle3D;
        std::vector<double> edgeLengths;
        for (int j = 0; j < n; j++) {
            double sumInnerAngle = PlushPatternGenerator::getSumInnerAngle(m_mesh, loop[(j-1+n) % n], loop[j]);
            innerAngle3D.push_back(sumInnerAngle);

            TriMesh::Point v1 = m_mesh->point(m_mesh->from_vertex_handle(loop[j]));
            TriMesh::Point v2 = m_mesh->point(m_mesh->to_vertex_handle(loop[j]));
            double length = (v1-v2).norm();
            edgeLengths.push_back(length);

            // Initialize parameters
            theta[j] = sumInnerAngle;
        }
        
        lmder_functor functor(n, edgeLengths, innerAngle3D);
        Eigen::LevenbergMarquardt<lmder_functor> lm(functor);
//        Eigen::NumericalDiff<lmder_functor> functor(n, edgeLengths, innerAngle3D);
//        Eigen::LevenbergMarquardt<Eigen::NumericalDiff<lmder_functor> > lm(functor);
        lm.parameters.maxfev = 400;
        lm.parameters.epsfcn = 1e-2;
        lm.minimize(theta);
        
        // Reconstruct flattened loop using optimized angles
        PolyLine loop2D;
        double phi_k = 0;
        double sumX = 0, sumY = 0;
        for (int k = 0; k <= n; k++) {
            if (k > 0) {
                phi_k += M_PI - theta[k-1];

                double lk = edgeLengths[k-1];
                
                sumX += lk * cos(phi_k);
                sumY += lk * sin(phi_k);
                
//                printf("%d | v: %d | theta: %.3lf | phi: %.3lf\n", k, m_mesh->from_vertex_handle(loop[k-1]).idx(), theta[k-1], phi_k);
            }
            // If k == 0, it will be placed at (0,0,0)
            TriMesh::Point p(sumX, sumY, 0);
            loop2D.add_point(p);
        }
        m_flattenedGraph.push_back(loop2D);
        setJobState((double)(i+1) / loops.size() * 100);
    }
    return true;
}

bool PlushPatternGenerator::packFlattenedGraph() {
    int n = m_flattenedGraph.size();
    double *minX = new double[n];
    double *maxX = new double[n];
    double *minY = new double[n];
    double *maxY = new double[n];
    double width = 0, height = 0;
    for (int i = 0; i < n; i++) {
        PolyLine &loop = m_flattenedGraph[i];
        for (size_t j = 0; j < loop.n_vertices(); j++) {
            TriMesh::Point p = loop.point(j);
            minX[i] = min(minX[i], p[0]);
            maxX[i] = max(maxX[i], p[0]);
            minY[i] = min(minY[i], p[1]);
            maxY[i] = max(maxY[i], p[1]);
        }
        width = max(width, maxX[i] - minX[i]);
        height += maxY[i] - minY[i];
    }
    
    double spacing = height / (5 * n);
    height += spacing * n -1;
    double offsetX = 0, offsetY = 0;
    for (int i = 0; i < n; i++) {
        PolyLine &loop = m_flattenedGraph[i];
        
        // Align to the right bound of previous loop
        if (i > 0) {
            offsetY += maxY[i-1] - minY[i-1] + spacing;
        }
        for (size_t j = 0; j < loop.n_vertices(); j++) {
            PolyLine::Point &p = loop.point(j);
            p[0] += -minX[i] - width/2;
            p[1] += -minY[i] + offsetY - height/2;
        }
    }
    
    delete[] minX;
    delete[] maxX;
    delete[] minY;
    delete[] maxY;
    return true;
}

/**
 * @brief This function divide a set of edges (with loops) into several loops.
 * @param loops Result loops will be stored here.
 * @param selectedVertices (This parameter should be compute internally) The set of intersection point of m_spanningTree.
 */
void PlushPatternGenerator::getLoops(std::vector< std::vector<HalfedgeHandle> > &loops, std::vector<VertexHandle> &selectedVertices) {
    std::set<HalfedgeHandle> isVisited;
    for (size_t i = 0; i < selectedVertices.size(); i++) {
        VertexHandle startV = selectedVertices[i];
        
        // Starting from a intersection point, search for its adjacenct edges
        for (TriMesh::VertexOHalfedgeIter voh_it = m_mesh->voh_iter(startV); voh_it; voh_it++) {
            // Skip visited edges and non-spanning-tree edges
            if (std::find(isVisited.begin(), isVisited.end(), *voh_it) != isVisited.end()
                ||  std::find(m_spanningTree.begin(), m_spanningTree.end(), m_mesh->edge_handle(*voh_it)) == m_spanningTree.end()) {
                continue;
            }
            
            // We got a non-visited edge, lets start searching from this edge
            HalfedgeHandle start_heh = *voh_it;
            HalfedgeHandle current_heh = start_heh;
            isVisited.insert(start_heh);
            
            std::vector<HalfedgeHandle> loop;
            loop.push_back(current_heh);
            do {
                assert(!m_mesh->is_boundary(current_heh));
                
                HalfedgeHandle next_heh = current_heh;
                // Find next edge which is in spanning tree
                do {
                    // While there are possibly many adjacent edges to choose, we use next->opposite to iterator through them
                    // so that we will visit them all clockwise without missing any one.
                    next_heh = m_mesh->opposite_halfedge_handle(
                                                                m_mesh->next_halfedge_handle(next_heh));
                    
                    // This means there are no other edges except current_heh (which is where you came from)
                    // In this case, just go back and treat this path as a dart
                    if (next_heh == current_heh) {
                        // Do nothing, it will be flipped later.
                    }
                } while (std::find(m_spanningTree.begin(), m_spanningTree.end(), m_mesh->edge_handle(next_heh)) == m_spanningTree.end());
                
                // Flip it to point to the next vertex
                current_heh = m_mesh->opposite_halfedge_handle(next_heh);
                isVisited.insert(current_heh);
                
                loop.push_back(current_heh);
            } while (current_heh != start_heh);
            
            // Delete last duplicated halfedge
            loop.pop_back();
            
            loops.push_back(loop);
        }
    }
}