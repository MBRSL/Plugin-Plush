#include "PlushPatternGenerator.hh"
#include <unsupported/Eigen/NonLinearOptimization>
#include <Eigen/Sparse>

#include <queue>

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
     Some internal parameters can be adjust. They are:
     lambdaTheta The value of Lagrange multiplier lambda theta.
     lambdaTheta_0x The value of Lagrange multiplier lambda theta 0x.
     lambdaTheta_0y The value of Lagrange multiplier lambda theta 0y.

     @param n Number of vertices(edges) in loop.
     @param edgeLengths Length of edge (v1, v2), (v2, v3),...,(vn-1, vn).
     @param innerAngle3D Inner angles for v1, v2,...,vn.
     */
    lmder_functor(const int n, const std::vector<double> &edgeLengths, const std::vector<double> &innerAngle3D)
        : Functor<double>(n, n), n(n), edgeLengths(edgeLengths), innerAngle3D(innerAngle3D) {
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

bool PlushPatternGenerator::calcFlattenedGraph()
{
    // Clear previous result
    m_flattenedGraph.clear();

    std::vector< std::vector<HalfedgeHandle> > boundaries;
    getBoundariesByEdges(boundaries, m_boundary);

    splitWithBoundary(boundaries);

    isJobCanceled = false;
    for (size_t i = 0; i < m_flattenedGraph.size(); i++) {
        if (isJobCanceled) {
            emit log(LOGINFO, "Flattening calculation canceled.");
            return false;
        }
        calcLPFB(m_flattenedGraph[i]);
        setJobState((double)(i+1) / m_flattenedGraph.size() * 100);
    }

    calcInteriorPoints();

    packFlattenedGraph();
    return true;
}

/// Split existing mesh into (new) connected components with boundary as separator
bool PlushPatternGenerator::splitWithBoundary(std::vector< std::vector<HalfedgeHandle> > &loops)
{
    for (size_t i = 0; i < loops.size(); i++) {
        TriMesh newMesh;
        newMesh.add_property(inverseMapping);
        
        std::vector<HalfedgeHandle> loop = loops[i];
        
        std::set<VertexHandle> isVisited;

        // Starting BFS from boundary halfedges
        std::set<FaceHandle> connectedFaces;
        std::queue<HalfedgeHandle> queue;
        
        // Old vertex to new vertex mapping function
        std::map<VertexHandle, VertexHandle> vertexToVertexMapping;
        
        // Because there will be more than one corresponding (new) vertex for a vertex on loop,
        // vertex-to-vertex is not a one-to-one function. We need a halfedge-to-vertices mapping for boundary halfedge.
        std::map< HalfedgeHandle, std::pair<VertexHandle, VertexHandle> > halfedgeToVerticesMapping;
        
        for (size_t j = 0; j < loop.size(); j++) {
            HalfedgeHandle boundary_heh = loop[j];
            queue.push(boundary_heh);

            VertexHandle v1 = m_mesh->from_vertex_handle(boundary_heh);
            VertexHandle v2 = m_mesh->to_vertex_handle(boundary_heh);
            VertexHandle newV1;
            if (j == 0) {
                newV1 = newMesh.add_vertex(m_mesh->point(v1));
                newMesh.property(inverseMapping, newV1) = v1;
            } else {
                newV1 = halfedgeToVerticesMapping[loop[j-1]].second;
            }
            
            VertexHandle newV2;
            if (j == loop.size()-1) {
                newV2 = halfedgeToVerticesMapping[loop[0]].first;
            } else {
                newV2 = newMesh.add_vertex(m_mesh->point(v2));
                newMesh.property(inverseMapping, newV2) = v2;
            }
            
            isVisited.insert(v1);
            
            halfedgeToVerticesMapping.emplace(boundary_heh, std::make_pair(newV1, newV2));
        }

        while (!queue.empty()) {
            HalfedgeHandle starting_heh = queue.front();
            queue.pop();

            VertexHandle v = m_mesh->to_vertex_handle(starting_heh);
            
            // Add new vertex only if starting_heh is not a boundary halfedge
            bool isBoundary = false;
            if (std::find(loop.begin(), loop.end(), starting_heh) == loop.end()) {
                VertexHandle newV = newMesh.add_vertex(m_mesh->point(v));
                vertexToVertexMapping.emplace(v, newV);
                newMesh.property(inverseMapping, newV) = v;
            } else {
                isBoundary = true;
            }

            // Loop until we meets boundary halfedge or starting halfedge
            HalfedgeHandle vih = m_mesh->opposite_halfedge_handle(m_mesh->next_halfedge_handle(starting_heh));
            while (vih != starting_heh && std::find(loop.begin(), loop.end(), m_mesh->opposite_halfedge_handle(vih)) == loop.end()) {
                connectedFaces.insert(m_mesh->face_handle(vih));
            
                VertexHandle nextV = m_mesh->from_vertex_handle(vih);
            
                if (isVisited.find(nextV) == isVisited.end()) {
                    isVisited.insert(nextV);
                    queue.push(m_mesh->opposite_halfedge_handle(vih));
                }
                
                if (isBoundary) {
                    // We don't care about the source, just use the dest vertex.
                    VertexHandle newV = halfedgeToVerticesMapping[starting_heh].second;
                    halfedgeToVerticesMapping.emplace(vih, std::make_pair(VertexHandle(), newV));
                }
                
                vih = m_mesh->opposite_halfedge_handle(m_mesh->next_halfedge_handle(vih));
            }
        }
        
        newMesh.request_face_colors();
        
        /**
         *          v3
         *         /  \
         *        /    \
         *      v1 ---- v2
         *          heh
         */
        for (std::set<FaceHandle>::iterator f_it = connectedFaces.begin(); f_it != connectedFaces.end(); f_it++) {
            HalfedgeHandle heh12 = m_mesh->halfedge_handle(*f_it);
            HalfedgeHandle heh23 = m_mesh->next_halfedge_handle(heh12);
            HalfedgeHandle heh31 = m_mesh->next_halfedge_handle(heh23);
            VertexHandle v1 = m_mesh->from_vertex_handle(heh12);
            VertexHandle v2 = m_mesh->to_vertex_handle(heh12);
            VertexHandle v3 = m_mesh->to_vertex_handle(heh23);

            std::vector<VertexHandle> face_vhs;
            VertexHandle newV1 = vertexToVertexMapping[v1];
            VertexHandle newV2 = vertexToVertexMapping[v2];
            VertexHandle newV3 = vertexToVertexMapping[v3];
            // Invalid vertex means they are on the boundary, we should use halfedgeToVerticesMapping to find correct vertex
            if (!newMesh.is_valid_handle(newV1)) {
                newV1 = halfedgeToVerticesMapping[heh31].second;
            }
            if (!newMesh.is_valid_handle(newV2)) {
                newV2 = halfedgeToVerticesMapping[heh12].second;
            }
            if (!newMesh.is_valid_handle(newV3)) {
                newV3 = halfedgeToVerticesMapping[heh23].second;
            }
            face_vhs.push_back(newV1);
            face_vhs.push_back(newV2);
            face_vhs.push_back(newV3);
            FaceHandle newF = newMesh.add_face(face_vhs);
            newMesh.set_color(newF, m_mesh->color(*f_it));
        }
        m_flattenedGraph.push_back(newMesh);
    }
    return true;
}

bool PlushPatternGenerator::calcInteriorPoints()
{
    for (size_t i = 0; i < m_flattenedGraph.size(); i++) {
        TriMesh &mesh = m_flattenedGraph[i];
        int n = mesh.n_vertices();
        
        std::vector<HalfedgeHandle> boundary;
        getBoundaryOfOpenedMesh(boundary, mesh);

        // Prepare for matrix
        std::vector< Eigen::Triplet<double> > tripletList;
        double *rowSum = new double[n];
        memset(rowSum, 0, sizeof(double)*n);

        /**     v3
         *     /  \
         *    /    \
         *  v1 ---> v2
         *    \    /
         *     \  /
         *      v4
         */
        for (EdgeIter e_it = mesh.edges_begin(); e_it != mesh.edges_end(); e_it++) {
            HalfedgeHandle heh12 = mesh.halfedge_handle(*e_it, 0);
            HalfedgeHandle heh21 = mesh.halfedge_handle(*e_it, 1);
            
            HalfedgeHandle heh23 = mesh.next_halfedge_handle(heh12);
            HalfedgeHandle heh14 = mesh.next_halfedge_handle(heh21);

            HalfedgeHandle heh31 = mesh.next_halfedge_handle(heh23);
            HalfedgeHandle heh41 = mesh.opposite_halfedge_handle(heh14);

            VertexHandle v1 = mesh.from_vertex_handle(heh12);
            VertexHandle v2 = mesh.to_vertex_handle(heh12);
            
            if (mesh.is_boundary(v1) && mesh.is_boundary(v2)) {
                continue;
            }
            
            double angle231 = mesh.calc_sector_angle(heh23);
            double angle142 = mesh.calc_sector_angle(heh14);
            
            double angle312 = mesh.calc_sector_angle(heh31);
            double angle412 = mesh.calc_sector_angle(heh41);

            if (angle312 < 0) {
                angle312 *= -1;
            }
            if (angle412 < 0) {
                angle412 *= -1;
            }

//            double cotR231 = tan(M_PI_2 - angle231);
//            double cotR142 = tan(M_PI_2 - angle142);
//            double cotSum = cotR231;
//            double cotSum = 1;
            double cotSum = (tan(angle312/2) + tan(angle412/2)) / mesh.calc_edge_length(heh12);
            assert(cotSum >= 0);
            
            if (!mesh.is_boundary(v1)) {
                tripletList.push_back(Eigen::Triplet<double>(v1.idx(), v2.idx(), cotSum));
                rowSum[v1.idx()] += cotSum;
            }
            if (!mesh.is_boundary(v2)) {
                tripletList.push_back(Eigen::Triplet<double>(v2.idx(), v1.idx(), cotSum));
                rowSum[v2.idx()] += cotSum;
            }
        }
        
        for (VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); v_it++) {
            VertexHandle v = *v_it;
            if (mesh.is_boundary(v)) {
                tripletList.push_back(Eigen::Triplet<double>(v.idx(), v.idx(), 1));
            } else {
                tripletList.push_back(Eigen::Triplet<double>(v.idx(), v.idx(), -rowSum[v.idx()]));
            }
        }
        
        // Filling matrix
        Eigen::SparseMatrix<double> M(n, n);
        M.setFromTriplets(tripletList.begin(), tripletList.end());
        
        Eigen::SparseLU< Eigen::SparseMatrix<double> > solver;
        M.makeCompressed();
        solver.compute(M);
        if (solver.info() != Eigen::Success)
        {
            std::cout << "Failed to decompose matrix " << std::endl;
            return false;
        }
        
        // Filling boundary condition
        Eigen::VectorXd Vx(n), Vy(n);
        Vx.setZero();
        Vy.setZero();
        for (size_t j = 0; j < boundary.size(); j++) {
            VertexHandle v = mesh.from_vertex_handle(boundary[j]);
            TriMesh::Point p = mesh.point(v);
            Vx[v.idx()] = p[0];
            Vy[v.idx()] = p[1];
        }

        Eigen::VectorXd Ux(n), Uy(n);
        Ux = solver.solve(Vx);
        Uy = solver.solve(Vy);
        
        // Assign planar coordinates to vertices
        for (VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); v_it++) {
            VertexHandle v = *v_it;
            TriMesh::Point &p = mesh.point(v);
            p[0] = Ux[v.idx()];
            p[1] = Uy[v.idx()];
            p[2] = 0;
        }
    }
    return true;
}

/**
 This function will calculate LPFB of given mesh with loops as boundary.
 The result is stored in planarCoordHandle property in boundary vertices.
 @param <#parameter#>
 @return <#retval#>
 @retval <#meaning#>
 */
bool PlushPatternGenerator::calcLPFB(TriMesh &mesh) {
    
    std::vector<HalfedgeHandle> boundary3D;
    getBoundaryOfOpenedMesh(boundary3D, mesh);
    int n = boundary3D.size();
    // Calculate inner angles for later use.
    Eigen::VectorXd theta(n);
    std::vector<double> innerAngle3D;
    std::vector<double> edgeLengths;
    for (int j = 0; j < n; j++) {
        double sumInnerAngle = PlushPatternGenerator::getSumInnerAngle(&mesh, boundary3D[(j-1+n) % n], boundary3D[j]);
        innerAngle3D.push_back(sumInnerAngle);

        TriMesh::Point v1 = mesh.point(mesh.from_vertex_handle(boundary3D[j]));
        TriMesh::Point v2 = mesh.point(mesh.to_vertex_handle(boundary3D[j]));
        double length = (v1-v2).norm();
        edgeLengths.push_back(length);

        // Initialize parameters
        theta[j] = sumInnerAngle;
    }
    
//    lmder_functor functor(n, edgeLengths, innerAngle3D);
//    Eigen::LevenbergMarquardt<lmder_functor> lm(functor);
        Eigen::NumericalDiff<lmder_functor> functor(n, edgeLengths, innerAngle3D);
        Eigen::LevenbergMarquardt<Eigen::NumericalDiff<lmder_functor> > lm(functor);
    lm.parameters.maxfev = 2000;
    lm.parameters.epsfcn = 1e-3;
    lm.minimize(theta);

    // Reconstruct flattened loop using optimized angles
    double phi_k = 0;
    double sumX = 0, sumY = 0;
    for (int k = 0; k <= n; k++) {
        if (k > 0) {
            phi_k += M_PI - theta[k-1];

            double lk = edgeLengths[k-1];
            
            sumX += lk * cos(phi_k);
            sumY += lk * sin(phi_k);
            
//            printf("%d | v: %d | theta: %.3lf | phi: %.3lf\n", k, m_mesh->from_vertex_handle(loop[k-1]).idx(), theta[k-1], phi_k);
        }
        // If k == 0, it will be placed at (0,0,0)
        VertexHandle v = mesh.from_vertex_handle(boundary3D[k % n]);
        TriMesh::Point &p = mesh.point(v);
        p[0] = sumX;
        p[1] = sumY;
        p[2] = 0;
    }
    return true;
}

bool PlushPatternGenerator::packFlattenedGraph() {
    size_t n = m_flattenedGraph.size();
    
    double *minX = new double[n];
    double *maxX = new double[n];
    double *minY = new double[n];
    double *maxY = new double[n];
    double *minZ = new double[n];
    double *maxZ = new double[n];
    double width = 0, height = 0, depth = 0;
    
    for (size_t i = 0; i < n; i++) {
        minX[i] = 1e9;
        minY[i] = 1e9;
        minZ[i] = 1e9;
        maxX[i] = -1e9;
        maxY[i] = -1e9;
        maxZ[i] = -1e9;
        
        TriMesh &mesh = m_flattenedGraph[i];
        for (VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); v_it++) {
            TriMesh::Point p = mesh.point(*v_it);
            minX[i] = min(minX[i], p[0]);
            maxX[i] = max(maxX[i], p[0]);
            minY[i] = min(minY[i], p[1]);
            maxY[i] = max(maxY[i], p[1]);
            minZ[i] = min(minZ[i], p[2]);
            maxZ[i] = max(maxZ[i], p[2]);
        }
        width = max(width, maxX[i] - minX[i]);
        depth = max(depth, maxZ[i] - minZ[i]);
        height += maxY[i] - minY[i];
    }
    
    double spacing = height / (5 * n);
    height += spacing * n -1;
    double offsetX = 0, offsetY = 0, offsetZ = 0;
    for (size_t i = 0; i < n; i++) {
        TriMesh &mesh = m_flattenedGraph[i];

        // Align to the right bound of previous loop
        if (i > 0) {
            offsetY += maxY[i-1] - minY[i-1] + spacing;
        }

        for (VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); v_it++) {
            VertexHandle v = *v_it;
            TriMesh::Point &p = mesh.point(v);

            p[0] += -minX[i] - width/2;
            p[1] += -minY[i] + offsetY - height/2;
            p[2] += -minZ[i] - depth/2;
        }
    }
    delete[] minX;
    delete[] maxX;
    delete[] minY;
    delete[] maxY;
    delete[] minZ;
    delete[] maxZ;
    return true;
}
