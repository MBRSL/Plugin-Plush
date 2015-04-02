#include "PlushPatternGenerator.hh"
#include <Eigen/Sparse>

#include <IpIpoptApplication.hpp>
#include <LPFB.hh>

bool PlushPatternGenerator::calcFlattenedGraph()
{
    std::vector<EdgeHandle> *seams = getSeams();
    if (!seams) {
        emit log(LOGERR, "Seams not ready. Call calcSeams() first.");
        return false;
    }
    
    // Clear previous result
    std::vector<TriMesh> *flattenedMeshes = &m_mesh->property(flattenedMeshesHandle);
    flattenedMeshes->clear();

    // Create sub meshes using boundaries
    splitWithBoundary(flattenedMeshes);

    isJobCanceled = false;
    for (size_t i = 0; i < flattenedMeshes->size(); i++) {
        if (isJobCanceled) {
            emit log(LOGINFO, "Flattening calculation canceled.");
            return false;
        }
        calcLPFB(&flattenedMeshes->at(i));
        setJobState((double)(i+1) / flattenedMeshes->size() * 100);
    }

    calcInteriorPoints(flattenedMeshes);

    packFlattenedGraph(flattenedMeshes);
    return true;
}

/**
 This function flattens non-boundary vertices using intrinsic parameterization (mean value coordinates.)
 @param flattenedMeshes Input mesh, its non-boundary vertices will be modified.
 @retval True on success. False if failed to solve linear system during intrinsic parameterization.
 */
bool PlushPatternGenerator::calcInteriorPoints(std::vector<TriMesh> *flattenedMeshes)
{
    for (size_t i = 0; i < flattenedMeshes->size(); i++) {
        TriMesh &mesh = flattenedMeshes->at(i);
        int n = mesh.n_vertices();
        
        // Prepare for matrix
        std::vector< Eigen::Triplet<double> > tripletList;
        double *rowSum = new double[n];
        memset(rowSum, 0, sizeof(double)*n);

        /**     v3
         *     /  \
         *    /    \
         *  v1⥫---⥬v2
         *    \    /
         *     \  /
         *      v4
         */
        for (EdgeIter e_it = mesh.edges_begin(); e_it != mesh.edges_end(); e_it++) {
            HalfedgeHandle heh12 = mesh.halfedge_handle(*e_it, 0);
            HalfedgeHandle heh21 = mesh.halfedge_handle(*e_it, 1);
            
            VertexHandle v1 = mesh.from_vertex_handle(heh12);
            VertexHandle v2 = mesh.to_vertex_handle(heh12);
            if (mesh.is_boundary(v1) && mesh.is_boundary(v2)) {
                continue;
            }

            HalfedgeHandle heh31 = mesh.prev_halfedge_handle(heh12);
            HalfedgeHandle heh42 = mesh.prev_halfedge_handle(heh21);

            double angle312 = mesh.calc_sector_angle(heh31);
            double angle214 = mesh.calc_sector_angle(heh21);
        
            double angle123 = mesh.calc_sector_angle(heh12);
            double angle421 = mesh.calc_sector_angle(heh42);


//            // Conformal mapping
//            double cotR231 = tan(M_PI_2 - angle231);
//            double cotR142 = tan(M_PI_2 - angle142);
//            double weight12 = cotR231 + cotR142;
//            double weight21 = weight12;
            
//            // Barycentrical mapping
//            double weight12 = 1;
//            double weight21 = 1;
            
            // Mean value mapping
            double wegith12 = (tan(angle312/2) + tan(angle214/2)) / mesh.calc_edge_length(heh12);
            double wegith21 = (tan(angle123/2) + tan(angle421/2)) / mesh.calc_edge_length(heh21);
        
            if (!mesh.is_boundary(v1)) {
                tripletList.push_back(Eigen::Triplet<double>(v1.idx(), v2.idx(), wegith12));
                rowSum[v1.idx()] += wegith12;
            }
            if (!mesh.is_boundary(v2)) {
                tripletList.push_back(Eigen::Triplet<double>(v2.idx(), v1.idx(), wegith21));
                rowSum[v2.idx()] += wegith21;
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
        M.makeCompressed();
        
        Eigen::SparseLU< Eigen::SparseMatrix<double> > solver;
        solver.compute(M);
        if (solver.info() != Eigen::Success)
        {
            std::cout << "Failed to decompose matrix " << std::endl;
            return false;
        }
        
        // Filling boundary condition
        std::vector< std::vector<HalfedgeHandle> > boundaries;
        getBoundaryOfOpenedMesh(boundaries, &mesh, false);

        Eigen::VectorXd Vx(n), Vy(n);
        Vx.setZero();
        Vy.setZero();
        for (size_t j = 0; j < boundaries.size(); j++) {
            std::vector<HalfedgeHandle> &boundary = boundaries[j];
            for (size_t k = 0; k < boundary.size(); k++) {
                VertexHandle v = mesh.from_vertex_handle(boundary[k]);
                TriMesh::Point p = mesh.point(v);
                Vx[v.idx()] = p[0];
                Vy[v.idx()] = p[1];
            }
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
bool PlushPatternGenerator::calcLPFB(TriMesh *mesh) {
    Ipopt::SmartPtr<Ipopt::TNLP> mynlp = new LPFB_NLP(mesh);
    Ipopt::SmartPtr<Ipopt::IpoptApplication> app = IpoptApplicationFactory();
    
    app->Options()->SetNumericValue("tol", 1e-10);
    app->Options()->SetStringValue("mu_strategy", "adaptive");
    app->Options()->SetStringValue("linear_solver", "ma57");
    // This let us not implementing eval_h()
    app->Options()->SetStringValue("hessian_approximation","limited-memory");
    
    // Intialize the IpoptApplication and process the options
    Ipopt::ApplicationReturnStatus status;
    status = app->Initialize();
    if (status != Ipopt::Solve_Succeeded) {
        printf("\n\n*** Error during initialization!\n");
        return false;
    }
    
    // Ask Ipopt to solve the problem
    status = app->OptimizeTNLP(mynlp);
    
    if (status == Ipopt::Solve_Succeeded) {
        printf("\n\n*** The problem solved!\n");
        return true;
    }
    else {
        printf("\n\n*** Some problem occured! See messages above\n");
        return false;
    }
    
    // As the SmartPtrs go out of scope, the reference count
    // will be decremented and the objects will automatically
    // be deleted.
}

bool PlushPatternGenerator::packFlattenedGraph(std::vector<TriMesh> *flattenedMeshes, const int nColumns) {
    size_t n = flattenedMeshes->size();
    
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
        
        TriMesh &mesh = flattenedMeshes->at(i);
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
        height = max(height, maxY[i] - minY[i]);
        depth = max(depth, maxZ[i] - minZ[i]);
    }
    
    double spacing = max(height, width)/5;
    double offsetX = 0, offsetY = 0, offsetZ = 0;
    for (size_t i = 0; i < n; i++) {
        TriMesh &mesh = flattenedMeshes->at(i);

        // For a new row, calculate row offset
        if (i > 1 && i % nColumns == 0) {
            double rowOffset = maxY[i-1];
            for (int j = 2; j <= nColumns; j++) {
                rowOffset = max(rowOffset, maxY[i-j] - minY[i-j]);
            }
            offsetY += rowOffset + spacing;
        }

        for (VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); v_it++) {
            VertexHandle v = *v_it;
            TriMesh::Point &p = mesh.point(v);

            p[0] += -minX[i] + width * (-0.5 + (i % nColumns) * spacing);
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
