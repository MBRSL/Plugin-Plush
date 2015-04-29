#include "PlushPatternGenerator.hh"
#include <Eigen/Sparse>

#include <IpIpoptApplication.hpp>
#include <LPFB.hh>

/*
 * v4 -------- v1
 *  \         / \
 *   \       /   \
 *    \     /     \
 *     \   /       \
 *      \ /         \
 *      v2 -------- v3
 */
void PlushPatternGenerator::calc_parameterization_weight_matrix(TriMesh *mesh,
                                                                Eigen::SparseMatrix<double> &M,
                                                                InteriorParameterizationMethod method) {
    int n = mesh->n_vertices();
    
    // Prepare for matrix
    std::vector< Eigen::Triplet<double> > tripletList;
    double *rowSum = new double[n];
    memset(rowSum, 0, sizeof(double)*n);
    
    for (EdgeHandle eh : mesh->edges()) {
        HalfedgeHandle heh12 = mesh->halfedge_handle(eh, 0);
        HalfedgeHandle heh21 = mesh->halfedge_handle(eh, 1);
        
        VertexHandle v1 = mesh->from_vertex_handle(heh12);
        VertexHandle v2 = mesh->to_vertex_handle(heh12);
        if (mesh->is_boundary(v1) && mesh->is_boundary(v2)) {
            continue;
        }
        
        HalfedgeHandle heh31 = mesh->prev_halfedge_handle(heh12);
        HalfedgeHandle heh42 = mesh->prev_halfedge_handle(heh21);
        
        double weight12 = 0, weight21 = 0;
        if (method == Barycentric) {
            weight12 = 1;
            weight21 = 1;
        } else if (method == MeanValue) {
            double angle312 = mesh->calc_sector_angle(heh31);
            double angle214 = mesh->calc_sector_angle(heh21);
            double angle123 = mesh->calc_sector_angle(heh12);
            double angle421 = mesh->calc_sector_angle(heh42);
            weight12 = (tan(angle312/2) + tan(angle214/2)) / mesh->calc_edge_length(heh12);
            weight21 = (tan(angle123/2) + tan(angle421/2)) / mesh->calc_edge_length(heh21);
        } else if (method == Conformal) {
            double angle231 = mesh->calc_sector_angle(mesh->next_halfedge_handle(heh12));
            double angle142 = mesh->calc_sector_angle(mesh->next_halfedge_handle(heh21));
            double cotR231 = tan(M_PI_2 - angle231);
            double cotR142 = tan(M_PI_2 - angle142);
            weight12 = cotR231 + cotR142;
            weight21 = weight12;
        } else {
            assert("Invalid InteriorParameterizationMethod");
            return;
        }
        
        if (!mesh->is_boundary(v1)) {
            tripletList.push_back(Eigen::Triplet<double>(v1.idx(), v2.idx(), weight12));
            rowSum[v1.idx()] += weight12;
        }
        if (!mesh->is_boundary(v2)) {
            tripletList.push_back(Eigen::Triplet<double>(v2.idx(), v1.idx(), weight21));
            rowSum[v2.idx()] += weight21;
        }
    }
    
    for (VertexHandle v : mesh->vertices()) {
        if (mesh->is_boundary(v)) {
            tripletList.push_back(Eigen::Triplet<double>(v.idx(), v.idx(), 1));
        } else {
            tripletList.push_back(Eigen::Triplet<double>(v.idx(), v.idx(), -rowSum[v.idx()]));
        }
    }
    // Filling matrix
    M.setFromTriplets(tripletList.begin(), tripletList.end());
    M.makeCompressed();
}

bool PlushPatternGenerator::calcFlattenedGraph()
{
    std::set<EdgeHandle> *seams = getSeams();
    if (!seams) {
        emit log(LOGERR, "Seams not ready. Call calcSeams() first.");
        return false;
    }
    
    // Clear previous result
    std::vector<TriMesh> *flattenedMeshes = &m_mesh->property(flattenedMeshesHandle);
    flattenedMeshes->clear();

    // Create sub meshes using boundaries
    splitWithBoundary(flattenedMeshes, seams);

    isJobCanceled = false;
    for (size_t i = 0; i < flattenedMeshes->size(); i++) {
        if (isJobCanceled) {
            emit log(LOGINFO, "Flattening calculation canceled.");
            return false;
        }
        std::map<VertexHandle, OpenMesh::Vec3d> boundaryPosition;
        calcLPFB(&flattenedMeshes->at(i), &boundaryPosition);
        calcInteriorPoints(&flattenedMeshes->at(i), &boundaryPosition);
        setJobState((double)(i+1) / flattenedMeshes->size() * 100);
    }

    
    calcDistortion(flattenedMeshes);

    packFlattenedGraph(flattenedMeshes);
    return true;
}

/**
 This function flattens non-boundary vertices using intrinsic parameterization
 @param flattenedMeshes Input mesh, its non-boundary vertices will be modified.
 @retval True on success. False if failed to solve linear system during intrinsic parameterization.
 */
bool PlushPatternGenerator::calcInteriorPoints(TriMesh *mesh, std::map<VertexHandle, OpenMesh::Vec3d> *boundaryPosition)
{
    int n = mesh->n_vertices();
    Eigen::SparseMatrix<double> M(n, n);
    calc_parameterization_weight_matrix(mesh, M, MeanValue);

    Eigen::SparseLU< Eigen::SparseMatrix<double> > solver;
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
    for (VertexIter v_it = mesh->vertices_begin(); v_it != mesh->vertices_end(); v_it++) {
        if (mesh->is_boundary(*v_it)) {
            OpenMesh::Vec3d p = boundaryPosition->at(*v_it);
            Vx[v_it->idx()] = p[0];
            Vy[v_it->idx()] = p[1];
        }
    }

    Eigen::VectorXd Ux(n), Uy(n);
    Ux = solver.solve(Vx);
    Uy = solver.solve(Vy);
    
    // Assign planar coordinates to vertices
    for (VertexIter v_it = mesh->vertices_begin(); v_it != mesh->vertices_end(); v_it++) {
        VertexHandle v = *v_it;
        TriMesh::Point &p = mesh->point(v);
        p[0] = Ux[v.idx()];
        p[1] = Uy[v.idx()];
        p[2] = 0;
    }
    return true;
}

/**
 This function will calculate LPFB of given mesh with loops as boundary.
 @param mesh The input mesh.
 @param boundaryPosition The result of boundary position will be stored here.
 @return <#retval#>
 @retval <#meaning#>
 */
bool PlushPatternGenerator::calcLPFB(TriMesh *mesh, std::map<VertexHandle, OpenMesh::Vec3d> *boundaryPosition) {
    Ipopt::SmartPtr<Ipopt::TNLP> mynlp = new LPFB_NLP(mesh, boundaryPosition);
    Ipopt::SmartPtr<Ipopt::IpoptApplication> app = IpoptApplicationFactory();
    
    app->Options()->SetNumericValue("tol", 1e-5);
    app->Options()->SetNumericValue("max_iter", 200);
    app->Options()->SetStringValue("mu_strategy", "adaptive");
    app->Options()->SetStringValue("linear_solver", "ma57");
    // This let us not implementing eval_h()
    app->Options()->SetStringValue("hessian_approximation","limited-memory");
//    app->Options()->SetStringValue("derivative_test","first-order");
    
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

void PlushPatternGenerator::calcDistortion(std::vector<TriMesh> *flattenedMeshes) {
    for (TriMesh &flattenedMesh : *flattenedMeshes) {
        calcDistortion(flattenedMesh);
    }
}

/**
 This function calculate the distortion between flattened meshes and original mesh.
 The result is stored at distortionHandle of m_mesh.
 @param <#parameter#>
 */
void PlushPatternGenerator::calcDistortion(TriMesh &flattenedMesh) {
    // For each face, calculate area difference
    // For each vertex, calculate averaged neighboring area difference
    OpenMesh::VPropHandleT<VertexHandle> inverseMapping = getInverseMappingHandle(&flattenedMesh);
    for (FaceHandle f : flattenedMesh.faces()) {
        HalfedgeHandle heh = flattenedMesh.halfedge_handle(f);
        double area = flattenedMesh.calc_sector_area(heh);

        TriMesh::ConstFaceVertexIter cfv_it = flattenedMesh.cfv_iter(f);
        VertexHandle original_v1 = get_original_handle(&flattenedMesh, *cfv_it++);
        VertexHandle original_v2 = get_original_handle(&flattenedMesh, *cfv_it++);
        VertexHandle original_v3 = get_original_handle(&flattenedMesh, *cfv_it++);

        FaceHandle original_f = get_original_handle(&flattenedMesh, f);
        HalfedgeHandle original_heh = m_mesh->halfedge_handle(original_f);
        double original_area = m_mesh->calc_sector_area(original_heh);

        double area_diff = (area-original_area) / original_area;
        m_mesh->property(distortionFHandle, original_f) = area_diff;
        m_mesh->property(distortionVHandle, original_v1) += area_diff;
        m_mesh->property(distortionVHandle, original_v2) += area_diff;
        m_mesh->property(distortionVHandle, original_v3) += area_diff;
    }
    
    // Do not divide distortionVHandle due to boundary incompleteness in sub mesh
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
    double *columnWidth = new double[nColumns];
    memset(columnWidth, 0, 4*sizeof(double));
    
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
        
        columnWidth[i % nColumns] = max(columnWidth[i % nColumns], maxX[i] - minX[i]);
    }
    
    double spacing = max(height, width)/5;
    double offsetX = 0, offsetY = 0, offsetZ = 0;
    double lastRowHeight = 0, rowHeight = 0;
    for (size_t i = 0; i < n; i++) {
        TriMesh &mesh = flattenedMeshes->at(i);

        // For a new row, calculate row offset
        if (i % nColumns == 0) {
            if (i > 0) {
                lastRowHeight = rowHeight;
                offsetY += lastRowHeight + spacing;

                rowHeight = 0;
                for (int j = 0; j < nColumns; j++) {
                    rowHeight = max(rowHeight, maxY[i+j] - minY[i+j]);
                }
            }
            offsetX = 0;
        } else {
            offsetX += columnWidth[(i-1) % nColumns] + spacing;
        }
        

        for (VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); v_it++) {
            VertexHandle v = *v_it;
            TriMesh::Point &p = mesh.point(v);

            p[0] += -minX[i] + offsetX + columnWidth[i % nColumns]/2 -(maxX[i] - minX[i])/2;
            p[1] += -minY[i] + offsetY + rowHeight/2 - (maxY[i] - minY[i])/2;
            p[2] += -minZ[i] - depth/2;
        }
    }
    delete[] minX;
    delete[] maxX;
    delete[] minY;
    delete[] maxY;
    delete[] minZ;
    delete[] maxZ;
    delete[] columnWidth;
    return true;
}
