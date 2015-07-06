#include "PlushPatternGenerator.hh"
#include <Eigen/Sparse>

#include <IpIpoptApplication.hpp>
#include <LPFB.hh>

bool PlushPatternGenerator::calcFlattenedGraph()
{
    std::set<EdgeHandle> &seams = m_mesh->property(seams_handle);
;
    if (seams.empty()) {
        emit log(LOGERR, "Seams not ready. Call calcSeams() first.");
        return false;
    }
    
    // Create sub meshes using boundaries
    m_mesh->property(flattenedSubMeshesHandle) = get_subMeshes_with_boundary(seams);

    isJobCanceled = false;
    int progress_counter = 0;

    SubMesh_graph &subMesh_graph = m_mesh->property(flattenedSubMeshesHandle);
    SubMesh_graph::vertex_iterator v_it, v_ite;
    for (boost::tie(v_it, v_ite) = boost::vertices(subMesh_graph); v_it != v_ite; v_it++) {
        if (isJobCanceled) {
            emit log(LOGINFO, "Flattening calculation canceled.");
            return false;
        }
        FilteredTriMesh &subMesh = boost::get(boost::vertex_owner, subMesh_graph, *v_it);

        std::map<HalfedgeHandle, OpenMesh::Vec3d> boundaryPosition;
        calcLPFB(subMesh, boundaryPosition);
        calcInteriorPoints(subMesh, boundaryPosition);
        calcDistortion(subMesh);

        progress_counter++;
        setJobState((double)progress_counter/ boost::num_vertices(subMesh_graph) * 100);
    }
    
    return true;
}

/**
 This function flattens non-boundary vertices using intrinsic parameterization
 @param mesh Input mesh, the flattened vertices position will be stored at mesh.flattened_position.
 @param boundaryPosition The given boundary position used for parameterization.
 @retval True on success. False if failed to solve linear system during intrinsic parameterization.
 */
bool PlushPatternGenerator::calcInteriorPoints(FilteredTriMesh &mesh, std::map<HalfedgeHandle, OpenMesh::Vec3d> &boundaryPosition)
{
    int n_duplicated_vertices = mesh.n_duplicated_vertices();

    Eigen::SparseMatrix<double> M = calc_parameterization_weight_matrix(mesh, MeanValue);

    Eigen::SparseLU< Eigen::SparseMatrix<double> > solver;
    solver.compute(M);
    assert(solver.info() == Eigen::Success && "Failed to decompose matrix ");
    
    // Filling boundary condition
    Eigen::VectorXd Vx(n_duplicated_vertices), Vy(n_duplicated_vertices);
    Vx.setZero();
    Vy.setZero();
    for (auto position_pair : boundaryPosition) {
        HalfedgeHandle heh = position_pair.first;
        TriMesh::Point &p = mesh.flattened_boundary_point(heh);
        p[0] = position_pair.second[0];
        p[1] = position_pair.second[1];
        p[2] = 0;
        
        int idx = mesh.get_boundary_vertex_idx(heh);
        assert(idx >= 0);
        Vx[idx] = p[0];
        Vy[idx] = p[1];
    }

    Eigen::VectorXd Ux(n_duplicated_vertices), Uy(n_duplicated_vertices);
    Ux = solver.solve(Vx);
    Uy = solver.solve(Vy);
    
    // Assign planar coordinates to vertices
    for (VertexHandle vh : mesh.vertices()) {
        if (!mesh.is_boundary(vh)) {
            int idx = mesh.get_non_boundary_vertex_idx(vh);
            assert(idx != -1);
            TriMesh::Point &p = mesh.flattened_non_boundary_point(vh);
            p[0] = Ux[idx];
            p[1] = Uy[idx];
            p[2] = 0;
        }
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
bool PlushPatternGenerator::calcLPFB(FilteredTriMesh &mesh, std::map<HalfedgeHandle, OpenMesh::Vec3d> &boundaryPosition) {
    Ipopt::SmartPtr<Ipopt::TNLP> mynlp = new LPFB_NLP(mesh, boundaryPosition);
    Ipopt::SmartPtr<Ipopt::IpoptApplication> app = IpoptApplicationFactory();
    
    app->Options()->SetNumericValue("tol", 1e-5);
    app->Options()->SetStringValue("mu_strategy", "adaptive");
    app->Options()->SetStringValue("linear_solver", "ma57");
    // This let us not implementing eval_h()
    app->Options()->SetStringValue("hessian_approximation","limited-memory");
//    app->Options()->SetStringValue("derivative_test","first-order");
    app->Options()->SetIntegerValue("print_level", 0);
    
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

void PlushPatternGenerator::calcDistortion(std::vector<FilteredTriMesh> &flattenedSubMeshes) {
    for (FilteredTriMesh &flattenedMesh : flattenedSubMeshes) {
        calcDistortion(flattenedMesh);
    }
}

/**
 This function calculate the distortion between flattened meshes and original mesh.
 The result is stored at distortionHandle of m_mesh.
 @param <#parameter#>
 */
void PlushPatternGenerator::calcDistortion(FilteredTriMesh &flattenedMesh) {
    // For each face, calculate area difference
    // For each vertex, calculate averaged neighboring area difference
    
    for (FaceHandle fh : flattenedMesh.faces()) {
        double original_area = m_mesh->calc_sector_area(flattenedMesh.halfedge_handle(fh));

        int count = 0;
        TriMesh::Point flattenedPoints[3];
        for (const HalfedgeHandle cfh : flattenedMesh.fh_range(fh)) {
            VertexHandle vh = flattenedMesh.from_vertex_handle(cfh);
            if (!flattenedMesh.is_boundary(vh)) {
                flattenedPoints[count++] = flattenedMesh.flattened_non_boundary_point(vh);
            } else {
                flattenedPoints[count++] = flattenedMesh.flattened_boundary_point(cfh);
            }
        }
        OpenMesh::Vec3d a = flattenedPoints[1] - flattenedPoints[0];
        OpenMesh::Vec3d b = flattenedPoints[2] - flattenedPoints[0];
        double flattened_area = (a % b).norm()/2;

        double area_diff = (flattened_area - original_area) / original_area;
        m_mesh->property(distortionFHandle, fh) = area_diff;
        
        flattenedMesh.max_distortion = max(flattenedMesh.max_distortion, abs(area_diff));
//        m_mesh->property(distortionVHandle, original_v1) += area_diff;
//        m_mesh->property(distortionVHandle, original_v2) += area_diff;
//        m_mesh->property(distortionVHandle, original_v3) += area_diff;
    }
    
    // Do not divide distortionVHandle due to boundary incompleteness in sub mesh
}
