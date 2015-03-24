#include "PlushPatternGenerator.hh"
#include <Eigen/Sparse>

#include <IpIpoptApplication.hpp>
#include <LPFB.hh>

#include <queue>

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
    Ipopt::SmartPtr<Ipopt::TNLP> mynlp = new LPFB_NLP(&mesh);
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

bool PlushPatternGenerator::packFlattenedGraph(const int nColumns) {
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
        height = max(height, maxY[i] - minY[i]);
        depth = max(depth, maxZ[i] - minZ[i]);
    }
    
    double spacing = max(height, width)/5;
    double offsetX = 0, offsetY = 0, offsetZ = 0;
    for (size_t i = 0; i < n; i++) {
        TriMesh &mesh = m_flattenedGraph[i];

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
