//
//  Geodesic.cc
//  OpenFlipper
//
//  Created by 饒亢 on 2015/1/13.
//
//
#include "PlushPatternGenerator.hh"

#include <Eigen/Sparse>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>
#include <CGAL/convex_hull_3.h>
#include <queue>

/**
 Select vertices with top (threshold)% curvature. Only select vertices curvature < 0 which are tip points.
 Maximum suppression (of nearby vertices) is applied to reduce the number of selection
 @param <#parameter#>
 @return <#retval#>
 @retval <#meaning#>
 */
bool PlushPatternGenerator::calcSelection(std::vector<VertexHandle> &targetVertices, double threshold) {
    int n = m_mesh->n_vertices();
    Eigen::SparseMatrix<double> M(n, n);
    calc_parameterization_weight_matrix(m_mesh, M, Conformal);
    
    Eigen::SparseLU< Eigen::SparseMatrix<double> > solver;
    solver.compute(M);
    if (solver.info() != Eigen::Success)
    {
        std::cout << "Failed to decompose matrix " << std::endl;
        return false;
    }
    
    // Filling boundary condition
    Eigen::VectorXd K_T(n), K_0(n);
    K_T.setZero();
    K_0.setZero();
    double sum_area = 0;
    double sum_K_0 = 0;
    for (VertexHandle v : m_mesh->vertices()) {
        double sumAngle = getSumInnerAngle(m_mesh, v);
        if (m_mesh->is_boundary(v)) {
            K_0[v.idx()] = M_PI - sumAngle;
        } else {
            K_0[v.idx()] = 2*M_PI - sumAngle;
        }
        sum_K_0 += K_0[v.idx()];
    }
    for (FaceHandle f : m_mesh->faces()) {
        HalfedgeHandle heh = m_mesh->halfedge_handle(f);
        sum_area += m_mesh->calc_sector_area(heh);
    }
    for (VertexHandle v : m_mesh->vertices()) {
        for (HalfedgeHandle cvih : m_mesh->vih_range(v)) {
            if (!m_mesh->is_boundary(cvih)) {
                K_T[v.idx()] += m_mesh->calc_sector_area(cvih);
            }
        }
        K_T[v.idx()] *= sum_K_0 / (3 * sum_area);
    }
    
    Eigen::VectorXd conformal_factors(n);
    conformal_factors = solver.solve(K_T - K_0);
    
    // For visualization
    for (VertexHandle v : m_mesh->vertices()) {
        m_mesh->property(distortionVHandle, v) = conformal_factors[v.idx()];
    }
    
    std::vector< std::pair<VertexHandle, double> > sorted_conformal_factors(conformal_factors.size());
    for (VertexHandle v : m_mesh->vertices()) {
        sorted_conformal_factors[v.idx()].first = v;
        sorted_conformal_factors[v.idx()].second = conformal_factors[v.idx()];
    }
//    auto comparator = [](const std::pair<VertexHandle, double> &left, const std::pair<VertexHandle, double> &right) ->bool {
//        return left.second > right.second;
//    };
//    std::sort(sorted_conformal_factors.begin(), sorted_conformal_factors.end(), comparator);
    
    // Maximum suppression
    for (size_t i = 0; i < sorted_conformal_factors.size() && targetVertices.size() < sqrt(m_mesh->n_vertices()); i++) {
        VertexHandle v = sorted_conformal_factors[i].first;
        double conformal_factor = sorted_conformal_factors[i].second;
        // Test for neigbor in 3-rings
        bool isMaximal = true;
        for (VertexHandle cvv : m_mesh->vv_range(v)) {
            for (VertexHandle cvv2 : m_mesh->vv_range(cvv)) {
                for (VertexHandle cvv3 : m_mesh->vv_range(cvv2)) {
                    if (m_mesh->property(distortionVHandle, cvv3) > conformal_factor) {
                        isMaximal = false;
                        break;
                    }
                }
            }
        }
        
        if (isMaximal) {
            targetVertices.push_back(v);
        }
    }
    return true;
}