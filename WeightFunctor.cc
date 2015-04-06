#include "WeightFunctor.hh"

WeightFunctor::WeightFunctor(TriMesh *mesh,
                             const Polyhedron &polyhedron,
                             boost_vertex_descriptor *currentV,
                             const boost::iterator_property_map<std::vector<boost_vertex_descriptor>::iterator, VertexIdPropertyMap>
                             &predecessor_pmap,
                             const EdgeIdPropertyMap &edgeIndexMap) :
m_mesh(mesh), m_polyhedron(polyhedron), m_predecessor_pmap(predecessor_pmap), m_edgeIndexMap(edgeIndexMap), m_currentV(currentV), m_maxEdgeLength(0) {
    for (EdgeIter e_it = m_mesh->edges_begin(); e_it != m_mesh->edges_end(); e_it++) {
        HalfedgeHandle he = m_mesh->halfedge_handle(*e_it, 0);
        TriMesh::Point p1 = m_mesh->point(m_mesh->from_vertex_handle(he));
        TriMesh::Point p2 = m_mesh->point(m_mesh->to_vertex_handle(he));
        
        m_maxEdgeLength = max(m_maxEdgeLength, (p1-p2).norm());
    }
    
    for (VertexIter v_it = m_mesh->vertices_begin(); v_it != m_mesh->vertices_end(); v_it++) {
        double curvature = m_mesh->property(PlushPatternGenerator::maxCurvatureHandle, *v_it);
    }
}

double WeightFunctor::distanceWeight(TriMesh::Point p1, TriMesh::Point p2) const {
    double weight = (p1-p2).norm()/m_maxEdgeLength;
    assert(0 <= weight && weight <= 1);
    return (p1-p2).norm()/m_maxEdgeLength;
}

double WeightFunctor::textureWeight(HalfedgeHandle he1, HalfedgeHandle he2) const {
    // if the two faces along this edge are different color, set weight of this edge to almost 0
    // we encourage path go through the boundary of different colors
    if (!m_mesh->is_boundary(he1) && !m_mesh->is_boundary(he2)
        &&  m_mesh->color(m_mesh->face_handle(he1)) != m_mesh->color(m_mesh->face_handle(he2))) {
        return 1e-9;
    } else {
        return 1;
    }
}

double WeightFunctor::curvatureWeight(VertexHandle v1, VertexHandle v2) const {
    double curvature1 = m_mesh->property(PlushPatternGenerator::maxCurvatureHandle, v1);
    double curvature2 = m_mesh->property(PlushPatternGenerator::maxCurvatureHandle, v2);

    // weight using curvature
    // encourage path with +1 curvature, not 0
    double weight = 1 / (1 + (abs(curvature1) + abs(curvature2)));
    assert(0 <= weight && weight <= 1);
    
    return weight;
}

/// This measures the direction similarity between edge and corresponding bones
double WeightFunctor::skeletonWeight(VertexHandle v1,
                                     VertexHandle v2,
                                     TriMesh::Point p1,
                                     TriMesh::Point p2) const {
    Skeleton *skeleton = m_mesh->property(PlushPatternGenerator::skeletonHandle);
    
    // Calculate corresponding averaged bone direction for each vertex
    OpenMesh::Vec3d avgBoneDirection(0,0,0);
    double sumOfSqr1 = 0, sumOfSqr2 = 0;
    int i = 0;
    for (std::vector<Bone>::iterator bone_it = skeleton->bones.begin(); bone_it != skeleton->bones.end(); bone_it++, i++) {
        OpenMesh::Vec3d dir = (bone_it->getA() - bone_it->getB()).normalize();
        double weight1 = m_mesh->property(PlushPatternGenerator::bonesWeightHandle, v1)[i];
        double weight2 = m_mesh->property(PlushPatternGenerator::bonesWeightHandle, v2)[i];
        
        assert(!isnan(weight1) && !isnan(weight2));
        avgBoneDirection += dir * (weight1 + weight2)/2;
        
        sumOfSqr1 += weight1 * weight1;
        sumOfSqr2 += weight2 * weight2;
    }
    
    // Calculate standard deviation
    double mean = 1.0/skeleton->bones.size();    // All coefficient sums to 1
    double std1 = sqrt(sumOfSqr1/skeleton->bones.size() - mean * mean);
    double std2 = sqrt(sumOfSqr2/skeleton->bones.size() - mean * mean);
    
    // If there are many influential bones (std is low) for this edge, it's probably near the body part.
    // In such case, skeleton information is not reliable. So we give it a lower influential.
    
    // If one std is low, then this edge is not reliable no matter the other one.
    double std = min(std1, std2);
    double maximumOfStd = sqrt(mean - mean * mean);
    double influential = pow(std/maximumOfStd, 6);
    
    double cosAngle = ((p1-p2)|avgBoneDirection)/((p1-p2).norm()*avgBoneDirection.norm());
    double weight = influential * (1-fabs(cosAngle));
    
    assert(0 <= weight && weight <= 1);
    return weight;
}

//double WeightFunctor::smoothnessWeight(boost_vertex_descriptor boost_v1,
//                                       boost_vertex_descriptor boost_v2,
//                                       TriMesh::Point p1,
//                                       TriMesh::Point p2) const {
//    boost_vertex_descriptor predecessor = m_predecessor_pmap[boost_v1];
//    if (predecessor == boost_v1) {
//        // we reach the begining, assign weight to 1 (max) to ensure that
//        // when combining path, total weight would not increase
//        return 1;
//    } else {
//        auto boost_p0 = predecessor->point();
//        TriMesh::Point p0(boost_p0[0],boost_p0[1],boost_p0[2]);
//        double cosAngle = ((p1-p0) | (p1-p2) / ((p1-p0).norm() * (p1-p2).norm()));
//        // cosAngle = [-1, 1], we need [0, 1]
//        return (cosAngle + 1) / 2;
//    }
//}

double WeightFunctor::operator()(boost_edge_descriptor e) const {
    
    int index = get(m_edgeIndexMap, e);
    EdgeHandle eh = m_mesh->edge_handle(index);
    
    boost_vertex_descriptor boost_v1 = boost::source(e, m_polyhedron);
    boost_vertex_descriptor boost_v2 = boost::target(e, m_polyhedron);
    
    //    assert (boost_v1 == *m_currentV || boost_v2 == *m_currentV);
    
    //    if (boost_v2 == *m_currentV) {
    //        boost_v2 = boost_v1;
    //        boost_v1 = *m_currentV;
    //    }
    
    HalfedgeHandle he1 = m_mesh->halfedge_handle(eh, 0);
    HalfedgeHandle he2 = m_mesh->halfedge_handle(eh, 1);
    
    VertexHandle v1 = m_mesh->from_vertex_handle(he1);
    VertexHandle v2 = m_mesh->to_vertex_handle(he1);
    
    // make sure the direction is correct with boost_vertex_descriptor
    if (v1.idx() != (int)boost_v1->id()) {
        VertexHandle tmp = v1;
        v1 = v2;
        v2 = tmp;
    }
    
    TriMesh::Point p1 = m_mesh->point(v1);
    TriMesh::Point p2 = m_mesh->point(v2);
    
    double coefficients = 0;
    double edgeWeight = 0;
    
    double distanceCoefficient = 1;
    double textureCoefficient = 1;
    double curvatureCoefficient = 5;
    double skeletonCoefficient = 1;
    coefficients += distanceCoefficient + textureCoefficient + curvatureCoefficient + skeletonCoefficient;
    
    if (m_mesh->property(PlushPatternGenerator::edgeWeightHandle, eh) >= 0) {
        edgeWeight = m_mesh->property(PlushPatternGenerator::edgeWeightHandle, eh);
    } else {
        edgeWeight += distanceCoefficient * distanceWeight(p1, p2);
        edgeWeight += textureCoefficient * textureWeight(he1, he2);
        edgeWeight += curvatureCoefficient * curvatureWeight(v1, v2);
        edgeWeight += skeletonCoefficient * skeletonWeight(v1, v2, p1, p2);
        
        m_mesh->property(PlushPatternGenerator::edgeWeightHandle, eh) = edgeWeight;
    }
    
    // re-calculate smoothness weight every time because it depends on path.
    // it can not be saved and reuse
//    double pathCoefficient = 0;
    double pathWeight = 0;
//    pathCoefficient * smoothnessWeight(boost_v1, boost_v2, p1, p2);
//    coefficients += pathCoefficient;
    
    return (edgeWeight + pathWeight) / coefficients;
}