#include "WeightFunctor.hh"
#include "PlushPatternGenerator.hh"

WeightFunctor::WeightFunctor(TriMesh *mesh,
                             VertexHandle &currentV,
                             const boost::iterator_property_map<std::vector<VertexHandle>::iterator, TriMesh_id_map>
                                 *predecessor_pmap,
                             double distanceCoefficient,
                             double textureCoefficient,
                             double curvatureCoefficient,
                             double skeletonCoefficient,
                             double pathCoefficient) :
m_mesh(mesh),
m_currentV(currentV),
m_predecessor_pmap(predecessor_pmap),
m_maxEdgeLength(0),
m_distanceCoefficient(distanceCoefficient),
m_textureCoefficient(textureCoefficient),
m_curvatureCoefficient(curvatureCoefficient),
m_skeletonCoefficient(skeletonCoefficient),
m_pathCoefficient(pathCoefficient)
{
    for (EdgeHandle eh : m_mesh->edges()) {
        m_maxEdgeLength = max(m_maxEdgeLength, m_mesh->calc_edge_length(eh));
    }
}

double WeightFunctor::distanceWeight(EdgeHandle eh) const {
    double weight = m_mesh->calc_edge_length(eh)/m_maxEdgeLength;
    assert(0 <= weight && weight <= 1);
    m_mesh->property(PlushPatternGenerator::geodesic_distance_weight_Handle, eh) = weight;

    return weight;
}

double WeightFunctor::textureWeight(HalfedgeHandle heh) const {
    // if the two faces along this edge are different color, set weight of this edge to almost 0
    // we encourage path go through the boundary of different colors
    if (PlushPatternGenerator::is_different_texture(m_mesh, m_mesh->edge_handle(heh))) {
        return 1e-9;
    } else {
        return 1;
    }
}

//double WeightFunctor::curvatureWeight(VertexHandle v1, VertexHandle v2) const {
//    double curvature1 = m_mesh->property(PlushPatternGenerator::maxCurvatureHandle, v1);
//    double curvature2 = m_mesh->property(PlushPatternGenerator::maxCurvatureHandle, v2);
//
//    // weight using curvature
//    // encourage path with +1 curvature, not 0
//    double weight = 1 / (1 + (abs(curvature1) + abs(curvature2)));
//    assert(0 <= weight && weight <= 1);
//    
//    return weight;
//}
double WeightFunctor::curvatureWeight(EdgeHandle eh) const {
    double weight = 1 / (1 + abs(m_mesh->calc_dihedral_angle(eh)));
    m_mesh->property(PlushPatternGenerator::curvature_weight_Handle, eh) = weight;

    return weight;
}

/// This measures the direction similarity between edge and corresponding bones
double WeightFunctor::skeletonWeight(EdgeHandle eh,
                                     VertexHandle v1,
                                     VertexHandle v2,
                                     TriMesh::Point p1,
                                     TriMesh::Point p2) const {
    Skeleton *skeleton = m_mesh->property(PlushPatternGenerator::skeletonHandle);
    if (!skeleton) {
        return 1;
    }
    
    // Calculate corresponding averaged bone direction for each vertex
    OpenMesh::Vec3d avgBoneDirection(0,0,0);
    int i = 0;
    for (std::vector<Bone>::iterator bone_it = skeleton->bones.begin(); bone_it != skeleton->bones.end(); bone_it++, i++) {
        OpenMesh::Vec3d dir = (bone_it->getA() - bone_it->getB()).normalize();
        double weight1 = m_mesh->property(PlushPatternGenerator::bonesWeightHandle, v1)[i];
        double weight2 = m_mesh->property(PlushPatternGenerator::bonesWeightHandle, v2)[i];
        
        assert(!isnan(weight1) && !isnan(weight2));
        avgBoneDirection += dir * (weight1 + weight2)/2;
    }
    avgBoneDirection.normalize();
    
    i = 0;
    double sumOfAngleDiffSqr = 0;
    for (std::vector<Bone>::iterator bone_it = skeleton->bones.begin(); bone_it != skeleton->bones.end(); bone_it++, i++) {
        OpenMesh::Vec3d dir = (bone_it->getA() - bone_it->getB()).normalize();
        double weight1 = m_mesh->property(PlushPatternGenerator::bonesWeightHandle, v1)[i];
        double weight2 = m_mesh->property(PlushPatternGenerator::bonesWeightHandle, v2)[i];
        
        assert(!isnan(weight1) && !isnan(weight2));
        double cosAngle = abs(dir | avgBoneDirection);
        sumOfAngleDiffSqr += cosAngle * (weight1 + weight2)/2;
    }
    // Prevent from floating precision error
    sumOfAngleDiffSqr = min(sumOfAngleDiffSqr, 1.0);
    
    // If there are many influential bones (std is low) for this edge, it's probably near the body part.
    // In such case, skeleton information is not reliable. So we give it a lower influential.
    double alpha = sumOfAngleDiffSqr;
    
    double cosAngle = ((p1-p2)|avgBoneDirection)/((p1-p2).norm()*avgBoneDirection.norm());
    double weight = 1 - alpha * fabs(cosAngle);
    m_mesh->property(PlushPatternGenerator::skeleton_direction_weight_Handle, eh) = weight;
    
    assert(0 <= weight && weight <= 1);
    return weight;
}

double WeightFunctor::smoothnessWeight(VertexHandle v1,
                                       VertexHandle v2,
                                       TriMesh::Point p1,
                                       TriMesh::Point p2) const {
    VertexHandle predecessor;
    if (m_predecessor_pmap != NULL) {
        predecessor = (*m_predecessor_pmap)[v1];
    } else if (m_mesh->is_valid_handle(m_prevV)) {
        predecessor = m_prevV;
    } else {
        assert("You should at least have valid m_predecessor_pmap or m_prevV to calculate smoothness weight");
    }
    if (predecessor == v1) {
        // we reach the begining, assign weight to 1 (max) to ensure that
        // when combining path, total weight would not increase
        return 1;
    } else {
        TriMesh::Point p0 = m_mesh->point(predecessor);
        double cosAngle = ((p1-p0) | (p1-p2) / ((p1-p0).norm() * (p1-p2).norm()));
        // cosAngle = [-1, 1], we need [0, 1]
        return (cosAngle + 1) / 2;
    }
}

double WeightFunctor::operator()(HalfedgeHandle heh) const {
    EdgeHandle eh = m_mesh->edge_handle(heh);
    VertexHandle v1 = m_mesh->from_vertex_handle(heh);
    VertexHandle v2 = m_mesh->to_vertex_handle(heh);

    if (m_mesh->is_valid_handle(m_currentV)) {
        assert (v1 == m_currentV || v2 == m_currentV);
        
        if (v2 == m_currentV) {
            v2 = v1;
            v1 = m_currentV;
        }
    }
    
    TriMesh::Point p1 = m_mesh->point(v1);
    TriMesh::Point p2 = m_mesh->point(v2);
    
    double coefficients = 0;
    double edgeWeight = 1;
    
    if (m_mesh->property(PlushPatternGenerator::edgeWeightHandle, eh) >= 0) {
        edgeWeight = m_mesh->property(PlushPatternGenerator::edgeWeightHandle, eh);
    } else {
        edgeWeight *= pow(distanceWeight(eh), m_distanceCoefficient);
//        edgeWeight *= pow(textureWeight(heh), m_textureCoefficient);
//        edgeWeight *= pow(curvatureWeight(v1, v2), m_curvatureCoefficient);
        edgeWeight *= pow(curvatureWeight(eh), m_curvatureCoefficient);
        edgeWeight *= pow(skeletonWeight(eh, v1, v2, p1, p2), m_skeletonCoefficient);
        
        m_mesh->property(PlushPatternGenerator::edgeWeightHandle, eh) = edgeWeight;
    }
    
    // re-calculate smoothness weight every time because it depends on path.
    // it can not be saved and reuse
    double pathWeight = 1;
    if (m_predecessor_pmap != NULL || m_mesh->is_valid_handle(m_prevV)) {
        pathWeight = pow(smoothnessWeight(v1, v2, p1, p2), m_pathCoefficient);
    }
    return edgeWeight * pathWeight;
}

/**
 Return the averaged sum of weights of given segment.
 @param <#parameter#>
 @return <#retval#>
 @retval <#meaning#>
 */
double WeightFunctor::operator()(std::vector<HalfedgeHandle> segment) {
    double sum = 0;
    for (HalfedgeHandle heh : segment) {
        m_currentV = m_mesh->from_vertex_handle(heh);
        
        sum += operator()(heh);
        
        m_prevV = m_currentV;
    }
    return sum/segment.size();
}