//
//  Geodesic.cc
//  OpenFlipper
//
//  Created by 饒亢 on 2015/1/13.
//
//

#include "PlushPatternGenerator.hh"
#include "WeightFunctor.hh"
#include "SuperDeform/Skeleton.hh"

#include <CGAL/boost/graph/dijkstra_shortest_paths.h>

#include <boost/property_map/function_property_map.hpp>

#include <QTextStream>

typedef boost::property_map<Polyhedron, boost::vertex_external_index_t>::type VertexIdPropertyMap;
typedef boost::property_map<Polyhedron, boost::edge_external_index_t>::type EdgeIdPropertyMap;

WeightFunctor::WeightFunctor(TriMesh *mesh,
              const Polyhedron &P,
              boost_vertex_descriptor *currentV,
              const boost::iterator_property_map<std::vector<boost_vertex_descriptor>::iterator, VertexIdPropertyMap>
                 &predecessor_pmap,
              const EdgeIdPropertyMap &edgeIndexMap) :
mesh(mesh), P(P), predecessor_pmap(predecessor_pmap), edgeIndexMap(edgeIndexMap), currentV(currentV), avgEdgeLength(0) {
    for (EdgeIter e_it = mesh->edges_begin(); e_it != mesh->edges_end(); e_it++) {
        HalfedgeHandle he = mesh->halfedge_handle(*e_it, 0);
        TriMesh::Point p1 = mesh->point(mesh->from_vertex_handle(he));
        TriMesh::Point p2 = mesh->point(mesh->to_vertex_handle(he));
        
        avgEdgeLength += (p1-p2).norm();
    }
    avgEdgeLength /= mesh->n_edges();
}

double WeightFunctor::distanceWeight(TriMesh::Point p1, TriMesh::Point p2) const {
    return (p1-p2).norm()/avgEdgeLength;
}
    
double WeightFunctor::textureWeight(HalfedgeHandle he1, HalfedgeHandle he2) const {
    // if the two faces along this edge are different color, set weight of this edge to almost 0
    // we encourage path go through the boundary of different colors
    if (!mesh->is_boundary(he1) && !mesh->is_boundary(he2)
        &&  mesh->color(mesh->face_handle(he1)) != mesh->color(mesh->face_handle(he2))) {
        return 1e-9;
    } else {
        return 1;
    }
}

double WeightFunctor::curvatureWeight(VertexHandle v1, VertexHandle v2) const {
    double curvature1 = mesh->property(PlushPatternGenerator::maxCurvatureHandle, v1);
    double curvature2 = mesh->property(PlushPatternGenerator::maxCurvatureHandle, v2);
    // clamp curvature to [-1, 1]
    curvature1 = fmax(fmin(curvature1, 1), -1);
    curvature2 = fmax(fmin(curvature2, 1), -1);

    // weight using curvature
    // encourage path with +1/-1 curvature, not 0
    return 1 - abs(curvature1 + curvature2)/2;
}

/// This measures the direction similarity between edge and corresponding bones
double WeightFunctor::skeletonWeight(VertexHandle v1,
                                     VertexHandle v2,
                                     TriMesh::Point p1,
                                     TriMesh::Point p2) const {
    Skeleton *skeleton = mesh->property(PlushPatternGenerator::skeletonHandle);

    // Calculate corresponding averaged bone direction for each vertex
    OpenMesh::Vec3d avgBoneDirection(0,0,0);
    double sumOfSqr1 = 0, sumOfSqr2 = 0;
    int i = 0;
    for (std::vector<Bone>::iterator bone_it = skeleton->bones.begin(); bone_it != skeleton->bones.end(); bone_it++, i++) {
        OpenMesh::Vec3d dir = (bone_it->getA() - bone_it->getB()).normalize();
        double weight1 = mesh->property(PlushPatternGenerator::bonesWeightHandle, v1)[i];
        double weight2 = mesh->property(PlushPatternGenerator::bonesWeightHandle, v2)[i];
        
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
    double weight = influential * (1-fabs(cosAngle)) + pow(1-influential, 3);
    return weight;
}

double WeightFunctor::smoothnessWeight(boost_vertex_descriptor boost_v1,
                        boost_vertex_descriptor boost_v2,
                        TriMesh::Point p1,
                        TriMesh::Point p2) const {
    boost_vertex_descriptor predecessor = predecessor_pmap[boost_v1];
    if (predecessor == boost_v1) {
        // we reach the begining, assign weight to 1 (max) to ensure that
        // when combining path, total weight would not increase
        return 1;
    } else {
        auto boost_p0 = predecessor->point();
        TriMesh::Point p0(boost_p0[0],boost_p0[1],boost_p0[2]);
        double cosAngle = ((p1-p0) | (p1-p2) / ((p1-p0).norm() * (p1-p2).norm()));
        // cosAngle = [-1, 1], we need [0, 1]
        return (cosAngle + 1) / 2;
    }
}
    
double WeightFunctor::operator()(boost_edge_descriptor e) const {
    
    int index = get(edgeIndexMap, e);
    EdgeHandle eh = mesh->edge_handle(index);

    boost_vertex_descriptor boost_v1 = boost::source(e, P);
    boost_vertex_descriptor boost_v2 = boost::target(e, P);
    
    assert (boost_v1 == *currentV || boost_v2 == *currentV);
    
    if (boost_v2 == *currentV) {
        boost_v2 = boost_v1;
        boost_v1 = *currentV;
    }

    HalfedgeHandle he1 = mesh->halfedge_handle(eh, 0);
    HalfedgeHandle he2 = mesh->halfedge_handle(eh, 1);
    
    VertexHandle v1 = mesh->from_vertex_handle(he1);
    VertexHandle v2 = mesh->to_vertex_handle(he1);
    
    // make sure the direction is correct with boost_vertex_descriptor
    if (v1.idx() != (int)boost_v1->id()) {
        VertexHandle tmp = v1;
        v1 = v2;
        v2 = tmp;
    }
    
    TriMesh::Point p1 = mesh->point(v1);
    TriMesh::Point p2 = mesh->point(v2);

    double coefficients = 0;
    double edgeWeight = 0;

    double distanceCoefficient = 1;
    double textureCoefficient = 1;
    double curvatureCoefficient = 1;
    double skeletonCoefficient = 1;
    coefficients += distanceCoefficient + textureCoefficient + curvatureCoefficient + skeletonCoefficient;

    if (mesh->property(PlushPatternGenerator::edgeWeightHandle, eh) >= 0) {
        edgeWeight += mesh->property(PlushPatternGenerator::edgeWeightHandle, eh);
    } else {
        edgeWeight += distanceCoefficient * distanceWeight(p1, p2);
        edgeWeight += textureCoefficient * textureWeight(he1, he2);
        edgeWeight += curvatureCoefficient * curvatureWeight(v1, v2);
        edgeWeight += skeletonCoefficient * skeletonWeight(v1, v2, p1, p2);
        
        mesh->property(PlushPatternGenerator::edgeWeightHandle, eh) = edgeWeight;
    }

    // re-calculate smoothness weight every time because it depends on path.
    // it can not be saved and reuse
    double pathCoefficient = 1;
    double pathWeight = pathCoefficient * smoothnessWeight(boost_v1, boost_v2, p1, p2);
    coefficients += pathCoefficient;
    
    return (edgeWeight + pathWeight) / coefficients;
}

void PlushPatternGenerator::calcGeodesic(std::vector<VertexHandle> targetVertices)
{
    // Prepare property maps for dijkstra algorithm
    VertexIdPropertyMap vertex_index_pmap = get(boost::vertex_external_index, m_polyhedron);
    EdgeIdPropertyMap edge_index_pmap = get(boost::edge_external_index, m_polyhedron);

    isJobCanceled = false;

    for (size_t i = 0; i < targetVertices.size(); i++) {
        if (isJobCanceled) {
            emit log(LOGINFO, "Geodesic calculation canceled.");
            return;
        }

        VertexHandle sourceHandle = targetVertices[i];
        
        std::vector<boost_vertex_descriptor> predecessor(boost::num_vertices(m_polyhedron));
        boost::iterator_property_map<std::vector<boost_vertex_descriptor>::iterator, VertexIdPropertyMap> predecessor_pmap(predecessor.begin(), vertex_index_pmap);

        std::vector<double> distance(boost::num_vertices(m_polyhedron));
        boost::iterator_property_map<std::vector<double>::iterator, VertexIdPropertyMap> distance_pmap(distance.begin(), vertex_index_pmap);
        
        boost_vertex_descriptor source = m_verticesMapping[targetVertices[i].idx()];
        
        // HACK: This is a local variable "shared" by both (functional)weightmap & visitor
        //       so that we can first assign it in visitor and then access it in weightmap
        boost_vertex_descriptor currentV = source;
        
        WeightFunctor weightFunctor(m_mesh,
                                    m_polyhedron,
                                    &currentV,
                                    predecessor_pmap,
                                    edge_index_pmap);
        auto weightmap = boost::make_function_property_map< boost_edge_descriptor,
                                                            double,
                                                            WeightFunctor > (weightFunctor);
        Dijkstra_visitor visitor(&currentV);
        
        // Calculate geodesic using Dijkstra
        boost::dijkstra_shortest_paths(m_polyhedron, source,
                                       vertex_index_map(vertex_index_pmap)
                                       .weight_map(weightmap)
                                       .distance_map(distance_pmap)
                                       .predecessor_map(predecessor_pmap)
                                       .visitor(visitor));
        
        // Save result into property
        std::map<std::pair<VertexHandle, VertexHandle>, double> &geodesicDistance = m_mesh->property(geodesicDistanceHandle);
        std::map<std::pair<VertexHandle, VertexHandle>, std::vector<VertexHandle> > &geodesicPath = m_mesh->property(geodesicPathHandle);
        for (size_t j = 0; j < targetVertices.size(); j++) {
        VertexHandle destHandle = targetVertices[j];
    //    for (VertexIter v_it = m_mesh->vertices_begin(); v_it != m_mesh->vertices_end(); v_it++) {
    //        VertexHandle destHandle = *v_it;
    //        int id = v_it->idx();
            if (destHandle == sourceHandle) {
                continue;
            }

            // Directions. from source to destination, or from destination to source
            std::pair<VertexHandle, VertexHandle> edgeSD = std::make_pair(sourceHandle, destHandle);
            std::pair<VertexHandle, VertexHandle> edgeDS = std::make_pair(destHandle, sourceHandle);
            
            // Back tracking
            std::vector<VertexHandle> path;
            for(boost_vertex_descriptor predecessor = m_verticesMapping[destHandle.idx()]; predecessor_pmap[predecessor] != predecessor;) {
                path.push_back(m_mesh->vertex_handle(predecessor->id()));
                predecessor = predecessor_pmap[predecessor];
            }
            
            if (path.size() <= 0) {
                // Unreachable
                return;
            } else {
                // Add source to path
                path.push_back(sourceHandle);
            }
            // It's reversed using back tracking, flip it back.
            std::reverse(path.begin(), path.end());

            // DEBUGGING
            // The cost is sometimes non symmetric, this should not happen.
    //        if (geodesicDistance.find(edgeSD) != geodesicDistance.end()) {
    //            double prevDist = geodesicDistance[edgeSD];
    //            if (abs(prevDist - distance_pmap[m_verticesMapping[id]]) > 1e-5) {
    //                double newDist = distance_pmap[m_verticesMapping[id]];
    //            }
    //        }
            
            // TODO:
            // The true distance from A to B should be identical with distance from B to A
            // Currently we only choose the direction with smaller cost
            double cost = distance_pmap[m_verticesMapping[destHandle.idx()]];
            
            // First, check if opposite direction is already calculated
            std::map<std::pair<VertexHandle, VertexHandle>, double>::iterator cost_reverse_it = geodesicDistance.find(edgeDS);
            if (cost_reverse_it != geodesicDistance.end()) {
                double cost_reverse = cost_reverse_it->second;

                if (cost_reverse < cost) {
                    // Copy cost and path from reverse direction
                    geodesicDistance.insert(std::make_pair(edgeSD, cost_reverse));
                    
                    std::map<std::pair<VertexHandle, VertexHandle>, std::vector<VertexHandle> >::iterator path_reverse_it = geodesicPath.find(edgeDS);
                    // Path should exist because they are calculated together
                    assert(path_reverse_it != geodesicPath.end());
                    geodesicPath.insert(std::make_pair(edgeSD, path_reverse_it->second));
                } else {
                    // Replace reverse direction with current cost and path
                    geodesicDistance.erase(edgeDS);
                    geodesicPath.erase(edgeDS);
                    
                    geodesicDistance.insert(std::make_pair(edgeSD, cost));
                    geodesicDistance.insert(std::make_pair(edgeDS, cost));
                    geodesicPath.insert(std::make_pair(edgeSD, path));
                    std::reverse(path.begin(), path.end());
                    geodesicPath.insert(std::make_pair(edgeDS, path));
                }
            }
            // If not, just save it into property
            else {
                geodesicDistance.insert(std::make_pair(edgeSD, cost));
                geodesicPath.insert(std::make_pair(edgeSD, path));
            }
        }
        // most of the time is spent on geodesic calculation
        int status = (double)(i+1)/targetVertices.size() * 100;
        emit setJobState(status);
    }
}

bool PlushPatternGenerator::saveGeodesic(std::vector<VertexHandle> selectedVertices) {
    QString geodesicFilename = m_meshName + ".geodesic";
    QFile file(geodesicFilename);
    if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) {
        emit log(LOGERR, QString("Error opening file: %1").arg(geodesicFilename));
        return false;
    }
    
    std::map<std::pair<VertexHandle, VertexHandle>, double> &geodesicDistance = m_mesh->property(geodesicDistanceHandle);
    std::map<std::pair<VertexHandle, VertexHandle>, std::vector<VertexHandle> > &geodesicPath = m_mesh->property(geodesicPathHandle);

    QTextStream out(&file);
    for (size_t i = 0; i < selectedVertices.size(); i++) {
        for (size_t j = i+1; j < selectedVertices.size(); j++) {
            VertexHandle sourceHandle = selectedVertices[i];
            VertexHandle destHandle = selectedVertices[j];

            // Directions. from source to destination, or from destination to source
            std::pair<VertexHandle, VertexHandle> edgeSD = std::make_pair(sourceHandle, destHandle);
            std::pair<VertexHandle, VertexHandle> edgeDS = std::make_pair(destHandle, sourceHandle);
            
            // source / dest / cost
            out << sourceHandle.idx() << " " << destHandle.idx() << " " << geodesicDistance[edgeSD] << endl;
            
            // num of Vertices on path/ v1 v2 ...
            std::vector<VertexHandle> path = geodesicPath[edgeSD];
            out << path.size();
            for (size_t pi = 0; pi < path.size(); pi++) {
                out << " " << path[pi].idx();
            }
            out << endl;
        }
    }
    file.close();
    return true;
}

bool PlushPatternGenerator::loadGeodesic() {
    QString geodesicFilename = m_meshName + ".geodesic";
    
    QFile file(geodesicFilename);
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
        emit log(LOGERR, QString("Error opening file: %1").arg(geodesicFilename));
        return false;
    }
    
    std::map<std::pair<VertexHandle, VertexHandle>, double> &geodesicDistance = m_mesh->property(geodesicDistanceHandle);
    std::map<std::pair<VertexHandle, VertexHandle>, std::vector<VertexHandle> > &geodesicPath = m_mesh->property(geodesicPathHandle);
    
    QTextStream fin(&file);
    while (!fin.atEnd()) {
        int sourceIdx, destIdx;
        double cost;
        fin >> sourceIdx >> destIdx >> cost;
        
        VertexHandle sourceHandle = m_mesh->vertex_handle(sourceIdx);
        VertexHandle destHandle = m_mesh->vertex_handle(destIdx);

        // Directions. from source to destination, or from destination to source
        std::pair<VertexHandle, VertexHandle> edgeSD = std::make_pair(sourceHandle, destHandle);
        std::pair<VertexHandle, VertexHandle> edgeDS = std::make_pair(destHandle, sourceHandle);
        
        geodesicDistance.erase(edgeSD);
        geodesicDistance.erase(edgeDS);
        geodesicDistance.emplace(edgeSD, cost);
        geodesicDistance.emplace(edgeDS, cost);

        // consume endl
        fin.readLine();
        QString pathStr = fin.readLine();
        QStringList vertices = pathStr.split(" ");
        // The first one is number of vertices on path
        int num = vertices[0].toInt();
        
        std::vector<VertexHandle> path;
        for (int i = 1; i <= num; i++) {
            path.push_back(m_mesh->vertex_handle(vertices[i].toInt()));
        }
        
        geodesicPath.erase(edgeSD);
        geodesicPath.erase(edgeDS);
        geodesicPath.emplace(edgeSD, path);
        reverse(path.begin(), path.end());
        geodesicPath.emplace(edgeDS, path);
    }
    file.close();
    return true;
}
