//
//  Geodesic.cc
//  OpenFlipper
//
//  Created by 饒亢 on 2015/1/13.
//
//

#include "PlushPatternGenerator.hh"
#include "SuperDeform/Skeleton.hh"

#include <CGAL/boost/graph/dijkstra_shortest_paths.h>

#include <boost/property_map/function_property_map.hpp>

typedef boost::property_map<Polyhedron, boost::vertex_external_index_t>::type VertexIdPropertyMap;
typedef boost::property_map<Polyhedron, boost::edge_external_index_t>::type EdgeIdPropertyMap;

class WeightFunctor {
private:
    TriMesh *mesh;
    const Polyhedron P;
    const boost::iterator_property_map<std::vector<boost_vertex_descriptor>::iterator, VertexIdPropertyMap>
        predecessor_pmap;
    const EdgeIdPropertyMap edgeIndexMap;
    boost_vertex_descriptor *currentV;

    double distanceWeight(TriMesh::Point p1, TriMesh::Point p2) const {
        return (p1-p2).norm();
    }
    
    double textureWeight(HalfedgeHandle he1, HalfedgeHandle he2) const {
        // if the two faces along this edge are different color, set weight of this edge to almost 0
        // we encourage path go through the boundary of different colors
        if (!mesh->is_boundary(he1) && !mesh->is_boundary(he2)
            &&  mesh->color(mesh->face_handle(he1)) != mesh->color(mesh->face_handle(he2))) {
            return 1e-9;
        } else {
            return 1;
        }
    }
    
    double curvatureWeight(VertexHandle v1, VertexHandle v2) const {
        double curvature1 = mesh->property(PlushPatternGenerator::maxCurvatureHandle, v1);
        double curvature2 = mesh->property(PlushPatternGenerator::maxCurvatureHandle, v2);
        // clamp curvature to [-1, 1]
        curvature1 = fmax(fmin(curvature1, 1), -1);
        curvature2 = fmax(fmin(curvature2, 1), -1);

        // weight using curvature
        // encourage path with +1/-1 curvature, not 0
        return 1 - abs(curvature1 + curvature2)/2;
    }
    
    double skeletonWeight(EdgeHandle eh,
                          VertexHandle v1,
                          VertexHandle v2,
                          TriMesh::Point p1,
                          TriMesh::Point p2) const {
        Skeleton *skeleton = mesh->property(PlushPatternGenerator::skeletonHandle);

        // Calculate corresponding averaged bone direction for each vertex
        OpenMesh::Vec3d avgBoneDirection(0,0,0);
        int i = 0;
        for (std::vector<Bone>::iterator bone_it = skeleton->bones.begin(); bone_it != skeleton->bones.end(); bone_it++, i++) {
            OpenMesh::Vec3d dir = (bone_it->getA() - bone_it->getB()).normalize();
            double weight1 = mesh->property(PlushPatternGenerator::bonesWeightHandle, v1)[i];
            double weight2 = mesh->property(PlushPatternGenerator::bonesWeightHandle, v2)[i];
            
            avgBoneDirection += dir * (weight1 + weight2)/2;
        }
        
        double weight = 1-fabs(((p1-p2)|avgBoneDirection)/((p1-p2).norm()*avgBoneDirection.norm()));
        return weight;
    }
    
    double smoothnessWeight(boost_vertex_descriptor boost_v1,
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
    
public:
    WeightFunctor(TriMesh *mesh,
                  const Polyhedron &P,
                  boost_vertex_descriptor *currentV,
                  const boost::iterator_property_map<std::vector<boost_vertex_descriptor>::iterator, VertexIdPropertyMap>
                    &predecessor_pmap,
                  const EdgeIdPropertyMap &edgeIndexMap) :
    mesh(mesh), P(P), predecessor_pmap(predecessor_pmap), edgeIndexMap(edgeIndexMap), currentV(currentV) {}
    
    double operator()(boost_edge_descriptor e) const {
        
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

        double edgeWeight = 0;
        double pathWeight = 0;
        if (mesh->property(PlushPatternGenerator::edgeWeightHandle, eh) >= 0) {
            edgeWeight += mesh->property(PlushPatternGenerator::edgeWeightHandle, eh);
        } else {
//            edgeWeight += distanceWeight(p1, p2);
//            edgeWeight += textureWeight(he1, he2);
//            edgeWeight += curvatureWeight(v1, v2);
            edgeWeight += skeletonWeight(eh, v1, v2, p1, p2);
            mesh->property(PlushPatternGenerator::edgeWeightHandle, eh) = edgeWeight;
        }

        // re-calculate smoothness weight every time because it depends on path.
        // it can not be saved and reuse
        pathWeight = smoothnessWeight(boost_v1, boost_v2, p1, p2);
        
        return edgeWeight + pathWeight;
    }
};

class Dijkstra_visitor : public boost::default_dijkstra_visitor {
private:
    boost_vertex_descriptor *currentV;
public:
    Dijkstra_visitor (boost_vertex_descriptor *currentV) : currentV(currentV) {}
    
    template <class Vertex, class Graph>
    void examine_vertex(Vertex u, Graph &g) {
        *currentV = u;
    }
};

void PlushPatternGenerator::calcGeodesic(VertexHandle sourceHandle, std::vector<int> targetVertices)
{
    // Prepare property maps for dijkstra algorithm
    VertexIdPropertyMap vertex_index_pmap = get(boost::vertex_external_index, m_polyhedron);
    EdgeIdPropertyMap edge_index_pmap = get(boost::edge_external_index, m_polyhedron);
    
    std::vector<boost_vertex_descriptor> predecessor(boost::num_vertices(m_polyhedron));
    boost::iterator_property_map<std::vector<boost_vertex_descriptor>::iterator, VertexIdPropertyMap> predecessor_pmap(predecessor.begin(), vertex_index_pmap);

    std::vector<double> distance(boost::num_vertices(m_polyhedron));
    boost::iterator_property_map<std::vector<double>::iterator, VertexIdPropertyMap> distance_pmap(distance.begin(), vertex_index_pmap);
    
    boost_vertex_descriptor source = m_verticesMapping[sourceHandle.idx()];
    
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
    std::map<std::pair<VertexHandle, VertexHandle>, std::vector<int> > &geodesicPath = m_mesh->property(geodesicPathHandle);
    for (size_t i = 0; i < targetVertices.size(); i++) {
        int id = targetVertices[i];
        VertexHandle destHandle = m_mesh->vertex_handle(id);
        if (destHandle == sourceHandle) {
            continue;
        }
    
        std::pair<VertexHandle, VertexHandle> edgeSD = std::make_pair(sourceHandle, destHandle);
        std::pair<VertexHandle, VertexHandle> edgeDS = std::make_pair(destHandle, sourceHandle);

        // clear previous result
        geodesicDistance.erase(edgeSD);
        geodesicDistance.erase(edgeDS);
        geodesicPath.erase(edgeSD);
        geodesicPath.erase(edgeDS);

        // back tracking
        std::vector<int> path;
        for(boost_vertex_descriptor predecessor = m_verticesMapping[id]; predecessor_pmap[predecessor] != predecessor;) {
            path.push_back(predecessor->id());
            predecessor = predecessor_pmap[predecessor];
        }
        
        // empty if unreachable
        if (path.size() > 0) {
            // add source to path
            path.push_back(sourceHandle.idx());
        }
        
        // write into property
        geodesicDistance.insert(std::make_pair(edgeSD, distance_pmap[m_verticesMapping[id]]));
        geodesicDistance.insert(std::make_pair(edgeDS, distance_pmap[m_verticesMapping[id]]));

        geodesicPath.insert(std::make_pair(edgeDS, path));
        std::reverse(path.begin(), path.end());
        geodesicPath.insert(std::make_pair(edgeSD, path));
    }
}
