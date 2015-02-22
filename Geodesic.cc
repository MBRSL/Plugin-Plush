//
//  Geodesic.cc
//  OpenFlipper
//
//  Created by 饒亢 on 2015/1/13.
//
//

#include "PlushPlugin.hh"
#include "CGAL_Polyhedron_builder.hh"
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
        double curvature1 = mesh->property(PlushPlugin::maxCurvatureHandle, v1);
        double curvature2 = mesh->property(PlushPlugin::maxCurvatureHandle, v2);
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
        Skeleton *skeleton = mesh->property(PlushPlugin::skeletonHandle);

        // Calculate corresponding averaged bone direction for each vertex
        OpenMesh::Vec3d avgBoneDirection(0,0,0);
        int i = 0;
        for (std::vector<Bone>::iterator bone_it = skeleton->bones.begin(); bone_it != skeleton->bones.end(); bone_it++, i++) {
            OpenMesh::Vec3d dir = (bone_it->getA() - bone_it->getB()).normalize();
            double weight1 = mesh->property(PlushPlugin::bonesWeightHandle, v1)[i];
            double weight2 = mesh->property(PlushPlugin::bonesWeightHandle, v2)[i];
            
            avgBoneDirection += dir * (weight1);
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
            // we reach the begining, no smoothness penalty
            return 0;
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
                  const boost::iterator_property_map<std::vector<boost_vertex_descriptor>::iterator, VertexIdPropertyMap>
                    &predecessor_pmap,
                  const EdgeIdPropertyMap &edgeIndexMap) :
    mesh(mesh), P(P), predecessor_pmap(predecessor_pmap), edgeIndexMap(edgeIndexMap) {}
    
    double operator()(boost_edge_descriptor e) const {
        boost_vertex_descriptor boost_v1 = boost::source(e, P);
        boost_vertex_descriptor boost_v2 = boost::target(e, P);
        
        int index = get(edgeIndexMap, e);
        EdgeHandle eh = mesh->edge_handle(index);
        HalfedgeHandle he1 = mesh->halfedge_handle(eh, 0);
        HalfedgeHandle he2 = mesh->halfedge_handle(eh, 1);
        
        VertexHandle v1 = mesh->from_vertex_handle(he1);
        VertexHandle v2 = mesh->to_vertex_handle(he1);
        
        TriMesh::Point p1 = mesh->point(v1);
        TriMesh::Point p2 = mesh->point(v2);

        double edgeWeight = 0;
        if (mesh->property(PlushPlugin::edgeWeightHandle, eh) >= 0) {
            edgeWeight += mesh->property(PlushPlugin::edgeWeightHandle, eh);
        } else {
//            totalWeight += distanceWeight(p1, p2);
//            totalWeight += textureWeight(he1, he2);
//            totalWeight += curvatureWeight(v1, v2);
            edgeWeight += skeletonWeight(eh, v1, v2, p1, p2);
            mesh->property(PlushPlugin::edgeWeightHandle, eh) = edgeWeight;
        }

        // re-calculate smoothness weight every time because it depends on path.
        // it can not be saved and reuse
        double pathWeight = smoothnessWeight(boost_v1, boost_v2, p1, p2);
        
        return edgeWeight + pathWeight;
    }
};


void PlushPlugin::calcGeodesic(TriMesh *mesh,
                               Polyhedron &P,
                               std::map<int, Polyhedron::Vertex_handle> &verticesMapping,
                               VertexHandle sourceHandle,
                               IdList targetVertices)
{
    std::map<std::pair<VertexHandle, VertexHandle>, double> &geodesicDistance = mesh->property(geodesicDistanceHandle);
    std::map<std::pair<VertexHandle, VertexHandle>, IdList> &geodesicPath = mesh->property(geodesicPathHandle);
    
    // associate indices to the vertices & pre-calculated distance (if available)
    std::vector<double> distance(boost::num_vertices(P));

    int index = 0;
    boost_vertex_iterator vit, ve;
    boost_vertex_descriptor source, dest;
    for(boost::tie(vit,ve)=boost::vertices(P); vit!=ve; ++vit ){
        boost_vertex_descriptor vd = *vit;
        if (index == sourceHandle.idx()) {
            source = vd;
        }
        vd->id() = index;
        
        // initialize with existing distance if available
        VertexHandle targetHandle = mesh->vertex_handle(index);
        std::map<std::pair<VertexHandle, VertexHandle>, double>::iterator found =
            geodesicDistance.find(std::make_pair(sourceHandle, targetHandle));
        if (found != geodesicDistance.end()) {
            distance[index] = found->second;
        }
        index++;
    }
    
    // prepare property maps for dijkstra algorithm
    VertexIdPropertyMap vertex_index_pmap = get(CGAL::vertex_external_index, P);
    EdgeIdPropertyMap edge_index_pmap = get(boost::edge_external_index, P);
    
    std::vector<boost_vertex_descriptor> predecessor(boost::num_vertices(P));
    boost::iterator_property_map<std::vector<boost_vertex_descriptor>::iterator, VertexIdPropertyMap> predecessor_pmap(predecessor.begin(), vertex_index_pmap);
    boost::iterator_property_map<std::vector<double>::iterator, VertexIdPropertyMap> distance_pmap(distance.begin(), vertex_index_pmap);
    WeightFunctor weightFunctor(mesh,
                                P,
                                predecessor_pmap,
                                edge_index_pmap);
    auto weightmap = boost::make_function_property_map< boost_edge_descriptor,
                                                        double,
                                                        WeightFunctor > (weightFunctor);
    
    boost::dijkstra_shortest_paths(P, source,
                                   vertex_index_map(vertex_index_pmap)
                                   .weight_map(weightmap)
                                   .distance_map(distance_pmap)
                                   .predecessor_map(predecessor_pmap));
    
    // Save result
    for(boost::tie(vit,ve)=boost::vertices(P); vit!=ve; ++vit ){
        VertexHandle destHandle = mesh->vertex_handle((*vit)->id());
        
        std::pair<VertexHandle, VertexHandle> edgeSD = std::make_pair(sourceHandle, destHandle);
        std::pair<VertexHandle, VertexHandle> edgeDS = std::make_pair(destHandle, sourceHandle);
        // clear previous result
        geodesicDistance.erase(edgeSD);
        geodesicDistance.erase(edgeDS);
        geodesicPath.erase(edgeSD);
        geodesicPath.erase(edgeDS);

        geodesicDistance.insert(std::make_pair(edgeSD, distance_pmap[*vit]));
        geodesicDistance.insert(std::make_pair(edgeDS, distance_pmap[*vit]));
        
        IdList path;
        for(boost_vertex_iterator predecessor_it = vit; predecessor_pmap[*predecessor_it] != *predecessor_it;) {
            path.push_back((*predecessor_it)->id());
            predecessor_it = predecessor_pmap[*predecessor_it];
        }
        
        // empty if unreachable
        if (path.size() > 0) {
            // add source to path
            path.push_back(sourceHandle.idx());
        }
        geodesicPath.insert(std::make_pair(edgeSD, path));
        std::reverse(path.begin(), path.end());
        geodesicPath.insert(std::make_pair(edgeDS, path));
    }
}

bool PlushPlugin::calcSpanningTree(QString _jobId,
                                   int meshId,
                                   std::vector<std::pair<IdList, double> > &result,
                                   IdList selectedVertices) {
    BaseObjectData *obj;
    PluginFunctions::getObject(meshId, obj);
    if (!obj->dataType(DATA_TRIANGLE_MESH)) {
        emit log(LOGERR, QString("Not a valid TriMesh of object %1").arg(QString::number(meshId)));
        return false;
    }
    
    isJobCanceled = false;
    
    TriMesh *mesh = PluginFunctions::triMesh(obj);
    QString meshName = QFileInfo(obj->name()).baseName();
    emit setJobDescription(_jobId, QString("Calculating curvature: %1").arg(meshName));
    
    std::map<std::pair<VertexHandle, VertexHandle>, double> &geodesicDistance = mesh->property(geodesicDistanceHandle);
    std::map<std::pair<VertexHandle, VertexHandle>, IdList> &geodesicPath = mesh->property(geodesicPathHandle);
    
    // convert mesh from OpenFlipper to CGAL
    Polyhedron P;
    std::map<int, boost_vertex_descriptor> verticesMapping;
    CGAL_Polyhedron_helper::convert_OpenMesh_to_CGAL(mesh, P, verticesMapping);

    // create a new graph with selectedVerices, calculate all paths between them
    int iterations = 0;
    int totalIterations = (selectedVertices.size()*(selectedVertices.size()-1)/2);
    for (size_t i = 0; i < selectedVertices.size(); i++) {
        for (size_t j = i+1; j < selectedVertices.size(); j++) {
            if (isJobCanceled) {
                emit log(LOGINFO, "Geodesic calculation canceled.");
                return false;
            }
            
            int sourceIdx = selectedVertices[i];
            int destIdx = selectedVertices[j];
            VertexHandle sourceHandle = mesh->vertex_handle(sourceIdx);
            VertexHandle destHandle = mesh->vertex_handle(destIdx);
            
            // if result is not found, caculate now
            std::map<std::pair<VertexHandle, VertexHandle>, double>::iterator distanceFound =
                geodesicDistance.find(std::make_pair(sourceHandle, destHandle));
            std::map<std::pair<VertexHandle, VertexHandle>, IdList>::iterator pathFound =
                geodesicPath.find(std::make_pair(sourceHandle, destHandle));
            
            if (distanceFound == geodesicDistance.end() || pathFound == geodesicPath.end()) {
                calcGeodesic(mesh, P, verticesMapping, sourceHandle, selectedVertices);
                distanceFound = geodesicDistance.find(std::make_pair(sourceHandle, destHandle));
                pathFound = geodesicPath.find(std::make_pair(sourceHandle, destHandle));
            }
            double cost = distanceFound->second;
            IdList path = pathFound->second;
            
            // assign edge & weight
            result.push_back(std::make_pair(path, cost));
            
            // most of the time is spent on geodesic calculation
            iterations++;
            int status = (double)iterations/totalIterations * 100;
            emit setJobState(_jobId, status);
        }
    }
    return true;
}