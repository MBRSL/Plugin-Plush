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
#include <boost/graph/kruskal_min_spanning_tree.hpp>

#include <boost/property_map/function_property_map.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/undirected_dfs.hpp>
#include <boost/cstdlib.hpp>

typedef boost::graph_traits<Polyhedron> GraphTraits;
typedef GraphTraits::vertex_descriptor boost_vertex_descriptor;
typedef GraphTraits::vertex_iterator boost_vertex_iterator;
typedef GraphTraits::edge_descriptor boost_edge_descriptor;
typedef GraphTraits::edge_iterator boost_edge_iterator;
typedef GraphTraits::halfedge_descriptor boost_halfedge_descriptor;
typedef GraphTraits::halfedge_iterator boost_halfedge_iterator;

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
        std::vector<OpenMesh::Vec3d> joints = mesh->property(PlushPlugin::skeletonJointsHandle);
        std::vector<Bone> bones = mesh->property(PlushPlugin::skeletonBonesHandle);

        // Calculate corresponding averaged bone direction for each vertex
        OpenMesh::Vec3d avgBoneDirection(0,0,0);
        int i = 0;
        for (std::vector<Bone>::iterator bone_it = bones.begin(); bone_it != bones.end(); bone_it++, i++) {
            OpenMesh::Vec3d dir = (bone_it->getA() - bone_it->getB()).normalize();
            double weight1 = mesh->property(PlushPlugin::skeletonBonesWeightHandle, v1)[i];
            double weight2 = mesh->property(PlushPlugin::skeletonBonesWeightHandle, v2)[i];
            
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

//double PlushPlugin::getEdgeWeight(TriMesh *mesh, int v1No, int v2No) {
//    VertexHandle v1 = mesh->vertex_handle(v1No);
//    VertexHandle v2 = mesh->vertex_handle(v2No);
//    return getEdgeWeight(mesh, v1, v2);
//}
//
//double PlushPlugin::getEdgeWeight(TriMesh *mesh, VertexHandle v1, VertexHandle v2) {
//    EdgeHandle eh;
//    getEdge(mesh, eh, v1, v2);
//    return getEdgeWeight(mesh, eh);
//}

//double PlushPlugin::getEdgeWeight(TriMesh *mesh, EdgeHandle eh) {
//    if (mesh->property(edgeWeightHandle, eh) >= 0) {
//        return mesh->property(edgeWeightHandle, eh);
//    } else {
//        HalfedgeHandle he1 = mesh->halfedge_handle(eh, 0);
//        HalfedgeHandle he2 = mesh->halfedge_handle(eh, 1);
//        
//        // if the two faces along this edge are different color, set weight of this edge to almost 0
//        // we encourage path go through the boundary of different colors
//        //    if (!mesh->is_boundary(he1) && !mesh->is_boundary(he2)
//        //        &&  mesh->color(mesh->face_handle(he1)) != mesh->color(mesh->face_handle(he2))) {
//        //        return 0;
//        //    } else {
//        VertexHandle v1 = mesh->from_vertex_handle(he1);
//        VertexHandle v2 = mesh->to_vertex_handle(he1);
//        TriMesh::Point p1 = mesh->point(v1);
//        TriMesh::Point p2 = mesh->point(v2);
//        OpenMesh::Vec3d vector_v = (p1-p2);
//        
//        std::vector<OpenMesh::Vec3d> joints = mesh->property(skeletonJointsHandle);
//        std::vector<Bone> bones = mesh->property(skeletonBonesHandle);
//        //    int closestIdx1 = mesh->property(closestSkeletonPointHandle, v1);
//        //    int closestIdx2 = mesh->property(closestSkeletonPointHandle, v2);
//        
//        // Calculate corresponding averaged bone direction for each vertex
//        
//        
//        // Search for adjacency points on skeleton
//        OpenMesh::Vec3d avgBoneDirection(0,0,0);
//        int i = 0;
//        for (std::vector<Bone>::iterator bone_it = bones.begin(); bone_it != bones.end(); bone_it++, i++) {
//            OpenMesh::Vec3d dir = (bone_it->getA() - bone_it->getB()).normalize();
//            double weight1 = mesh->property(skeletonBonesWeightHandle, v1)[i];
//            double weight2 = mesh->property(skeletonBonesWeightHandle, v2)[i];
//            
//            avgBoneDirection += dir * (weight1);
//        }
//        double weight = 1-fabs((vector_v|avgBoneDirection)/(vector_v.norm()*avgBoneDirection.norm()));
//        mesh->property(edgeWeightHandle, eh) = weight;
//        return weight;
//        //        return std::sqrt((p1-p2).norm());
//        //            double curvature1 = mesh->property(maxCurvatureHandle, v1);
//        //            double curvature2 = mesh->property(maxCurvatureHandle, v2);
//        //            // clamp curvature to [-1, 1]
//        //            curvature1 = fmax(fmin(curvature1, 1), -1);
//        //            curvature2 = fmax(fmin(curvature2, 1), -1);
//        //
//        //            // weight using curvature
//        //            // encourage path with +1/-1 curvature, not 0
//        //            return 1 - abs(curvature1 + curvature2)/2;
//        //    }
//    }
//}

void PlushPlugin::calcGeodesic(TriMesh *mesh, VertexHandle sourceHandle)
{
    // convert mesh from OpenFlipper to CGAL
    Polyhedron P;
    translate_mesh_from_OpenMesh_to_CGAL(mesh, &P);
    
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
            geodesicDistance->find(std::make_pair(sourceHandle, targetHandle));
        if (found != geodesicDistance->end()) {
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
//                                   .weight_map(weight_pmap)
                                   .weight_map(weightmap)
                                   .distance_map(distance_pmap)
                                   .predecessor_map(predecessor_pmap));
    
    // Save result
    for(boost::tie(vit,ve)=boost::vertices(P); vit!=ve; ++vit ){
        VertexHandle destHandle = mesh->vertex_handle((*vit)->id());
        
        std::pair<VertexHandle, VertexHandle> edgeSD = std::make_pair(sourceHandle, destHandle);
        std::pair<VertexHandle, VertexHandle> edgeDS = std::make_pair(destHandle, sourceHandle);
        // clear previous result
        geodesicDistance->erase(edgeSD);
        geodesicDistance->erase(edgeDS);
        geodesicPath->erase(edgeSD);
        geodesicPath->erase(edgeDS);

        geodesicDistance->insert(std::make_pair(edgeSD, distance_pmap[*vit]));
        geodesicDistance->insert(std::make_pair(edgeDS, distance_pmap[*vit]));
        
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
        geodesicPath->insert(std::make_pair(edgeSD, path));
        std::reverse(path.begin(), path.end());
        geodesicPath->insert(std::make_pair(edgeDS, path));
    }
}

template <class Vertex, class Graph>
struct DFS_Visitor : public boost::dfs_visitor<>
{
    Vertex source, dest;
    std::vector<Vertex> *path;
    bool isTerminated;
    
    DFS_Visitor(std::vector<Vertex> *path, Vertex source, Vertex dest) : source(source), dest(dest), path(path) {
        isTerminated = false;
    }
    
    void discover_vertex(Vertex u, Graph g) {
        if (!isTerminated) {
            path->push_back(u);
        }
        
        if (u == dest) {
            isTerminated = true;
        }
    }
    void finish_vertex(Vertex u, Graph g) {
        if (!isTerminated) {
            path->pop_back();
        }
    }
};

//void PlushPlugin::findPath(TriMesh *mesh, std::set<EdgeHandle> &spanningTree, std::vector<VertexHandle> &result, VertexHandle sourceHandle, VertexHandle destHandle) {
//    typedef boost::adjacency_list< boost::vecS, boost::vecS, boost::undirectedS,
//    boost::property<boost::vertex_index_t, int>,
//    boost::property<boost::edge_color_t, boost::default_color_type> > graph_t;
//    typedef boost::graph_traits<graph_t>::vertex_descriptor vertex_t;
//    
//    int verticesCount = 0;
//    for (VertexIter v_it = mesh->vertices_begin(); v_it != mesh->vertices_end(); v_it++, verticesCount++);
//    
//    graph_t g(verticesCount);
//    
//    for (std::set<EdgeHandle>::iterator e_it = spanningTree.begin(); e_it != spanningTree.end(); e_it++) {
//        HalfedgeHandle heh = mesh->halfedge_handle(*e_it, 0);
//        VertexHandle v1 = mesh->from_vertex_handle(heh);
//        VertexHandle v2 = mesh->to_vertex_handle(heh);
//        
//        boost::add_edge(v1.idx(), v2.idx(), g);
//    }
//    
//    std::vector<vertex_t> path;
//    vertex_t source = boost::vertex(sourceHandle.idx(), g);
//    vertex_t dest = boost::vertex(destHandle.idx(), g);
//    DFS_Visitor<vertex_t, graph_t> vis(&path, source, dest);
//    undirected_dfs(g,
//                   boost::root_vertex(vertex_t(sourceHandle.idx()))
//                   .visitor(vis)
//                   .edge_color_map(get(boost::edge_color, g)));
//    for (std::vector<vertex_t>::iterator v_it = vis.path->begin(); v_it != vis.path->end(); v_it++) {
//        int id = boost::get(boost::vertex_index, g, *v_it);
//        VertexHandle vh = mesh->vertex_handle(id);
//        result.push_back(vh);
//    }
//}
bool PlushPlugin::calcSpanningTree(QString _jobId, int meshId, std::set<EdgeHandle> &result, IdList selectedVertices) {
    return calcSpanningTree(_jobId, meshId, result, selectedVertices, (selectedVertices.size()-1+selectedVertices.size()-2)*(selectedVertices.size()-1)/2);
}

bool PlushPlugin::calcSpanningTree(QString _jobId, int meshId, std::set<EdgeHandle> &result, IdList selectedVertices, int edges) {
    typedef boost::adjacency_list < boost::vecS, boost::vecS, boost::undirectedS,
        boost::no_property, boost::property < boost::edge_weight_t, int > > Graph;
    typedef boost::graph_traits < Graph >::edge_descriptor Edge;
    typedef boost::graph_traits < Graph >::vertex_descriptor Vertex;
//    typedef std::pair<int, int> E;
    
//    std::vector<E> edges;
//    std::vector<double> weights;
    
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
    
    std::vector<std::pair<IdList, double> > distance;
    
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
            std::map<std::pair<VertexHandle, VertexHandle>, double>::iterator found =
            geodesicDistance->find(std::make_pair(sourceHandle, destHandle));
            if (found == geodesicDistance->end()) {
                calcGeodesic(mesh, sourceHandle);
                found = geodesicDistance->find(std::make_pair(sourceHandle, destHandle));
            }
            double cost = found->second;
            
            std::map<std::pair<VertexHandle, VertexHandle>, IdList>::iterator found2 =
            geodesicPath->find(std::make_pair(sourceHandle, destHandle));
            if (found2 == geodesicPath->end()) {
                calcGeodesic(mesh, sourceHandle);
                found2 = geodesicPath->find(std::make_pair(sourceHandle, destHandle));
            }
            IdList path = found2->second;
            
            // assign edge & weight
            distance.push_back(std::make_pair(path, cost/path.size()));
            
            // most of the time is spent on geodesic calculation
            iterations++;
            int status = (double)iterations/totalIterations * 100;
            emit setJobState(_jobId, status);
        }
    }
    
    struct Comparator {
        bool operator() (std::pair<IdList, double> a,
                         std::pair<IdList, double> b) {
            return a.second < b.second;
        }
    } comparator;
    
    std::sort(distance.begin(), distance.end(), comparator);
    
    int count = 0;
    for (std::vector<std::pair<IdList, double> >::iterator it = distance.begin(); count < edges && it != distance.end(); it++, count++) {
        IdList path = it->first;
        IdList::iterator vIdx_it = path.begin();
        int prevIdx = *vIdx_it;
        vIdx_it++;
        for (; vIdx_it != path.end(); vIdx_it++) {
            EdgeHandle eh;
            assert(getEdge(mesh, eh, prevIdx, *vIdx_it));
            result.insert(eh);
            prevIdx = *vIdx_it;
        }
    }

    return true;
    
//    // create a new graph with selectedVerices, calculate all paths between them
//    for (size_t i = 0; i < selectedVertices.size(); i++) {
//        for (size_t j = i+1; j < selectedVertices.size(); j++) {
//            int sourceIdx = selectedVertices[i];
//            int destIdx = selectedVertices[j];
//            VertexHandle sourceHandle = mesh->vertex_handle(sourceIdx);
//            VertexHandle destHandle = mesh->vertex_handle(destIdx);
//            
//            // if result is not found, caculate now
//            std::map<std::pair<VertexHandle, VertexHandle>, double>::iterator found =
//                geodesicDistance->find(std::make_pair(sourceHandle, destHandle));
//            if (found == geodesicDistance->end()) {
//                calcGeodesic(mesh, sourceHandle);
//                found = geodesicDistance->find(std::make_pair(sourceHandle, destHandle));
//            }
//            
//            // assign edge & weight
//            double distance = found->second;
//            edges.push_back(std::make_pair(i, j));
//            weights.push_back(distance);
//        }
//    }
//    
//    Graph g(edges.begin(), edges.end(), weights.begin(), selectedVertices.size());
//    
//    std::list<Edge> mst;
//    boost::kruskal_minimum_spanning_tree(g,
//                                         std::back_inserter(mst));
//    auto index = get(boost::vertex_index, g);
//    for (std::list<Edge>::iterator it = mst.begin(); it != mst.end(); it++) {
//        int i = get(index, boost::source(*it, g));
//        int j = get(index, boost::target(*it, g));
//        VertexHandle sourceHandle = mesh->vertex_handle(selectedVertices[i]);
//        VertexHandle destHandle = mesh->vertex_handle(selectedVertices[j]);
//        
//        std::map<std::pair<VertexHandle, VertexHandle>, IdList>::iterator found =
//            geodesicPath->find(std::make_pair(sourceHandle, destHandle));
//        
//        // geodesicPath should exist before calculating spanning tree.
//        assert(found != geodesicPath->end());
//        
//        IdList path = found->second;
//        IdList::iterator vIdx_it = path.begin();
//        int prevIdx = *vIdx_it;
//        vIdx_it++;
//        for (; vIdx_it != path.end(); vIdx_it++) {
//            EdgeHandle eh;
//            assert(getEdge(mesh, eh, prevIdx, *vIdx_it));
//            result.insert(eh);
//            prevIdx = *vIdx_it;
//        }
//    }
}