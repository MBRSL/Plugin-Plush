//
//  Geodesic.cc
//  OpenFlipper
//
//  Created by 饒亢 on 2015/1/13.
//
//

#include "PlushPlugin.hh"
#include "CGAL_Polyhedron_builder.hh"

#include <CGAL/boost/graph/dijkstra_shortest_paths.h>
#include <boost/graph/kruskal_min_spanning_tree.hpp>

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
    
    // assign edge weight
    std::vector<double> weight(boost::num_edges(P));
    
    index = 0;
    boost_edge_iterator eit, eit_end;
    for (boost::tie(eit, eit_end) = boost::edges(P); eit != eit_end; ++eit) {
        boost_halfedge_descriptor hed = halfedge(*eit, P);
        hed->id() = index;
        
        weight[index++] = getEdgeWeight(mesh, hed->vertex()->id(), hed->opposite()->vertex()->id());
    }
    
    // prepare property maps for dijkstra algorithm
    typedef boost::property_map<Polyhedron, boost::vertex_external_index_t>::type VertexIdPropertyMap;
    VertexIdPropertyMap vertex_index_pmap = get(CGAL::vertex_external_index, P);
    
    typedef boost::property_map<Polyhedron, boost::edge_external_index_t>::type EdgeIdPropertyMap;
    EdgeIdPropertyMap edge_index_pmap = get(boost::edge_external_index, P);
    
    std::vector<boost_vertex_descriptor> predecessor(boost::num_vertices(P));
    boost::iterator_property_map<std::vector<boost_vertex_descriptor>::iterator, VertexIdPropertyMap> predecessor_pmap(predecessor.begin(), vertex_index_pmap);
    boost::iterator_property_map<std::vector<double>::iterator, VertexIdPropertyMap> distance_pmap(distance.begin(), vertex_index_pmap);
    boost::iterator_property_map<std::vector<double>::iterator, EdgeIdPropertyMap> weight_pmap(weight.begin(), edge_index_pmap);
    
    boost::dijkstra_shortest_paths(P, source,
                                   vertex_index_map(vertex_index_pmap)
                                   .weight_map(weight_pmap)
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
void PlushPlugin::calcSpanningTree(TriMesh *mesh, int meshId, std::set<EdgeHandle> &result, IdList selectedVertices) {
    calcSpanningTree(mesh, meshId, result, selectedVertices, (selectedVertices.size()-1+selectedVertices.size()-2)*(selectedVertices.size()-1)/2);
}

void PlushPlugin::calcSpanningTree(TriMesh *mesh, int meshId, std::set<EdgeHandle> &result, IdList selectedVertices, int edges) {
    typedef boost::adjacency_list < boost::vecS, boost::vecS, boost::undirectedS,
        boost::no_property, boost::property < boost::edge_weight_t, int > > Graph;
    typedef boost::graph_traits < Graph >::edge_descriptor Edge;
    typedef boost::graph_traits < Graph >::vertex_descriptor Vertex;
//    typedef std::pair<int, int> E;
    
//    std::vector<E> edges;
//    std::vector<double> weights;
    std::vector<std::pair<IdList, double> > distance;
    
    // create a new graph with selectedVerices, calculate all paths between them
    for (size_t i = 0; i < selectedVertices.size(); i++) {
        for (size_t j = i+1; j < selectedVertices.size(); j++) {
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
            distance.push_back(std::make_pair(path, cost));
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
double PlushPlugin::getEdgeWeight(TriMesh *mesh, int v1No, int v2No) {
    VertexHandle v1 = mesh->vertex_handle(v1No);
    VertexHandle v2 = mesh->vertex_handle(v2No);
    return getEdgeWeight(mesh, v1, v2);
}

double PlushPlugin::getEdgeWeight(TriMesh *mesh, VertexHandle v1, VertexHandle v2) {
    EdgeHandle eh;
    getEdge(mesh, eh, v1, v2);
    return getEdgeWeight(mesh, eh);
}

double PlushPlugin::getEdgeWeight(TriMesh *mesh, EdgeHandle eh) {
    HalfedgeHandle he1 = mesh->halfedge_handle(eh, 0);
    HalfedgeHandle he2 = mesh->halfedge_handle(eh, 1);
    
    // if the two faces along this edge are different color, set weight of this edge to almost 0
    // we encourage path go through the boundary of different colors
//    if (!mesh->is_boundary(he1) && !mesh->is_boundary(he2)
//        &&  mesh->color(mesh->face_handle(he1)) != mesh->color(mesh->face_handle(he2))) {
//        return 0;
//    } else {
        const TriMesh::Point p1 = mesh->point(mesh->from_vertex_handle(he1));
        const TriMesh::Point p2 = mesh->point(mesh->to_vertex_handle(he1));
        return std::sqrt((p1-p2).norm());
        //            double curvature1 = mesh->property(maxCurvatureHandle, v1);
        //            double curvature2 = mesh->property(maxCurvatureHandle, v2);
        //            // clamp curvature to [-1, 1]
        //            curvature1 = fmax(fmin(curvature1, 1), -1);
        //            curvature2 = fmax(fmin(curvature2, 1), -1);
        //
        //            // weight using curvature
        //            // encourage path with +1/-1 curvature, not 0
        //            return 1 - abs(curvature1 + curvature2)/2;
//    }
}