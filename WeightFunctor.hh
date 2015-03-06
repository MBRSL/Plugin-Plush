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
    double avgEdgeLength;
    
    double distanceWeight(TriMesh::Point p1, TriMesh::Point p2) const;
    double textureWeight(HalfedgeHandle he1, HalfedgeHandle he2) const;
    
    double curvatureWeight(VertexHandle v1, VertexHandle v2) const;
    
    double skeletonWeight(VertexHandle v1,
                          VertexHandle v2,
                          TriMesh::Point p1,
                          TriMesh::Point p2) const;
    
    double smoothnessWeight(boost_vertex_descriptor boost_v1,
                            boost_vertex_descriptor boost_v2,
                            TriMesh::Point p1,
                            TriMesh::Point p2) const;
    
public:
    WeightFunctor(TriMesh *mesh,
                  const Polyhedron &P,
                  boost_vertex_descriptor *currentV,
                  const boost::iterator_property_map<std::vector<boost_vertex_descriptor>::iterator, VertexIdPropertyMap>
                  &predecessor_pmap,
                  const EdgeIdPropertyMap &edgeIndexMap);
    
    double operator()(boost_edge_descriptor e) const;
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