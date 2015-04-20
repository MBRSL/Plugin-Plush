#include "OpenMesh_Boost_Wrapper.hh"
#include "PlushPatternGenerator.hh"

#include <boost/graph/dijkstra_shortest_paths.hpp>

class WeightFunctor {
private:
    TriMesh *m_mesh;
    const boost::iterator_property_map<std::vector<VertexHandle>::iterator, TriMesh_id_map>
    m_predecessor_pmap;
    VertexHandle *m_currentV;
    double m_maxEdgeLength;
    
    double distanceWeight(TriMesh::Point p1, TriMesh::Point p2) const;
    double textureWeight(HalfedgeHandle heh) const;
    
    double curvatureWeight(VertexHandle v1, VertexHandle v2) const;
    
    double skeletonWeight(VertexHandle v1,
                          VertexHandle v2,
                          TriMesh::Point p1,
                          TriMesh::Point p2) const;
    
    double smoothnessWeight(VertexHandle v1,
                            VertexHandle v2,
                            TriMesh::Point p1,
                            TriMesh::Point p2) const;
    
public:
    WeightFunctor(TriMesh *mesh,
                  VertexHandle *currentV,
                  const boost::iterator_property_map<std::vector<VertexHandle>::iterator, TriMesh_id_map>
                  &predecessor_pmap);
    
    double operator()(HalfedgeHandle heh) const;
};

class Dijkstra_visitor : public boost::default_dijkstra_visitor {
private:
    VertexHandle *currentV;
public:
    Dijkstra_visitor (VertexHandle *currentV) : currentV(currentV) {}
    
    template <class Vertex, class Graph>
    void examine_vertex(Vertex u, Graph &g) {
        *currentV = u;
    }
};