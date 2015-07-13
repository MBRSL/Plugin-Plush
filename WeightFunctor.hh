#include "OpenMesh_Boost_Wrapper.hh"
#include "PlushPatternGenerator.hh"

#include <boost/graph/dijkstra_shortest_paths.hpp>

class WeightFunctor {
private:
    TriMesh *m_mesh;
    VertexHandle &m_currentV;
    VertexHandle m_prevV;
    const boost::iterator_property_map<std::vector<VertexHandle>::iterator, TriMesh_id_map>
        *m_predecessor_pmap;
    double m_maxEdgeLength;
    
    double m_distanceCoefficient;
    double m_textureCoefficient;
    double m_curvatureCoefficient;
    double m_skeletonCoefficient;
    double m_pathCoefficient;

    double distanceWeight(EdgeHandle eh) const;
    double textureWeight(HalfedgeHandle heh) const;
    
//    double curvatureWeight(VertexHandle v1, VertexHandle v2) const;
    double curvatureWeight(EdgeHandle eh) const;
    
    double skeletonWeight(EdgeHandle eh,
                          VertexHandle v1,
                          VertexHandle v2,
                          TriMesh::Point p1,
                          TriMesh::Point p2) const;
    
    double smoothnessWeight(VertexHandle v1,
                            VertexHandle v2,
                            TriMesh::Point p1,
                            TriMesh::Point p2) const;
    
public:
    WeightFunctor(TriMesh *mesh,
                  VertexHandle &currentV,
                  const boost::iterator_property_map<std::vector<VertexHandle>::iterator, TriMesh_id_map>
                  *predecessor_pmap,
                  double distanceCoefficient,
                  double textureCoefficient,
                  double curvatureCoefficient,
                  double skeletonCoefficient,
                  double pathCoefficient);
    
    double operator()(HalfedgeHandle heh) const;
    double operator()(std::vector<HalfedgeHandle> segment);
};

class Dijkstra_visitor : public boost::default_dijkstra_visitor {
private:
    VertexHandle &currentV;
public:
    Dijkstra_visitor (VertexHandle &currentV) : currentV(currentV) {}
    
    template <class Vertex, class Graph>
    void examine_vertex(Vertex u, Graph &g) {
        currentV = u;
    }
};
