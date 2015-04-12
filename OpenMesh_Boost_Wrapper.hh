#ifndef OPENMESH_BOOST_WRAPPER_HH
#define OPENMESH_BOOST_WRAPPER_HH

#include <ObjectTypes/TriangleMesh/TriangleMesh.hh>

#include <boost/config.hpp>
#include <boost/iterator/iterator_facade.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>

namespace boost {
    template<>
    struct graph_traits< TriMesh > {
        typedef typename TriMesh::Point vertex_property_type;
        typedef TriMesh::VertexHandle vertex_descriptor;
        typedef TriMesh::HalfedgeHandle edge_descriptor;
        typedef TriMesh::FaceHandle face_descriptor;
        // iterator typedefs...
        
        struct graph_traversal_category :   public virtual boost::bidirectional_graph_tag,
        public virtual boost::vertex_list_graph_tag,
        public virtual boost::edge_list_graph_tag
        {};
        typedef bidirectional_tag directed_category;
        typedef disallow_parallel_edge_tag edge_parallel_category;
        typedef graph_traversal_category traversal_category;
        
        typedef int degree_size_type;
        typedef int vertices_size_type;
        typedef int edges_size_type;
        
        typedef TriMesh::ConstVertexOHalfedgeIter out_edge_iterator;
        typedef TriMesh::ConstVertexIHalfedgeIter in_edge_iterator;
        typedef TriMesh::ConstHalfedgeIter edge_iterator;
        typedef TriMesh::ConstVertexIter vertex_iterator;
    };

    template<class Graph>
    typename boost::graph_traits< Graph >::vertices_size_type
    num_vertices(const Graph& g)
    {
        return g.n_vertices();
    }
    
    template<class Graph>
    typename boost::graph_traits< Graph >::edges_size_type
    num_edges(const Graph& g)
    {
        return g.n_halfedges() / 2;
    }
    
    template<class Graph>
    typename boost::graph_traits< Graph >::degree_size_type
    degree(typename boost::graph_traits< Graph >::vertex_descriptor v
           , const Graph& g)
    {
        return g.valence(v);
    }
    
    template<class Graph>
    typename boost::graph_traits< Graph >::degree_size_type
    out_degree(typename boost::graph_traits< Graph >::vertex_descriptor v
               , const Graph& g)
    {
        return g.valence(v);
    }
    
    template<class Graph>
    typename boost::graph_traits< Graph >::degree_size_type
    in_degree(typename boost::graph_traits< Graph >::vertex_descriptor v
              , const Graph& g)
    {
        return g.valence(v);
    }
    
    template<class Graph>
    typename boost::graph_traits< Graph >::vertex_descriptor
    source(
           typename boost::graph_traits< Graph >::edge_descriptor e,
           const Graph& g)
    {
        return g.from_vertex_handle(e);
    }
    
    template<class Graph>
    typename boost::graph_traits< Graph >::vertex_descriptor
    target(
           typename boost::graph_traits< Graph >::edge_descriptor e,
           const Graph& g)
    {
        return g.to_vertex_handle(e);
    }
    
    template<class Graph>
    std::pair<
    typename boost::graph_traits< Graph >::edge_descriptor
    , bool>
    edge(typename boost::graph_traits< Graph >::vertex_descriptor u
         , typename boost::graph_traits< Graph >::vertex_descriptor v
         , const Graph& g)
    {
        typedef boost::graph_traits< Graph > Traits;
        typedef typename Traits::edge_descriptor edge;
        
        // circulate around the inedges of u
        for (typename Graph::ConstVertexIHalfedgeIter vih_it = g.cvih_begin(u); vih_it; vih_it++) {
            if(*vih_it == v) {
                return std::make_pair(edge(g.opposite_halfedge_handle(*vih_it)), true);
            }
        }
        return std::make_pair(edge(), false);
    }
    
    template<class Graph>
    inline std::pair<typename boost::graph_traits< Graph >::vertex_iterator
    ,typename boost::graph_traits< Graph >::vertex_iterator
    >
    vertices( const Graph& g)
    {
        return std::make_pair( g.vertices_begin(), g.vertices_end());
    }
    
    template<class Graph>
    inline std::pair<typename boost::graph_traits< Graph >::edge_iterator
    ,typename boost::graph_traits< Graph >::edge_iterator
    >
    edges( const Graph& g)
    {
        return std::make_pair( g.halfedges_begin(), g.halfedges_end() );
    }
    
    template<class Graph>
    inline std::pair<typename boost::graph_traits< Graph >::in_edge_iterator
    ,typename boost::graph_traits< Graph >::in_edge_iterator
    >
    in_edges( typename boost::graph_traits< Graph >::vertex_descriptor u
             , const Graph& g)
    {
        return std::make_pair( g.cvih_begin(u), g.cvih_end(u) );
    }
    
    template<class Graph>
    inline std::pair<typename boost::graph_traits< Graph >::out_edge_iterator
    ,typename boost::graph_traits< Graph >::out_edge_iterator
    >
    out_edges( typename boost::graph_traits< Graph >::vertex_descriptor u
              , const Graph& g)
    {
        return std::make_pair( g.cvoh_begin(u), g.cvoh_end(u) );
    }
    
    //===========================================================================
    // property maps for TriMesh
    
    class TriMesh_Vertices_id_map
    : public put_get_helper<int, TriMesh_Vertices_id_map>
    {
    public:
        typedef readable_property_map_tag category;
        typedef int value_type;
        typedef int reference;
        typedef TriMesh::VertexHandle key_type;
        TriMesh_Vertices_id_map() { }
        template <class T>
        long operator[](T x) const { return x.idx(); }
    };

    class TriMesh_Halfedges_id_map
    : public put_get_helper<int, TriMesh_Halfedges_id_map>
    {
    public:
        typedef readable_property_map_tag category;
        typedef int value_type;
        typedef int reference;
        typedef TriMesh::HalfedgeHandle key_type;
        TriMesh_Halfedges_id_map() { }
        template <class T>
        long operator[](T x) const { return x.idx(); }
    };
    
    inline TriMesh_Vertices_id_map
    get(vertex_index_t, const TriMesh& g) {
        return TriMesh_Vertices_id_map();
    }

    inline TriMesh_Halfedges_id_map
    get(edge_index_t, const TriMesh& g) {
        return TriMesh_Halfedges_id_map();
    }
} // namespace boost

#endif