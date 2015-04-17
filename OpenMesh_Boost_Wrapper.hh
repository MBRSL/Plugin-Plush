#ifndef OPENMESH_BOOST_WRAPPER_HH
#define OPENMESH_BOOST_WRAPPER_HH

#include <ObjectTypes/TriangleMesh/TriangleMesh.hh>

#include <boost/config.hpp>
#include <boost/graph/graph_concepts.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>

namespace boost {
    template<>
    struct graph_traits< TriMesh > {
        // Graph
        typedef TriMesh::VertexHandle vertex_descriptor;
//        typedef TriMesh::HalfedgeHandle edge_descriptor;
        typedef TriMesh::EdgeHandle edge_descriptor;
        
        typedef undirected_tag directed_category;
        typedef disallow_parallel_edge_tag edge_parallel_category;
        
        static vertex_descriptor null_vertex() {
            return vertex_descriptor();
        }

        // VertexListGraph
        struct graph_traversal_category :
            public virtual boost::bidirectional_graph_tag,
            public virtual boost::vertex_list_graph_tag,
            public virtual boost::edge_list_graph_tag
        {};
        
        typedef graph_traversal_category traversal_category;
        typedef TriMesh::ConstVertexIter vertex_iterator;
        typedef int vertices_size_type;
        
        // IncidenceGraph
        typedef TriMesh::ConstVertexEdgeIter out_edge_iterator;
        typedef int degree_size_type;
        
        // BidirectionalGraph
        typedef TriMesh::ConstVertexEdgeIter in_edge_iterator;

        // EdgeListGraph
        typedef TriMesh::ConstEdgeIter edge_iterator;
        typedef int edges_size_type;
    };
    
    template<>
    struct vertex_property_type<TriMesh>
    {
        typedef TriMesh::Point type;
    };
}
//===========================================================================
// property maps for TriMesh

// id maps
class TriMesh_id_map
: public boost::put_get_helper<std::size_t, TriMesh_id_map>
{
public:
    typedef boost::readable_property_map_tag category;
    typedef std::size_t value_type;
    typedef std::size_t reference;
    typedef OpenMesh::BaseHandle key_type;
    TriMesh_id_map() { }
    template <class T>
    long operator[](T x) const { return x.idx(); }
};

inline
TriMesh_id_map
get(boost::vertex_index_t, const TriMesh& g) {
    return TriMesh_id_map();
}

inline
TriMesh_id_map
get(boost::edge_index_t, const TriMesh& g) {
    return TriMesh_id_map();
}

// Tags
template <class Tag>
struct TriMesh_property_map { };

template <>
struct TriMesh_property_map<boost::vertex_index_t> {
    struct bind_ {
        typedef TriMesh_id_map type;
        typedef TriMesh_id_map const_type;
    };
};

template <>
struct TriMesh_property_map<boost::edge_index_t> {
    struct bind_ {
        typedef TriMesh_id_map type;
        typedef TriMesh_id_map const_type;
    };
};

// Convert custom property maps to boost::property_map
template<class Tag>
struct boost::property_map<TriMesh, Tag> {
    typedef typename TriMesh_property_map<Tag>::bind_ map_gen;
    typedef typename map_gen::type type;
    typedef typename map_gen::const_type const_type;
};

//
// operations
//
inline boost::graph_traits< TriMesh >::vertices_size_type
num_vertices(const TriMesh& g)
{
    return g.n_vertices();
}

inline boost::graph_traits< TriMesh >::edges_size_type
num_edges(const TriMesh& g)
{
    return g.n_halfedges() / 2;
}

inline boost::graph_traits< TriMesh >::degree_size_type
degree(boost::graph_traits< TriMesh >::vertex_descriptor v
       , const TriMesh& g)
{
    return g.valence(v);
}

inline boost::graph_traits< TriMesh >::degree_size_type
out_degree(boost::graph_traits< TriMesh >::vertex_descriptor v
           , const TriMesh& g)
{
    return g.valence(v);
}

inline boost::graph_traits< TriMesh >::degree_size_type
in_degree(boost::graph_traits< TriMesh >::vertex_descriptor v
          , const TriMesh& g)
{
    return g.valence(v);
}

//inline boost::graph_traits< TriMesh >::vertex_descriptor
//source(
//       boost::graph_traits< TriMesh >::edge_descriptor e,
//       const TriMesh& g)
//{
//    return g.from_vertex_handle(e);
//}
//
//inline boost::graph_traits< TriMesh >::vertex_descriptor
//target(
//       boost::graph_traits< TriMesh >::edge_descriptor e,
//       const TriMesh& g)
//{
//    return g.to_vertex_handle(e);
//}

inline boost::graph_traits< TriMesh >::vertex_descriptor
source(
       boost::graph_traits< TriMesh >::edge_descriptor e,
       const TriMesh& g)
{
    return g.from_vertex_handle(g.halfedge_handle(e,0));
}

inline boost::graph_traits< TriMesh >::vertex_descriptor
target(
       boost::graph_traits< TriMesh >::edge_descriptor e,
       const TriMesh& g)
{
    return g.to_vertex_handle(g.halfedge_handle(e,0));
}

//inline std::pair<boost::graph_traits< TriMesh >::edge_descriptor, bool>
//edge(boost::graph_traits< TriMesh >::vertex_descriptor u
//     , boost::graph_traits< TriMesh >::vertex_descriptor v
//     , const TriMesh& g)
//{
//    // circulate around the inedges of u
//    for (TriMesh::ConstVertexOHalfedgeIter voh_it = g.cvoh_begin(u); voh_it; voh_it++) {
//        if(*voh_it == v) {
//            return std::make_pair(*voh_it, true);
//        }
//    }
//    return std::make_pair(TriMesh::HalfedgeHandle(), false);
//}

//inline std::pair<boost::graph_traits< TriMesh >::edge_descriptor, bool>
//edge(boost::graph_traits< TriMesh >::vertex_descriptor u
//     , boost::graph_traits< TriMesh >::vertex_descriptor v
//     , const TriMesh& g)
//{
//    // circulate around the inedges of u
//    for (TriMesh::ConstVertexOHalfedgeIter cvoh_it = g.cvoh_begin(u); cvoh_it; cvoh_it++) {
//        if(g.to_vertex_handle(*cvoh_it) == v) {
//            return std::make_pair(g.edge_handle(*cvoh_it), true);
//        }
//    }
//    return std::make_pair(TriMesh::EdgeHandle(), false);
//}

inline
std::pair<boost::graph_traits< TriMesh >::vertex_iterator
,boost::graph_traits< TriMesh >::vertex_iterator
>
vertices( const TriMesh& g)
{
    return std::make_pair( g.vertices_begin(), g.vertices_end());
}

inline
std::pair<boost::graph_traits< TriMesh >::edge_iterator
,boost::graph_traits< TriMesh >::edge_iterator
>
edges( const TriMesh& g)
{
    return std::make_pair( g.edges_begin(), g.edges_end() );
}

inline
std::pair<
boost::graph_traits< TriMesh >::in_edge_iterator,
boost::graph_traits< TriMesh >::in_edge_iterator>
in_edges( boost::graph_traits< TriMesh >::vertex_descriptor u
         , const TriMesh& g)
{
    return std::make_pair( g.cve_begin(u), g.cve_end(u) );
}

inline
std::pair<
boost::graph_traits< TriMesh >::out_edge_iterator,
boost::graph_traits< TriMesh >::out_edge_iterator>
out_edges( boost::graph_traits< TriMesh >::vertex_descriptor u
          , const TriMesh& g)
{
    for (auto e_it = g.cve_begin(u); e_it != g.cve_end(u); e_it++) {
        assert(target(*e_it, g) == u);
    }
    return std::make_pair( g.cve_begin(u), g.cve_end(u) );
}

BOOST_CONCEPT_ASSERT((boost::BidirectionalGraphConcept<TriMesh>));
BOOST_CONCEPT_ASSERT((boost::VertexListGraphConcept<TriMesh>));
BOOST_CONCEPT_ASSERT((boost::EdgeListGraphConcept<TriMesh>));
#endif