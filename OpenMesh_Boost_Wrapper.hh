#ifndef OPENMESH_BOOST_WRAPPER_HH
#define OPENMESH_BOOST_WRAPPER_HH

#include <ObjectTypes/TriangleMesh/TriangleMesh.hh>

#include <boost/config.hpp>
#include <boost/iterator/iterator_facade.hpp>
#include <boost/graph/graph_concepts.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>

namespace boost {
    template<>
    struct graph_traits< TriMesh > {
        // Graph
        typedef TriMesh::VertexHandle vertex_descriptor;
        typedef TriMesh::HalfedgeHandle edge_descriptor;
        
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
        typedef TriMesh::ConstVertexOHalfedgeIter out_edge_iterator;
        typedef int degree_size_type;
        
        // BidirectionalGraph
        typedef TriMesh::ConstVertexIHalfedgeIter in_edge_iterator;
        
        // EdgeListGraph
        //        typedef TriMesh::ConstHalfedgeIter edge_iterator;
        class edge_iterator
        : public boost::iterator_facade<edge_iterator,
        edge_descriptor,
        bidirectional_traversal_tag,
        const edge_descriptor&>
        {
        public:
            edge_iterator() : m_g(NULL) {}
            
            edge_iterator(const TriMesh* g,
                          bool isEnd)
            : m_g(g) {
                if (!isEnd) {
                    m_eh_iter = m_g->edges_begin();
                    m_current_heh = m_g->halfedge_handle(*m_eh_iter, 0);
                } else {
                    m_eh_iter = m_g->edges_end();
                }
            }
            
        private:
            const edge_descriptor& dereference() const {
                return m_current_heh;
            }
            
            bool equal(const edge_iterator& other) const {
                return *other == m_current_heh;
            }
            
            void increment() {
                m_eh_iter++;
                if (m_eh_iter != m_g->edges_end()) {
                    m_current_heh = m_g->halfedge_handle(*m_eh_iter, 0);
                } else {
                    m_current_heh.invalidate();
                }
            }
            void decrement() {
                m_eh_iter--;
                if (m_eh_iter != m_g->edges_end()) {
                    m_current_heh = m_g->halfedge_handle(*m_eh_iter, 0);
                } else {
                    m_current_heh.invalidate();
                }
            }
            
            const TriMesh *m_g;
            TriMesh::EdgeIter m_eh_iter;
            
            // To prevent from returning local variable
            edge_descriptor m_current_heh;
            
            friend class iterator_core_access;
        };
        
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
    return g.n_edges();
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

inline boost::graph_traits< TriMesh >::vertex_descriptor
source(
       boost::graph_traits< TriMesh >::edge_descriptor e,
       const TriMesh& g)
{
    return g.from_vertex_handle(e);
}

inline boost::graph_traits< TriMesh >::vertex_descriptor
target(
       boost::graph_traits< TriMesh >::edge_descriptor e,
       const TriMesh& g)
{
    return g.to_vertex_handle(e);
}

//inline boost::graph_traits< TriMesh >::vertex_descriptor
//source(
//       boost::graph_traits< TriMesh >::edge_descriptor e,
//       const TriMesh& g)
//{
////    return g.from_vertex_handle(g.halfedge_handle(e,0));
//    return e.source();
//}
//
//inline boost::graph_traits< TriMesh >::vertex_descriptor
//target(
//       boost::graph_traits< TriMesh >::edge_descriptor e,
//       const TriMesh& g)
//{
////    return g.to_vertex_handle(g.halfedge_handle(e,0));
//    return e.target();
//}

//inline std::pair<boost::graph_traits< TriMesh >::edge_descriptor, bool>
//edge(boost::graph_traits< TriMesh >::vertex_descriptor u
//     , boost::graph_traits< TriMesh >::vertex_descriptor v
//     , const TriMesh& g)
//{
//    typedef boost::graph_traits<TriMesh>::edge_descriptor edge_descriptor;
//    // circulate around the inedges of u
//    for (TriMesh::HalfedgeHandle voh : g.voh_range(u)) {
//        if(g.to_vertex_handle(voh) == v) {
//            edge_descriptor e(&g);
//            e.m_heh = voh;
//            e.m_eh = g.edge_handle(voh);
//            return std::make_pair(e, true);
//        }
//    }
//    return std::make_pair(edge_descriptor(NULL), false);
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
    typedef typename boost::graph_traits< TriMesh >::edge_iterator Iter;
    Iter ve_it(&g, false);
    Iter ve_ite(&g, true);
    
    return std::make_pair( ve_it, ve_ite );
}

inline
std::pair<
boost::graph_traits< TriMesh >::in_edge_iterator,
boost::graph_traits< TriMesh >::in_edge_iterator>
in_edges( boost::graph_traits< TriMesh >::vertex_descriptor u
         , const TriMesh& g)
{
    return std::make_pair( g.cvih_begin(u), g.cvih_end(u) );
}

inline
std::pair<
boost::graph_traits< TriMesh >::out_edge_iterator,
boost::graph_traits< TriMesh >::out_edge_iterator>
out_edges( boost::graph_traits< TriMesh >::vertex_descriptor u
          , const TriMesh& g)
{
    // Boost require the ouput of out_edges(u, g) to be (u, v) not (v, u)
    for (auto e_it = g.cvoh_begin(u); e_it != g.cvoh_end(u); e_it++) {
        assert(source(*e_it, g) == u);
    }
    
    return std::make_pair( g.cvoh_begin(u), g.cvoh_end(u) );
}

BOOST_CONCEPT_ASSERT((boost::BidirectionalGraphConcept<TriMesh>));
BOOST_CONCEPT_ASSERT((boost::VertexListGraphConcept<TriMesh>));
BOOST_CONCEPT_ASSERT((boost::EdgeListGraphConcept<TriMesh>));
#endif