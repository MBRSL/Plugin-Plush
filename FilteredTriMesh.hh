#ifndef TRIMESH_FILTER_HH
#define TRIMESH_FILTER_HH

#include "common.hh"

#include <ObjectTypes/TriangleMesh/TriangleMesh.hh>
#include <boost/range/adaptor/filtered.hpp>

// True means it's contained in this sub-mesh
struct Predicate {
    std::vector<bool> m_vertex_filter;
    std::vector<bool> m_edge_filter;
    std::vector<bool> m_halfedge_filter;
    std::vector<bool> m_face_filter;
    
    bool operator()(VertexHandle vh) const {
        return m_vertex_filter[vh.idx()];
    }
    bool operator()(EdgeHandle eh) const {
        return m_edge_filter[eh.idx()];
    }
    bool operator()(HalfedgeHandle heh) const {
        return m_halfedge_filter[heh.idx()];
    }
    bool operator()(FaceHandle fh) const {
        return m_face_filter[fh.idx()];
    }
};

struct Boundary_predicate {
    std::vector<bool> m_boundary_vertex_filter;
    std::vector<bool> m_boundary_edge_filter;
    bool operator()(VertexHandle vh) const {
        return m_boundary_vertex_filter[vh.idx()];
    }
    bool operator()(EdgeHandle eh) const {
        return m_boundary_edge_filter[eh.idx()];
    }
};

template <class Handle, class Predicate>
struct Range {
    typedef boost::filter_iterator<Predicate, Handle> FilterIter;
    FilterIter m_begin;
    FilterIter m_end;
    
    Range(FilterIter begin, FilterIter end) : m_begin(begin), m_end(end) {}
    FilterIter begin() {
        return m_begin;
    }
    FilterIter end() {
        return m_end;
    }
};

/**
 * @brief This class stores essential information of a submesh of trimesh.
 */
class FilteredTriMesh
{
private:
    TriMesh *m_mesh;
    Predicate predicate;
    Boundary_predicate boundary_predicate;
    
    /* @brief Given a halfedge handle, return the flattened position of "from" vertex of this halfedge.
     * We use halfedge handle to store the flattened position of boundary vertex because a boundary vertex in the original mesh
     * may map to two vertices after split. In such case, vertex-to-vertex mapping is not one-to-one.
     * We use halfedge-to-vertex mapping which is one-to-one.
     */
    std::map<HalfedgeHandle, int> boundary_vertex_mapping;
    std::map<VertexHandle, int> non_boundary_vertex_mapping;
    
    /// Stores the flattened position with the indices of mapping above
    std::vector<TriMesh::Point> flattened_position;

    HalfedgeHandle get_nearby_boundary_halfedge(HalfedgeHandle heh);
public:
    std::vector<bool> subMesh_id;
    
    /// You shouldn't use this constructor unless you just need to allocate a variable and assigned it later.
    FilteredTriMesh() {}
    FilteredTriMesh(TriMesh *mesh,
                    std::vector<VertexHandle> vertices,
                    std::vector<EdgeHandle> edges,
                    std::vector<FaceHandle> faces,
                    std::vector<EdgeHandle> boundary_edges);
    
    bool is_valid(VertexHandle vh) const;
    bool is_valid(EdgeHandle eh) const;
    bool is_valid(HalfedgeHandle heh) const;
    bool is_valid(FaceHandle fh) const;
    bool is_boundary(VertexHandle vh) const;
    bool is_boundary(EdgeHandle eh) const;
    bool is_boundary(HalfedgeHandle heh) const;

    /* For the T-shape boundary like darts, we called the extrusion part dual boundary.
     * Each of them will be split into two edges after flattening.
     * However, the original mesh contains only one edge. So we need to take care of such edges
     */
    ///@{
    bool is_dual_boundary(EdgeHandle eh) const;
    bool is_dual_boundary(HalfedgeHandle eh) const;
    ///@}

    /// The filtered iterators of original operations
    ///@{
    unsigned long n_vertices() const;
    
    Range<VertexIter, Predicate> vertices() const;
    Range<EdgeIter, Predicate> edges() const;
    Range<HalfedgeIter, Predicate> halfedges() const;
    Range<FaceIter, Predicate> faces() const;
    Range<TriMesh::ConstFaceVertexIter, Predicate> fv_range(FaceHandle fh) const;
    Range<TriMesh::ConstFaceHalfedgeIter, Predicate> fh_range(FaceHandle fh) const;
    
    Range<VertexIter, Boundary_predicate> boundary_vertices() const;
    Range<EdgeIter, Boundary_predicate> boundary_edges() const;
    ///@}

    /// Filtered functions
    ///@{
    VertexHandle vertex_handle(int idx);
    
    EdgeHandle edge_handle(HalfedgeHandle heh) const;
    
    HalfedgeHandle halfedge_handle(EdgeHandle eh, int idx) const;
    HalfedgeHandle halfedge_handle(FaceHandle fh) const;
    HalfedgeHandle opposite_halfedge_handle(HalfedgeHandle heh) const;
    HalfedgeHandle next_halfedge_handle(HalfedgeHandle heh) const;
    HalfedgeHandle next_boundary_halfedge_handle(HalfedgeHandle heh) const;
    HalfedgeHandle prev_halfedge_handle(HalfedgeHandle heh) const;
    HalfedgeHandle prev_boundary_halfedge_handle(HalfedgeHandle heh) const;
    
    VertexHandle from_vertex_handle(HalfedgeHandle heh) const;
    VertexHandle to_vertex_handle(HalfedgeHandle heh) const;
    
    size_t valence(VertexHandle vh) const;
    TriMesh::Point& point(VertexHandle vh);
    const TriMesh::Color& color(FaceHandle fh) const;
    double calc_edge_length(EdgeHandle eh) const;
    double calc_edge_length(HalfedgeHandle heh) const;
    double calc_sector_angle(HalfedgeHandle heh) const;
    double calc_sector_area(HalfedgeHandle heh) const;
    ///@}
    
    /// Additional operations
    ///@{
    int n_duplicated_vertices();

    TriMesh::Point& flattened_non_boundary_point(VertexHandle vh);
    TriMesh::Point& flattened_boundary_point(HalfedgeHandle heh);

    double getSumInnerAngle(HalfedgeHandle from, HalfedgeHandle to);
    
    int get_non_boundary_vertex_idx(VertexHandle vh);
    int get_boundary_vertex_idx(HalfedgeHandle heh);
    ///@}
    
    /// Use this if you want to use operations for original unfiltered mesh
    TriMesh *mesh() const;
};
#endif //TRIMESH_FILTER_HH
