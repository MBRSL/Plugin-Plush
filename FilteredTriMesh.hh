#ifndef TRIMESH_FILTER_HH
#define TRIMESH_FILTER_HH

#include "common.hh"

#include <ObjectTypes/TriangleMesh/TriangleMesh.hh>

#include <boost/range/adaptor/filtered.hpp>
#include <boost/dynamic_bitset.hpp>

// True means it's contained in this sub-mesh
struct Predicate {
    boost::dynamic_bitset<> m_vertex_filter;
    boost::dynamic_bitset<> m_edge_filter;
    boost::dynamic_bitset<> m_halfedge_filter;
    boost::dynamic_bitset<> m_face_filter;
    
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
    boost::dynamic_bitset<> m_boundary_vertex_filter;
    boost::dynamic_bitset<> m_boundary_edge_filter;
    bool operator()(VertexHandle vh) const {
        return m_boundary_vertex_filter[vh.idx()];
    }
    bool operator()(EdgeHandle eh) const {
        return m_boundary_edge_filter[eh.idx()];
    }
};

template <class Iter, class Pred>
struct FilteredRange {
    typedef boost::filter_iterator<const Pred&, Iter> FilterIter;
    const Pred &m_pred;
    FilterIter m_begin;
    FilterIter m_end;
    
    FilteredRange(const Pred &pred, Iter begin, Iter end) :
        m_pred(pred),
        m_begin(boost::make_filter_iterator<const Pred&>(pred, begin, end)),
        m_end(boost::make_filter_iterator<const Pred&>(pred, end, end))
    {
    }
    FilteredRange(FilteredRange&& other) noexcept :
    m_pred(other.m_pred),
    m_begin(std::move(other.m_begin)),
    m_end(std::move(other.m_end))
    {
    }
    FilteredRange& operator=(FilteredRange&& other) noexcept
    {
        if (this != &other) {
            m_pred = other.m_pred;
            m_begin = std::move(other.m_begin);
            m_end = std::move(other.m_end);
        }
        return *this;
    }

    
    FilterIter& begin() {
        return m_begin;
    }
    FilterIter& end() {
        return m_end;
    }
};

template <class Iter>
struct Range {
    Iter m_begin;
    Iter m_end;

    Range(Iter begin, Iter end) : m_begin(begin), m_end(end) {}
    Iter begin() {
        return m_begin;
    }
    Iter end() {
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
    /* Use this to record which sub-meshes do this patch have merged
     * A patch with bits [1,0,0,0,0,....] mean it contains #1 sub-mesh
     * A patch with bits [0,1,1,0,0,....] mean it merged #2 and #3 sub-meshes
     */
    boost::dynamic_bitset<> merged_subMesh_idx;
    /* Use this to record which seam segment do this patch have merged
     * A patch with bits [1,0,0,0,0,....] mean it merged #1 seam
     * A patch with bits [0,1,1,0,0,....] mean it merged #2 and #3 seams
     */
    boost::dynamic_bitset<> merged_seam_idx;
    /* Use this to record which edge handle do this patch have merged
     * A merged edge should not be added to boundary anymore
     * A patch with bits [1,0,0,0,0,....] mean it merged #1 edge
     * A patch with bits [0,1,1,0,0,....] mean it merged #2 and #3 edges
     */
    boost::dynamic_bitset<> merged_edge_idx;
    
    // The face area difference between original sub-mesh & flattened sub-mesh
    double max_distortion = -1;
    
    /// You shouldn't use this constructor unless you just need to allocate a variable and assigned it later.
    FilteredTriMesh() {}
    FilteredTriMesh(TriMesh *mesh,
                    std::set<FaceHandle> &faces,
                    std::set<EdgeHandle> &boundary_edges);

    /// Copy constructor (backward compatible with old boost code)
    FilteredTriMesh(const FilteredTriMesh& other);
    FilteredTriMesh& operator=(const FilteredTriMesh& other);
    /// Move constructor
    FilteredTriMesh(FilteredTriMesh&& other) noexcept;
    FilteredTriMesh& operator=(FilteredTriMesh&& other) noexcept;
    
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
    
    FilteredRange<VertexIter, Predicate> vertices() const;
    FilteredRange<EdgeIter, Predicate> edges() const;
    FilteredRange<HalfedgeIter, Predicate> halfedges() const;
    FilteredRange<FaceIter, Predicate> faces() const;
    /// There seems to be some problem about filter_iterator with circulator
    /// You should use *(iter.base()) instead of *iter
    ///@{
//    FilteredRange<TriMesh::ConstFaceVertexIter, Predicate> fv_range(FaceHandle fh) const;
//    FilteredRange<TriMesh::ConstFaceHalfedgeIter, Predicate> fh_range(FaceHandle fh) const;
//    boost::filter_iterator<Predicate, TriMesh::ConstFaceVertexIter> cfv_begin(FaceHandle fh) const;
//    boost::filter_iterator<Predicate, TriMesh::ConstFaceVertexIter> cfv_end(FaceHandle fh) const;
    
    /// For face-halfedge circulator, if a face is valid, then its all related halfedges are valid
    /// Thus we don't need to filter them
    Range<TriMesh::ConstFaceHalfedgeIter> fh_range(FaceHandle fh) const;
    ///@}
    
    FilteredRange<VertexIter, Boundary_predicate> boundary_vertices() const;
    FilteredRange<EdgeIter, Boundary_predicate> boundary_edges() const;
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
