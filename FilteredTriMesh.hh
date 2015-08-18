#ifndef TRIMESH_FILTER_HH
#define TRIMESH_FILTER_HH

#include "common.hh"

#include <ObjectTypes/TriangleMesh/TriangleMesh.hh>

#include <boost/iterator/iterator_facade.hpp>
#include <boost/dynamic_bitset.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/map.hpp>

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

class FilteredVertexIter;
class FilteredEdgeIter;
class FilteredHalfedgeIter;
class FilteredFaceIter;
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
    
    friend class boost::serialization::access;
    template <class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & predicate.m_vertex_filter;
        ar & predicate.m_edge_filter;
        ar & predicate.m_halfedge_filter;
        ar & predicate.m_face_filter;
        ar & boundary_predicate.m_boundary_vertex_filter;
        ar & boundary_predicate.m_boundary_edge_filter;
        ar & merged_subMesh_idx;
        ar & merged_seam_idx;
        ar & merged_edge_idx;
        ar & boundary_vertex_mapping;
        ar & non_boundary_vertex_mapping;
        ar & flattened_position;
        ar & max_distortion;
        ar & seam_score;
        ar & n_merged_seams;
    }
public:
    /* Use this to record which sub-meshes do this patch have merged
     * A patch with bits [1,0,0,0,0,....] mean it contains #1 sub-mesh
     * A patch with bits [0,1,1,0,0,....] mean it merged #2 and #3 sub-meshes
     */
    typedef boost::dynamic_bitset<> Patch_idx;
    /* Use this to record which seam segment do this patch have merged
     * A patch with bits [1,0,0,0,0,....] mean it merged #1 seam
     * A patch with bits [0,1,1,0,0,....] mean it merged #2 and #3 seams
     */
    typedef boost::dynamic_bitset<> Seam_idx;

    Patch_idx merged_subMesh_idx;
    Seam_idx merged_seam_idx;
    /* Use this to record which edge handle do this patch have merged
     * A merged edge should not be added to boundary anymore
     * A patch with bits [1,0,0,0,0,....] mean it merged #1 edge
     * A patch with bits [0,1,1,0,0,....] mean it merged #2 and #3 edges
     */
    boost::dynamic_bitset<> merged_edge_idx;
    
    // The face area difference between original sub-mesh & flattened sub-mesh
    double max_distortion = -1;

    // The seam score by evaluating weightFunction (lower is better)
    double seam_score = -1;
    
    int n_merged_seams = 0;

    /// You shouldn't use this constructor unless you just need to allocate a variable and assigned it later.
    FilteredTriMesh() {}
    FilteredTriMesh(TriMesh *mesh,
                    boost::archive::text_iarchive &ia);
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
    
    Range<FilteredVertexIter> vertices() const;
    Range<FilteredEdgeIter> edges() const;
    Range<FilteredHalfedgeIter> halfedges() const;
    Range<FilteredFaceIter> faces() const;
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
    
    Range<FilteredVertexIter> boundary_vertices() const;
    Range<FilteredEdgeIter> boundary_edges() const;
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

template <class HandleType, class Derived>
class FilteredIterator :
public boost::iterator_facade<  FilteredIterator<HandleType, Derived>,
                                HandleType,
                                boost::forward_traversal_tag,
                                HandleType>
{
public:
    FilteredIterator(const FilteredTriMesh &mesh,
                     const boost::dynamic_bitset<> &filter,
                     size_t pos = 0) :
    m_mesh(mesh), m_filter(filter), m_pos(pos) {
    }
    
protected:
    friend class boost::iterator_core_access;
    
    void increment() {
        m_pos = m_filter.find_next(m_pos);
    }
    
    bool equal(FilteredIterator<HandleType, Derived> const& other) const {
        auto handle1 = this->dereference();
        auto handle2 = *other;
        return handle1 == handle2;
    }
    
    HandleType dereference() const {
        auto handle = static_cast<const Derived*>(this)->dereference_impl();
        return handle;
    }
    
    const FilteredTriMesh &m_mesh;
    const boost::dynamic_bitset<> &m_filter;
    size_t m_pos;
};

class FilteredVertexIter : public FilteredIterator<VertexHandle, FilteredVertexIter> {
public:
    FilteredVertexIter(const FilteredTriMesh &mesh,
                       const boost::dynamic_bitset<> &filter,
                       size_t pos) :
    FilteredIterator(mesh, filter, pos) {
    }
    VertexHandle dereference_impl() const {
        return m_mesh.mesh()->vertex_handle(m_pos);
    }
};
class FilteredEdgeIter : public FilteredIterator<EdgeHandle, FilteredEdgeIter> {
public:
    FilteredEdgeIter(const FilteredTriMesh &mesh,
                         const boost::dynamic_bitset<> &filter,
                         size_t pos) :
    FilteredIterator(mesh, filter, pos) {
    }
    EdgeHandle dereference_impl() const {
        return m_mesh.mesh()->edge_handle(m_pos);
    }
};
class FilteredHalfedgeIter : public FilteredIterator<HalfedgeHandle, FilteredHalfedgeIter> {
public:
    FilteredHalfedgeIter(const FilteredTriMesh &mesh,
                             const boost::dynamic_bitset<> &filter,
                             size_t pos) :
    FilteredIterator(mesh, filter, pos) {
    }
    HalfedgeHandle dereference_impl() const {
        return m_mesh.mesh()->halfedge_handle(m_pos);
    }
};
class FilteredFaceIter : public FilteredIterator<FaceHandle, FilteredFaceIter> {
public:
    FilteredFaceIter(const FilteredTriMesh &mesh,
                         const boost::dynamic_bitset<> &filter,
                         size_t pos) :
    FilteredIterator(mesh, filter, pos) {
    }
    FaceHandle dereference_impl() const {
        return m_mesh.mesh()->face_handle(m_pos);
    }
};

namespace boost {
    namespace serialization {
        template <typename Ar, typename Block, typename Alloc>
        void save(Ar& ar, dynamic_bitset<Block, Alloc> const& bs, unsigned) {
            size_t num_bits = bs.size();
            std::vector<Block> blocks(bs.num_blocks());
            to_block_range(bs, blocks.begin());
            
            ar & num_bits & blocks;
        }
        
        template <typename Ar, typename Block, typename Alloc>
        void load(Ar& ar, dynamic_bitset<Block, Alloc>& bs, unsigned) {
            size_t num_bits;
            std::vector<Block> blocks;
            ar & num_bits & blocks;
            
            bs.resize(num_bits);
            from_block_range(blocks.begin(), blocks.end(), bs);
            bs.resize(num_bits);
        }
        
        template <typename Ar, typename Block, typename Alloc>
        void serialize(Ar& ar, dynamic_bitset<Block, Alloc>& bs, unsigned version) {
            split_free(ar, bs, version);
        }
        
        template <typename Ar>
        void serialize(Ar& ar, VertexHandle& handle, unsigned version) {
            int idx = handle.idx();
            ar & idx;
            handle = VertexHandle(idx);
        }

        template <typename Ar>
        void serialize(Ar& ar, EdgeHandle& handle, unsigned version) {
            int idx = handle.idx();
            ar & idx;
            handle = EdgeHandle(idx);
        }

        template <typename Ar>
        void serialize(Ar& ar, HalfedgeHandle& handle, unsigned version) {
            int idx = handle.idx();
            ar & idx;
            handle = HalfedgeHandle(idx);
        }
        
        template <typename Ar>
        void serialize(Ar& ar, FaceHandle& handle, unsigned version) {
            int idx = handle.idx();
            ar & idx;
            handle = FaceHandle(idx);
        }

        template <typename Ar>
        void serialize(Ar& ar, OpenMesh::Vec3d& vec, unsigned version) {
            ar & vec[0];
            ar & vec[1];
            ar & vec[2];
        }
    }
}
#endif //TRIMESH_FILTER_HH
