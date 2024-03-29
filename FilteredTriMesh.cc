#include "FilteredTriMesh.hh"

/**
 Given a halfedge pointing "from" a boundary vertex, find its nearby boundary halfedge.
 @param <#parameter#>
 @return Return the nearby boundary halfedge.
 @retval <#meaning#>
 */
HalfedgeHandle FilteredTriMesh::get_nearby_boundary_halfedge(HalfedgeHandle heh) {
    VertexHandle vh = m_mesh->from_vertex_handle(heh);
    assert(is_boundary(vh));
    
    HalfedgeHandle current_heh = heh;
    do {
        current_heh = m_mesh->next_halfedge_handle(m_mesh->opposite_halfedge_handle(current_heh));
    } while (current_heh != heh && (!is_valid(current_heh) || !is_boundary(current_heh)));
    assert(current_heh != heh && "No boundary halfedge is found?");
    
    return current_heh;
}

FilteredTriMesh::FilteredTriMesh(TriMesh *mesh, boost::archive::text_iarchive &ia) :
m_mesh(mesh) {
    ia >> *this;
}

FilteredTriMesh::FilteredTriMesh(TriMesh *mesh,
                                 std::set<FaceHandle> &faces,
                                 std::set<EdgeHandle> &boundary_edges) : m_mesh(mesh) {
    predicate.m_vertex_filter.resize(m_mesh->n_vertices(), false);
    predicate.m_edge_filter.resize(m_mesh->n_edges(), false);
    predicate.m_halfedge_filter.resize(m_mesh->n_halfedges(), false);
    predicate.m_face_filter.resize(m_mesh->n_faces(), false);

    boundary_predicate.m_boundary_vertex_filter.resize(m_mesh->n_vertices(), false);
    boundary_predicate.m_boundary_edge_filter.resize(m_mesh->n_edges(), false);
    
    for (FaceHandle fh : faces) {
        for (VertexHandle vh : m_mesh->fv_range(fh)) {
            predicate.m_vertex_filter[vh.idx()] = true;
        }
        for (EdgeHandle eh : m_mesh->fe_range(fh)) {
            predicate.m_edge_filter[eh.idx()] = true;
            
            HalfedgeHandle heh1 = m_mesh->halfedge_handle(eh, 0);
            HalfedgeHandle heh2 = m_mesh->halfedge_handle(eh, 1);
            predicate.m_halfedge_filter[heh1.idx()] = true;
            predicate.m_halfedge_filter[heh2.idx()] = true;
        }
        predicate.m_face_filter[fh.idx()] = true;
    }
    
    // boundary
    for (EdgeHandle eh : boundary_edges) {
        boundary_predicate.m_boundary_edge_filter[eh.idx()] = true;

        HalfedgeHandle heh = m_mesh->halfedge_handle(eh, 0);
        VertexHandle vh1 = m_mesh->from_vertex_handle(heh);
        VertexHandle vh2 = m_mesh->to_vertex_handle(heh);
        boundary_predicate.m_boundary_vertex_filter[vh1.idx()] = true;
        boundary_predicate.m_boundary_vertex_filter[vh2.idx()] = true;
    }
    int idx_counter = 0;
    size_t pos = predicate.m_vertex_filter.find_first();
    while (pos != boost::dynamic_bitset<>::npos) {
        VertexHandle vh = m_mesh->vertex_handle(pos);
        if (is_valid(vh) && !is_boundary(vh)) {
            non_boundary_vertex_mapping.emplace(vh, idx_counter++);
        }
        pos = predicate.m_vertex_filter.find_next(pos);
    }
    pos = predicate.m_edge_filter.find_first();
    while (pos != boost::dynamic_bitset<>::npos) {
        EdgeHandle eh = m_mesh->edge_handle(pos);
        HalfedgeHandle heh1 = m_mesh->halfedge_handle(eh, 0);
        if (is_boundary(heh1)) {
            boundary_vertex_mapping.emplace(heh1, idx_counter++);
        }
        HalfedgeHandle heh2 = m_mesh->halfedge_handle(eh, 1);
        if (is_boundary(heh2)) {
            boundary_vertex_mapping.emplace(heh2, idx_counter++);
        }
        pos = predicate.m_edge_filter.find_next(pos);
    }
    flattened_position.resize(n_duplicated_vertices());
    
    merged_edge_idx.resize(m_mesh->n_edges());
}

FilteredTriMesh::FilteredTriMesh(FilteredTriMesh&& other) noexcept :
m_mesh(other.m_mesh),
predicate(std::move(other.predicate)),
boundary_predicate(std::move(other.boundary_predicate)),
boundary_vertex_mapping(std::move(other.boundary_vertex_mapping)),
non_boundary_vertex_mapping(std::move(other.non_boundary_vertex_mapping)),
flattened_position(std::move(other.flattened_position)),
merged_subMesh_idx(std::move(other.merged_subMesh_idx)),
merged_seam_idx(std::move(other.merged_seam_idx)),
merged_edge_idx(std::move(other.merged_edge_idx)),
max_distortion(other.max_distortion),
seam_score(other.seam_score),
n_merged_seams(other.n_merged_seams)
{
}

FilteredTriMesh& FilteredTriMesh::operator=(FilteredTriMesh&& other) noexcept {
    if (this != &other) {
        m_mesh = other.m_mesh;
        other.m_mesh = nullptr;

        predicate.m_edge_filter = std::move(other.predicate.m_edge_filter);
        predicate.m_face_filter = std::move(other.predicate.m_face_filter);
        predicate.m_halfedge_filter = std::move(other.predicate.m_halfedge_filter);
        predicate.m_vertex_filter = std::move(other.predicate.m_vertex_filter);
        
        boundary_predicate.m_boundary_edge_filter = std::move(other.boundary_predicate.m_boundary_edge_filter);
        boundary_predicate.m_boundary_vertex_filter = std::move(other.boundary_predicate.m_boundary_vertex_filter);
        
        non_boundary_vertex_mapping = std::move(other.non_boundary_vertex_mapping);
        boundary_vertex_mapping = std::move(other.boundary_vertex_mapping);
        
        flattened_position = std::move(other.flattened_position);
        
        merged_subMesh_idx = std::move(other.merged_subMesh_idx);
        merged_seam_idx = std::move(other.merged_seam_idx);
        merged_edge_idx = std::move(other.merged_edge_idx);
        
        max_distortion = other.max_distortion;
        seam_score = other.seam_score;
        n_merged_seams = other.n_merged_seams;
    }
    return *this;
}

FilteredTriMesh::FilteredTriMesh(const FilteredTriMesh& other) :
m_mesh(other.m_mesh),
predicate(other.predicate),
boundary_predicate(other.boundary_predicate),
boundary_vertex_mapping(other.boundary_vertex_mapping),
non_boundary_vertex_mapping(other.non_boundary_vertex_mapping),
flattened_position(other.flattened_position),
merged_subMesh_idx(other.merged_subMesh_idx),
merged_seam_idx(other.merged_seam_idx),
merged_edge_idx(other.merged_edge_idx),
max_distortion(other.max_distortion),
seam_score(other.seam_score),
n_merged_seams(other.n_merged_seams)
{
}

FilteredTriMesh& FilteredTriMesh::operator=(const FilteredTriMesh& other) {
    m_mesh = other.m_mesh;
    predicate = other.predicate;
    boundary_predicate = other.boundary_predicate;
    boundary_vertex_mapping = other.boundary_vertex_mapping;
    non_boundary_vertex_mapping = other.non_boundary_vertex_mapping;
    flattened_position = other.flattened_position;
    merged_subMesh_idx = other.merged_subMesh_idx;
    merged_seam_idx = other.merged_seam_idx;
    merged_edge_idx = other.merged_edge_idx;
    max_distortion = other.max_distortion;
    seam_score = other.seam_score;
    n_merged_seams = other.n_merged_seams;
    return *this;
}

bool FilteredTriMesh::is_valid(VertexHandle vh) const {
    return vh.is_valid() && predicate.m_vertex_filter[vh.idx()];
}

bool FilteredTriMesh::is_valid(EdgeHandle eh) const {
    return eh.is_valid() && predicate.m_edge_filter[eh.idx()];
}

bool FilteredTriMesh::is_valid(HalfedgeHandle heh) const {
    return heh.is_valid() && predicate.m_halfedge_filter[heh.idx()];
}

bool FilteredTriMesh::is_valid(FaceHandle fh) const {
    return fh.is_valid() && predicate.m_face_filter[fh.idx()];
}

bool FilteredTriMesh::is_boundary(VertexHandle vh) const {
    assert(is_valid(vh));
    return boundary_predicate.m_boundary_vertex_filter[vh.idx()];
}
bool FilteredTriMesh::is_boundary(EdgeHandle eh) const {
    assert(is_valid(eh));
    return boundary_predicate.m_boundary_edge_filter[eh.idx()];
}
bool FilteredTriMesh::is_boundary(HalfedgeHandle heh) const {
    assert(is_valid(heh));
    if (!is_boundary(m_mesh->edge_handle(heh))) {
        return false;
    } else if (is_valid(m_mesh->face_handle(heh)) && !is_valid(m_mesh->opposite_face_handle(heh))) {
        return false;
    } else {
        return true;
    }
}
bool FilteredTriMesh::is_dual_boundary(EdgeHandle eh) const {
    assert(is_valid(eh));
    HalfedgeHandle heh = m_mesh->halfedge_handle(eh, 0);
    if (is_boundary(eh) && is_valid(m_mesh->face_handle(heh)) && is_valid(m_mesh->opposite_face_handle(heh))) {
        return true;
    } else {
        return false;
    }
}
bool FilteredTriMesh::is_dual_boundary(HalfedgeHandle heh) const {
    assert(is_valid(heh));
    if (is_boundary(m_mesh->edge_handle(heh)) && is_valid(m_mesh->face_handle(heh)) && is_valid(m_mesh->opposite_face_handle(heh))) {
        return true;
    } else {
        return false;
    }
}

unsigned long FilteredTriMesh::n_vertices() const {
    unsigned long count = 0;
    for (size_t i = 0; i < m_mesh->n_vertices(); i++) {
        if (predicate.m_vertex_filter[i]) {
            count++;
        }
    }
    return count;
}

Range<FilteredVertexIter> FilteredTriMesh::vertices() const {
    FilteredVertexIter v_it_begin(*this, predicate.m_vertex_filter, predicate.m_vertex_filter.find_first());
    FilteredVertexIter v_it_end(*this, predicate.m_vertex_filter, m_mesh->n_vertices());

    return Range<FilteredVertexIter>(v_it_begin, v_it_end);
}

Range<FilteredEdgeIter> FilteredTriMesh::edges() const {
    FilteredEdgeIter e_it_begin(*this, predicate.m_edge_filter, predicate.m_edge_filter.find_first());
    FilteredEdgeIter e_it_end(*this, predicate.m_edge_filter, m_mesh->n_edges());
    
    return Range<FilteredEdgeIter>(e_it_begin, e_it_end);
}

Range<FilteredHalfedgeIter> FilteredTriMesh::halfedges() const {
    FilteredHalfedgeIter he_it_begin(*this, predicate.m_halfedge_filter, predicate.m_halfedge_filter.find_first());
    FilteredHalfedgeIter he_it_end(*this, predicate.m_halfedge_filter, m_mesh->n_halfedges());
    
    return Range<FilteredHalfedgeIter>(he_it_begin, he_it_end);
}

Range<FilteredFaceIter> FilteredTriMesh::faces() const {
    FilteredFaceIter f_it_begin(*this, predicate.m_face_filter, predicate.m_face_filter.find_first());
    FilteredFaceIter f_it_end(*this, predicate.m_face_filter, m_mesh->n_faces());
    
    return Range<FilteredFaceIter>(f_it_begin, f_it_end);
}

Range<TriMesh::ConstFaceHalfedgeIter> FilteredTriMesh::fh_range(FaceHandle fh) const {
    return Range<TriMesh::ConstFaceHalfedgeIter>( m_mesh->cfh_begin(fh),
                                                  m_mesh->cfh_end(fh));
}

Range<FilteredVertexIter> FilteredTriMesh::boundary_vertices() const {
    FilteredVertexIter v_it_begin(*this, boundary_predicate.m_boundary_vertex_filter, boundary_predicate.m_boundary_vertex_filter.find_first());
    FilteredVertexIter v_it_end(*this, boundary_predicate.m_boundary_vertex_filter, m_mesh->n_faces());
    
    return Range<FilteredVertexIter>(v_it_begin, v_it_end);
}

Range<FilteredEdgeIter> FilteredTriMesh::boundary_edges() const {
    FilteredEdgeIter e_it_begin(*this, boundary_predicate.m_boundary_edge_filter, boundary_predicate.m_boundary_edge_filter.find_first());
    FilteredEdgeIter e_it_end(*this, boundary_predicate.m_boundary_edge_filter, m_mesh->n_edges());
    
    return Range<FilteredEdgeIter>(e_it_begin, e_it_end);
}

///@}

/// Filtered functions
///@{
VertexHandle FilteredTriMesh::vertex_handle(int idx) {
    VertexHandle vh = m_mesh->vertex_handle(idx);
    assert(is_valid(vh));
    return vh;
}

EdgeHandle FilteredTriMesh::edge_handle(HalfedgeHandle heh) const {
    assert(is_valid(heh));
    return m_mesh->edge_handle(heh);
}

HalfedgeHandle FilteredTriMesh::halfedge_handle(EdgeHandle eh, int idx) const {
    assert(is_valid(eh));
    return m_mesh->halfedge_handle(eh, idx);
}

HalfedgeHandle FilteredTriMesh::halfedge_handle(FaceHandle fh) const {
    assert(is_valid(fh));
    return m_mesh->halfedge_handle(fh);
}

HalfedgeHandle FilteredTriMesh::opposite_halfedge_handle(HalfedgeHandle heh) const {
    assert(is_valid(heh));
    return m_mesh->opposite_halfedge_handle(heh);
}

/// Only support non-boundary halfedges. For boundary halfedges, use next_boundary_halfedge_handle(HalfedgeHandle heh)
HalfedgeHandle FilteredTriMesh::next_halfedge_handle(HalfedgeHandle heh) const {
    assert(is_valid(heh));
    assert(!is_boundary(heh) || is_dual_boundary(heh));
    return m_mesh->next_halfedge_handle(heh);
}

HalfedgeHandle FilteredTriMesh::next_boundary_halfedge_handle(HalfedgeHandle heh) const {
    assert(is_valid(heh));
    assert(is_boundary(heh));

    // Find the adjacent boundary halfedge
    HalfedgeHandle next_heh = heh;
    do {
        next_heh = m_mesh->prev_halfedge_handle(m_mesh->opposite_halfedge_handle(next_heh));
    } while (next_heh != heh && (!is_valid(next_heh) || !is_boundary(m_mesh->edge_handle(next_heh))));
    
    if (!is_boundary(m_mesh->opposite_halfedge_handle(next_heh))) {
        assert(next_heh != heh && "Invalid sub-mesh, the next halfedge of the given halfedge is itself.");
    }
    // The next halfedge is reversed during rotating, flip it.
    next_heh = m_mesh->opposite_halfedge_handle(next_heh);
    if (!is_boundary(next_heh)) {
        is_boundary(next_heh);
    }
    assert(is_boundary(next_heh));
    return next_heh;
}

/// Only support non-boundary halfedges. For boundary halfedges, use prev_boundary_halfedge_handle(HalfedgeHandle heh)
HalfedgeHandle FilteredTriMesh::prev_halfedge_handle(HalfedgeHandle heh) const {
    assert(is_valid(heh));
    assert(!is_boundary(heh) || is_dual_boundary(heh));
    return m_mesh->prev_halfedge_handle(heh);
}

HalfedgeHandle FilteredTriMesh::prev_boundary_halfedge_handle(HalfedgeHandle heh) const {
    assert(is_valid(heh));
    assert(is_boundary(heh));
    
    // Find the adjacent boundary halfedge
    HalfedgeHandle prev_heh = heh;
    do {
        prev_heh = m_mesh->next_halfedge_handle(m_mesh->opposite_halfedge_handle(prev_heh));
    } while (prev_heh != heh && (!is_valid(prev_heh) || !is_boundary(m_mesh->edge_handle(prev_heh))));
    
    if (!is_boundary(m_mesh->opposite_halfedge_handle(prev_heh))) {
        assert(prev_heh != heh && "Invalid sub-mesh, the prev halfedge of the given halfedge is itself.");
    }
    // The prev halfedge is reversed during rotating, flip it.
    prev_heh = m_mesh->opposite_halfedge_handle(prev_heh);
    assert(is_boundary(prev_heh));
    return prev_heh;
}

VertexHandle FilteredTriMesh::from_vertex_handle(HalfedgeHandle heh) const {
    assert(is_valid(heh));
    return m_mesh->from_vertex_handle(heh);
}

VertexHandle FilteredTriMesh::to_vertex_handle(HalfedgeHandle heh) const {
    assert(is_valid(heh));
    return m_mesh->to_vertex_handle(heh);
}

size_t FilteredTriMesh::valence(VertexHandle vh) const {
    size_t count = 0;
    for (VertexHandle vh_neighbor : m_mesh->vv_range(vh)) {
        if (is_valid(vh_neighbor)) {
            count++;
        }
    }
    return count;
}

TriMesh::Point& FilteredTriMesh::point(VertexHandle vh) {
    assert(is_valid(vh));
    return m_mesh->point(vh);
}

const TriMesh::Color& FilteredTriMesh::color(FaceHandle fh) const {
    assert(is_valid(fh));
    return m_mesh->color(fh);
}

double FilteredTriMesh::calc_edge_length(EdgeHandle eh) const {
    assert(is_valid(eh));
    return m_mesh->calc_edge_length(eh);
}

double FilteredTriMesh::calc_edge_length(HalfedgeHandle heh) const {
    assert(is_valid(heh));
    return m_mesh->calc_edge_length(heh);
}

double FilteredTriMesh::calc_sector_angle(HalfedgeHandle heh) const {
    assert(is_valid(heh));
    return m_mesh->calc_sector_angle(heh);
}

double FilteredTriMesh::calc_sector_area(HalfedgeHandle heh) const {
    assert(is_valid(heh));
    return m_mesh->calc_sector_area(heh);
}
///@}

/// Additional operations
///@{
TriMesh::Point& FilteredTriMesh::flattened_non_boundary_point(VertexHandle vh) {
    int idx = get_non_boundary_vertex_idx(vh);
    assert(idx >= 0);
    return flattened_position[idx];
}

TriMesh::Point& FilteredTriMesh::flattened_boundary_point(HalfedgeHandle heh) {
    int idx = get_boundary_vertex_idx(heh);
    assert(idx >= 0);
    return flattened_position[idx];
}

double FilteredTriMesh::getSumInnerAngle(HalfedgeHandle from, HalfedgeHandle to) {
    assert(is_valid(from));
    assert(is_valid(to));
    assert(to_vertex_handle(from) == from_vertex_handle(to));
    
    // We calculate both half sum and full sum (a circle) of angle
    double sum_angle = 0, half_sum_angle = 0;
    HalfedgeHandle current_heh = from, target_heh = to;
    do {
        sum_angle += m_mesh->calc_sector_angle(current_heh);
        HalfedgeHandle next_heh = m_mesh->opposite_halfedge_handle(
                                                                   m_mesh->next_halfedge_handle(current_heh));
        if (m_mesh->edge_handle(next_heh) == m_mesh->edge_handle(to)) {
            half_sum_angle = sum_angle;
        }
        current_heh = next_heh;
    } while (current_heh != from);
    
    // output depends on whether "from" is boundary or not
    if (is_boundary(from) && !is_dual_boundary(from)) {
        return sum_angle - half_sum_angle;
    } else {
        return half_sum_angle;
    }
}

/**
 Return the index of the given non-boundary vertex in this filtered mesh.
 This index is different from vh.idx() since this mesh skipped some vertices.
 @param <#parameter#>
 @return The index of the given vertex handle.
 @retval <#meaning#>
 */
int FilteredTriMesh::get_non_boundary_vertex_idx(VertexHandle vh) {
    assert(!is_boundary(vh));
    assert(non_boundary_vertex_mapping.find(vh) != non_boundary_vertex_mapping.end());
    return non_boundary_vertex_mapping[vh];
}

/**
 Return the "from" vertex index of the given boundary halfedge in this filtered mesh.
 This index is different from from_vertex_handle(vh).idx() since this mesh skipped some vertices.
 @param <#parameter#>
 @return The "from" vertex index of the given halfedge handle.
 @retval <#meaning#>
 */
int FilteredTriMesh::get_boundary_vertex_idx(HalfedgeHandle heh) {
    assert(is_valid(heh));
    
    // Make sure the "from" vertex is a boundary vertex
    VertexHandle vh = m_mesh->from_vertex_handle(heh);
    assert(is_boundary(vh));
    
    if (boundary_vertex_mapping.find(heh) != boundary_vertex_mapping.end()) {
        return boundary_vertex_mapping[heh];
    } else {
        // If not found, this is not a boundary halfedge. However, one of the neighboring halfedges is.
        // Just search through them.
        // boundary vertex on non-boundary halfedge, we need to search for neighbor boundary halfedge
        HalfedgeHandle nearby_boundary_halfedge = get_nearby_boundary_halfedge(heh);
        
        assert(boundary_vertex_mapping.find(nearby_boundary_halfedge) != boundary_vertex_mapping.end());
        return boundary_vertex_mapping[nearby_boundary_halfedge];
    }
    return -1;
}

int FilteredTriMesh::n_duplicated_vertices() {
    return non_boundary_vertex_mapping.size() + boundary_vertex_mapping.size();
}

TriMesh* FilteredTriMesh::mesh() const {
    return m_mesh;
}
