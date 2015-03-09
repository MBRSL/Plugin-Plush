#include "PlushPatternGenerator.hh"
#include <Eigen/Sparse>

/**
 * @brief This function divide a set of edges (with loops) into several loops.
 * @param loops Result loops will be stored here.
 * @param selectedVertices (This parameter should be compute internally) The set of intersection point of m_spanningTree.
 */
void PlushPatternGenerator::getLoops(std::vector< std::vector<HalfedgeHandle> > &loops, std::vector<VertexHandle> &selectedVertices) {
    std::set<HalfedgeHandle> isVisited;
    for (size_t i = 0; i < selectedVertices.size(); i++) {
        VertexHandle startV = selectedVertices[i];
        
        // Starting from a intersection point, search for its adjacenct edges
        for (TriMesh::VertexOHalfedgeIter voh_it = m_mesh->voh_iter(startV); voh_it; voh_it++) {
            // Skip visited edges and non-spanning-tree edges
            if (std::find(isVisited.begin(), isVisited.end(), *voh_it) != isVisited.end()
                ||  std::find(m_spanningTree.begin(), m_spanningTree.end(), m_mesh->edge_handle(*voh_it)) == m_spanningTree.end()) {
                continue;
            }
            
            // We got a non-visited edge, lets start searching from this edge
            HalfedgeHandle start_heh = *voh_it;
            HalfedgeHandle current_heh = start_heh;
            isVisited.insert(start_heh);
            
            std::vector<HalfedgeHandle> loop;
            loop.push_back(current_heh);
            do {
                assert(!m_mesh->is_boundary(current_heh));
                
                HalfedgeHandle next_heh = current_heh;
                // Find next edge which is in spanning tree
                do {
                    // While there are possibly many adjacent edges to choose, we use next->opposite to iterator through them
                    // so that we will visit them all clockwise without missing any one.
                    next_heh = m_mesh->opposite_halfedge_handle(
                                                                m_mesh->next_halfedge_handle(next_heh));
                    
                    // This means there are no other edges except current_heh (which is where you came from)
                    // In this case, just go back and treat this path as a dart
                    if (next_heh == current_heh) {
                        // Do nothing, it will be flipped later.
                    }
                } while (std::find(m_spanningTree.begin(), m_spanningTree.end(), m_mesh->edge_handle(next_heh)) == m_spanningTree.end());
                
                // Flip it to point to the next vertex
                current_heh = m_mesh->opposite_halfedge_handle(next_heh);
                isVisited.insert(current_heh);
                
                loop.push_back(current_heh);
            } while (current_heh != start_heh);
            
            loops.push_back(loop);
        }
    }
}