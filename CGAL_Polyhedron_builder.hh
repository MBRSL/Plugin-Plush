//
//  CGAL_Polyhedron.h
//  OpenFlipper
//
//  Created by 饒亢 on 2015/1/8.
//
//

#ifndef PLUSHPLUGIN_CGAL_POLYHEDRON
#define PLUSHPLUGIN_CGAL_POLYHEDRON

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/iterator.h>

#include <ObjectTypes/TriangleMesh/TriangleMesh.hh>

typedef CGAL::Simple_cartesian<double>     Kernel;
typedef CGAL::Polyhedron_3<Kernel, CGAL::Polyhedron_items_with_id_3>         Polyhedron;

// A modifier creating a triangle with the incremental builder.
template<class HDS>
class CGAL_Polyhedron_builder : public CGAL::Modifier_base<HDS> {
private:
    
public:
    std::vector<double> coords;
    std::vector<int>    tris;
    
    void operator()( HDS& hds) {
        typedef typename HDS::Vertex   Vertex;
        typedef typename Vertex::Point Point;
        
        
        // create a cgal incremental builder
        CGAL::Polyhedron_incremental_builder_3<HDS> B( hds, true);
        B.begin_surface( coords.size()/3, tris.size()/3 );
        
        // add the polyhedron vertices
        for( int i=0; i<(int)coords.size(); i+=3 ){
            B.add_vertex( Point( coords[i+0], coords[i+1], coords[i+2] ) );
        }
        
        // add the polyhedron triangles
        for( int i=0; i<(int)tris.size(); i+=3 ){
            B.begin_facet();
            B.add_vertex_to_facet( tris[i+0] );
            B.add_vertex_to_facet( tris[i+1] );
            B.add_vertex_to_facet( tris[i+2] );
            B.end_facet();
        }
        
        // finish up the surface
        B.end_surface();
    }
};

void translate_mesh_from_OpenMesh_to_CGAL(TriMesh *mesh, Polyhedron *P) {
    CGAL_Polyhedron_builder<Polyhedron::HalfedgeDS> builder;
    
    for (VertexIter v_it = mesh->vertices_begin(); v_it != mesh->vertices_end(); v_it++)
    {
        TriMesh::Point p = mesh->point(*v_it);
        builder.coords.push_back(p[0]);
        builder.coords.push_back(p[1]);
        builder.coords.push_back(p[2]);
    }
    for (FaceIter f_it = mesh->faces_begin(); f_it != mesh->faces_end(); f_it++)
    {
        for (TriMesh::FaceVertexIter fv_it = mesh->fv_iter(*f_it); fv_it; fv_it++)
        {
            builder.tris.push_back(fv_it->idx());
        }
    }
    P->delegate(builder);
}
#endif