#include "PlushPatternGenerator.hh"

OpenMesh::VPropHandleT<double> PlushPatternGenerator::minCurvatureHandle;
OpenMesh::VPropHandleT<double> PlushPatternGenerator::maxCurvatureHandle;
OpenMesh::VPropHandleT<OpenMesh::Vec3d> PlushPatternGenerator::minCurvatureDirectionHandle;
OpenMesh::VPropHandleT<OpenMesh::Vec3d> PlushPatternGenerator::maxCurvatureDirectionHandle;

OpenMesh::EPropHandleT<double> PlushPatternGenerator::edgeWeightHandle;
OpenMesh::MPropHandleT< std::map<std::pair<VertexHandle, VertexHandle>, double> > PlushPatternGenerator::geodesicDistanceHandle;
OpenMesh::MPropHandleT< std::map<std::pair<VertexHandle, VertexHandle>, std::vector<VertexHandle> > > PlushPatternGenerator::geodesicPathHandle;

OpenMesh::MPropHandleT<Skeleton*> PlushPatternGenerator::skeletonHandle;
OpenMesh::VPropHandleT<double*> PlushPatternGenerator::bonesWeightHandle;

PlushPatternGenerator::PlushPatternGenerator(TriMesh *mesh, QString meshName) : m_mesh(mesh), m_meshName(meshName) {
    // convert m_mesh from OpenFlipper to CGAL
    CGAL_Polyhedron_helper::convert_OpenMesh_to_CGAL(m_mesh, m_polyhedron, m_verticesMapping);
    
    // setup properties
    initProperties();
    
    loadCurvature();
    
    loadSkeleton();
    
    loadBoneWeight();
}

PlushPatternGenerator::~PlushPatternGenerator() {
    uninitProperties();
}

void PlushPatternGenerator::initProperties() {
    m_mesh->add_property(minCurvatureHandle, "Min Curvature");
    m_mesh->add_property(maxCurvatureHandle, "Max Curvature");
    m_mesh->add_property(minCurvatureDirectionHandle, "Min curvature direction");
    m_mesh->add_property(maxCurvatureDirectionHandle, "Max curvature direction");
    
    m_mesh->add_property(edgeWeightHandle, "Edge weight");
    m_mesh->add_property(geodesicDistanceHandle, "Geodesic distance between vertices pair");
    m_mesh->add_property(geodesicPathHandle, "Geodesic path between vertices pair");
    
    m_mesh->add_property(skeletonHandle, "Skeleton");
    m_mesh->add_property(bonesWeightHandle, "Bone weights for each vertex");
    
    //    for (VertexIter v_it = mesh->vertices_begin(); v_it != mesh->vertices_end(); v_it++) {
    //    }
    for (EdgeIter e_it = m_mesh->edges_begin(); e_it != m_mesh->edges_end(); e_it++) {
        m_mesh->property(edgeWeightHandle, *e_it) = -1;
    }
}

void PlushPatternGenerator::uninitProperties() {
    if (m_mesh->property(skeletonHandle) != NULL) {
        delete m_mesh->property(skeletonHandle);
    }
    for (VertexIter v_it = m_mesh->vertices_begin(); v_it != m_mesh->vertices_end(); v_it++) {
        if (m_mesh->property(bonesWeightHandle, *v_it) != NULL) {
            delete[] m_mesh->property(bonesWeightHandle, *v_it);
        }
    }
}

void PlushPatternGenerator::cancelJob() {
    isJobCanceled = true;
}

bool PlushPatternGenerator::getEdge(TriMesh *mesh, EdgeHandle &_eh, VertexHandle v1, VertexHandle v2) {
    return getEdge(mesh, _eh, v1.idx(), v2.idx());
}

bool PlushPatternGenerator::getEdge(TriMesh *mesh, EdgeHandle &_eh, int v1No, int v2No) {
    VertexHandle v1 = mesh->vertex_handle(v1No);
    VertexHandle v2 = mesh->vertex_handle(v2No);
    for (TriMesh::VertexEdgeIter ve_it = mesh->ve_iter(v1); ve_it; ve_it++)
    {
        HalfedgeHandle heh = mesh->halfedge_handle(*ve_it, 0);
        if (mesh->from_vertex_handle(heh) == v2
            ||  mesh->to_vertex_handle(heh) == v2)
        {
            _eh = *ve_it;
            return true;
        }
    }
    return false;
}
