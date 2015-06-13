#include "PlushPatternGenerator.hh"

OpenMesh::EPropHandleT<int> PlushPatternGenerator::segment_no_handle;

OpenMesh::MPropHandleT< std::set<VertexHandle> > PlushPatternGenerator::intersection_points_handle;

OpenMesh::MPropHandleT< std::map< std::set<HalfedgeHandle>, double> > PlushPatternGenerator::joint_boundary_area_handle;
OpenMesh::MPropHandleT< std::map< std::set<HalfedgeHandle>, double> > PlushPatternGenerator::joint_boundary_distortion_handle;
OpenMesh::VPropHandleT<double> PlushPatternGenerator::meanCurvatureHandle;
OpenMesh::VPropHandleT<double> PlushPatternGenerator::gaussianCurvatureHandle;

OpenMesh::EPropHandleT<double> PlushPatternGenerator::edgeWeightHandle;
OpenMesh::MPropHandleT< std::map<std::pair<VertexHandle, VertexHandle>, double> > PlushPatternGenerator::geodesicDistanceHandle;
OpenMesh::MPropHandleT< std::map<std::pair<VertexHandle, VertexHandle>, std::vector<VertexHandle> > > PlushPatternGenerator::geodesicPathHandle;

OpenMesh::MPropHandleT<Skeleton*> PlushPatternGenerator::skeletonHandle;
OpenMesh::VPropHandleT<double*> PlushPatternGenerator::bonesWeightHandle;

OpenMesh::MPropHandleT< std::vector<TriMesh> > PlushPatternGenerator::flattenedMeshesHandle;
OpenMesh::FPropHandleT<double> PlushPatternGenerator::distortionFHandle;
OpenMesh::VPropHandleT<double> PlushPatternGenerator::distortionVHandle;

OpenMesh::VPropHandleT<VertexHandle> PlushPatternGenerator::getInverseMappingHandle(TriMesh *mesh) {
    OpenMesh::VPropHandleT<VertexHandle> inverse_mappingHandle;
    if (!mesh->get_property_handle(inverse_mappingHandle, "inverse_mapping")) {
        mesh->add_property(inverse_mappingHandle, "inverse_mapping");
        for (VertexHandle v : mesh->vertices()) {
            mesh->property(inverse_mappingHandle, v).invalidate();
        }
    }
    return inverse_mappingHandle;
}

OpenMesh::MPropHandleT< std::set<EdgeHandle> > PlushPatternGenerator::getSeamsHandle(TriMesh *mesh) {
    OpenMesh::MPropHandleT< std::set<EdgeHandle> > seamsHandle;
    if(!mesh->get_property_handle(seamsHandle, "seams")) {
        mesh->add_property(seamsHandle, "seams");
    }
    return seamsHandle;
}

std::set<EdgeHandle>* PlushPatternGenerator::getSeams() {
    OpenMesh::MPropHandleT< std::set<EdgeHandle> > seamsHandle = getSeamsHandle(m_mesh);
    if (m_mesh->property(seamsHandle).size() == 0) {
        return NULL;
    } else {
        return &m_mesh->property(seamsHandle);
    }
}

PlushPatternGenerator::PlushPatternGenerator(TriMesh *mesh, QString meshName) : m_mesh(mesh), m_meshName(meshName) {
    // setup properties
    initProperties();
    
    loadCurvature();
    
    loadSkeleton();
    
    loadBoneWeight();
    
    loadGeodesic();
    
//    load_seams();
}

PlushPatternGenerator::~PlushPatternGenerator() {
    uninitProperties();
}

void PlushPatternGenerator::initProperties() {
    m_mesh->add_property(segment_no_handle, "Segment number");
    
    m_mesh->add_property(intersection_points_handle, "Set of intersection point of union of seams");
    m_mesh->add_property(joint_boundary_area_handle, "Sum of area for joint boundary");
    m_mesh->add_property(joint_boundary_distortion_handle, "Distortion for joint boundary");
    
    m_mesh->add_property(meanCurvatureHandle, "Mean Curvature");
    m_mesh->add_property(gaussianCurvatureHandle, "Gaussian Curvature");
    
    m_mesh->add_property(edgeWeightHandle, "Edge weight");
    m_mesh->add_property(geodesicDistanceHandle, "Geodesic distance between vertices pair");
    m_mesh->add_property(geodesicPathHandle, "Geodesic path between vertices pair");
    
    m_mesh->add_property(skeletonHandle, "Skeleton");
    m_mesh->add_property(bonesWeightHandle, "Bone weights for each vertex");
    
    m_mesh->add_property(flattenedMeshesHandle, "flattened meshes");

    m_mesh->add_property(distortionVHandle, "Flattening distortion of vertices");
    m_mesh->add_property(distortionFHandle, "Flattening distortion of faces");
    
    for (VertexHandle v : m_mesh->vertices()) {
        m_mesh->property(distortionVHandle, v) = 0;
        m_mesh->property(meanCurvatureHandle, v) = 0;
        m_mesh->property(gaussianCurvatureHandle, v) = 0;
        m_mesh->property(bonesWeightHandle, v) = NULL;
    }

    m_mesh->request_edge_colors();
    for (EdgeHandle eh : m_mesh->edges()) {
        m_mesh->property(edgeWeightHandle, eh) = -1;
        m_mesh->property(segment_no_handle, eh) = 0;
        m_mesh->set_color(eh, TriMesh::Color(1,1,1,0));
    }
    
    for (FaceHandle f : m_mesh->faces()) {
        m_mesh->property(distortionFHandle, f) = 0;
    }
}

void PlushPatternGenerator::uninitProperties() {
    m_mesh->remove_property(segment_no_handle);
    
    m_mesh->remove_property(intersection_points_handle);

    m_mesh->remove_property(joint_boundary_area_handle);
    m_mesh->remove_property(joint_boundary_distortion_handle);
    
    m_mesh->remove_property(meanCurvatureHandle);
    m_mesh->remove_property(gaussianCurvatureHandle);
    
    m_mesh->remove_property(edgeWeightHandle);
    m_mesh->remove_property(geodesicDistanceHandle);
    m_mesh->remove_property(geodesicPathHandle);
    
    if (m_mesh->property(skeletonHandle) != NULL) {
        delete m_mesh->property(skeletonHandle);
    }
    for (VertexIter v_it = m_mesh->vertices_begin(); v_it != m_mesh->vertices_end(); v_it++) {
        if (m_mesh->property(bonesWeightHandle, *v_it) != NULL) {
            delete[] m_mesh->property(bonesWeightHandle, *v_it);
        }
    }
    m_mesh->remove_property(skeletonHandle);
    m_mesh->remove_property(bonesWeightHandle);
    
    OpenMesh::MPropHandleT< std::vector<EdgeHandle> > seamsHandle;
    if(!m_mesh->get_property_handle(seamsHandle, "seams")) {
        m_mesh->remove_property(seamsHandle);
    }
    m_mesh->remove_property(flattenedMeshesHandle);

    m_mesh->remove_property(distortionVHandle);
    m_mesh->remove_property(distortionFHandle);
    
}

void PlushPatternGenerator::cancelJob() {
    isJobCanceled = true;
}

bool PlushPatternGenerator::getHalfedge(const TriMesh *mesh, HalfedgeHandle &heh, VertexHandle from, VertexHandle to) {
    for (TriMesh::ConstVertexEdgeIter cve_it = mesh->cve_iter(from); cve_it; cve_it++)
    {
        // Return halfedge with correct direction
        HalfedgeHandle _heh = mesh->halfedge_handle(*cve_it, 0);
        if (mesh->from_vertex_handle(_heh) == from
            &&  mesh->to_vertex_handle(_heh) == to)
        {
            heh = _heh;
            return true;
        } else if (mesh->from_vertex_handle(_heh) == to
                   &&  mesh->to_vertex_handle(_heh) == from) {
            heh = mesh->halfedge_handle(*cve_it, 1);
            return true;
        }
    }
    return false;
}

bool PlushPatternGenerator::getHalfedge(const TriMesh *mesh, HalfedgeHandle &heh, int fromNo, int toNo) {
    VertexHandle from = mesh->vertex_handle(fromNo);
    VertexHandle to = mesh->vertex_handle(toNo);
    return getHalfedge(mesh, heh, from, to);
}

bool PlushPatternGenerator::getEdge(const TriMesh *mesh, EdgeHandle &eh, VertexHandle v1, VertexHandle v2) {
    HalfedgeHandle _heh;
    bool result = getHalfedge(mesh, _heh, v1, v2);
    if (result) {
        eh = mesh->edge_handle(_heh);
        return true;
    } else {
        return false;
    }
}

bool PlushPatternGenerator::getEdge(const TriMesh *mesh, EdgeHandle &eh, int v1No, int v2No) {
    VertexHandle v1 = mesh->vertex_handle(v1No);
    VertexHandle v2 = mesh->vertex_handle(v2No);
    return getEdge(mesh, eh, v1, v2);
}

bool PlushPatternGenerator::getFace(const TriMesh *mesh, FaceHandle &fh, VertexHandle v1, VertexHandle v2, VertexHandle v3) {
    for (TriMesh::ConstVertexFaceIter cvf_it = mesh->cvf_begin(v1); cvf_it; cvf_it++) {
        HalfedgeHandle heh = mesh->halfedge_handle(*cvf_it);
        VertexHandle _v1 = mesh->from_vertex_handle(heh);
        VertexHandle _v2 = mesh->to_vertex_handle(heh);
        VertexHandle _v3 = mesh->to_vertex_handle(mesh->next_halfedge_handle(heh));

        if ((_v1 == v1 && _v2 == v2 && _v3 == v3)
        ||  (_v1 == v1 && _v2 == v3 && _v3 == v2)
        ||  (_v1 == v2 && _v2 == v1 && _v3 == v3)
        ||  (_v1 == v2 && _v2 == v3 && _v3 == v1)
        ||  (_v1 == v3 && _v2 == v1 && _v3 == v2)
        ||  (_v1 == v3 && _v2 == v2 && _v3 == v1)) {
            fh = *cvf_it;
            return true;
        }
    }
    
    // Not found
    return false;
}

std::vector<TriMesh>* PlushPatternGenerator::getFlattenedMeshes() {
    if (!flattenedMeshesHandle.is_valid() || m_mesh->property(flattenedMeshesHandle).size() == 0) {
        return NULL;
    } else {
        return &m_mesh->property(flattenedMeshesHandle);
    }
}

VertexHandle PlushPatternGenerator::get_original_handle(TriMesh *mesh, const VertexHandle vh) const {
    if (mesh == m_mesh) {
        return vh;
    } else {
        OpenMesh::VPropHandleT<VertexHandle> inverse_mapping = getInverseMappingHandle(mesh);
        return mesh->property(inverse_mapping, vh);
    }
}

EdgeHandle PlushPatternGenerator::get_original_handle(TriMesh *mesh, const EdgeHandle eh) const {
    if (mesh == m_mesh) {
        return eh;
    } else {
        OpenMesh::VPropHandleT<VertexHandle> inverse_mapping = getInverseMappingHandle(mesh);
        HalfedgeHandle heh = mesh->halfedge_handle(eh, 0);
        VertexHandle v1 = mesh->property(inverse_mapping, mesh->from_vertex_handle(heh));
        VertexHandle v2 = mesh->property(inverse_mapping, mesh->to_vertex_handle(heh));
        EdgeHandle original_eh;
        bool found = getEdge(m_mesh, original_eh, v1, v2);
        assert(found);
        return original_eh;
    }
}

HalfedgeHandle PlushPatternGenerator::get_original_handle(TriMesh *mesh, const HalfedgeHandle heh) const {
    if (mesh == m_mesh) {
        return heh;
    } else {
        OpenMesh::VPropHandleT<VertexHandle> inverse_mapping = getInverseMappingHandle(mesh);
        VertexHandle v1 = mesh->property(inverse_mapping, mesh->from_vertex_handle(heh));
        VertexHandle v2 = mesh->property(inverse_mapping, mesh->to_vertex_handle(heh));
        HalfedgeHandle original_heh;
        bool found = getHalfedge(m_mesh, original_heh, v1, v2);
        assert(found);
        return original_heh;
    }
}

FaceHandle PlushPatternGenerator::get_original_handle(TriMesh *mesh, const FaceHandle fh) const {
    if (mesh == m_mesh) {
        return fh;
    } else {
        OpenMesh::VPropHandleT<VertexHandle> inverse_mapping = getInverseMappingHandle(mesh);
        TriMesh::ConstFaceVertexIter cfv_it = mesh->cfv_iter(fh);
        VertexHandle v1 = mesh->property(inverse_mapping, *cfv_it++);
        VertexHandle v2 = mesh->property(inverse_mapping, *cfv_it++);
        VertexHandle v3 = mesh->property(inverse_mapping, *cfv_it++);
        FaceHandle original_fh;
        bool found = getFace(m_mesh, original_fh, v1, v2, v3);
        assert(found);
        return original_fh;
    }
}