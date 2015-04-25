#ifndef PLUSH_PATTERN_GENERATOR_HH
#define PLUSH_PATTERN_GENERATOR_HH

// This one must be included before any other OpenFlipper headers, otherwise it will cause weird runtime errors.
#include <OpenFlipper/BasePlugin/BaseInterface.hh>

#include "SuperDeform/Skeleton.hh"

#include "Common.hh"

#include <Eigen/sparse>

#include <MeshTools/MeshSelectionT.hh>

#include <QObject>

/**
 * @brief This class contains algorithms and utilities for plush pattern generation.
 */

class PlushPatternGenerator : public QObject
{
    Q_OBJECT
signals:
    
    /// Emitted when error occurred, or other informations.
    void log(int logType, QString msg);
    
    /// Use this signal to inform user current job progress.
    void setJobState(int val);
    
public:
    /// Used for log
    ///@{
    static const int LOGERR = 0;
    static const int LOGINFO = 1;
    ///@}
    
    /// @name Curvature property handle
    ///@{
    /// Curvature value/direction of each vertex. It's empty before curvature calculation.
    static OpenMesh::VPropHandleT<double> maxCurvatureHandle;
    ///@}
    
    /// @name Geodesic distance property handle
    ///@{
    /// Geodesic distance/path from one vertex to another. It's empty before geodesic calculation.
    static OpenMesh::EPropHandleT<double> edgeWeightHandle;
    static OpenMesh::MPropHandleT< std::map<std::pair<VertexHandle, VertexHandle>, double> > geodesicDistanceHandle;
    static OpenMesh::MPropHandleT< std::map<std::pair<VertexHandle, VertexHandle>, std::vector<VertexHandle> > > geodesicPathHandle;
    ///@}
    
    /// @name Skeleton property handle
    ///@{
    
    /// The skeleton of this mesh. It's built during mesh loading unless error occurs
    static OpenMesh::MPropHandleT<Skeleton*> skeletonHandle;
    /// The weight of corresponding bones for each vertex
    static OpenMesh::VPropHandleT<double*> bonesWeightHandle;
    ///@}
    
    /// @name Flattening-related handle
    ///@{
    /// This stores flattened meshes.
    static OpenMesh::MPropHandleT< std::vector<TriMesh> > flattenedMeshesHandle;
    /// Distortion indicator for visualization
    static OpenMesh::FPropHandleT<double> distortionFHandle;
    static OpenMesh::VPropHandleT<double> distortionVHandle;
    /// Mapping flattened vertex handle back to original vertex handle
    static OpenMesh::VPropHandleT<VertexHandle> getInverseMappingHandle(TriMesh *mesh);
    static OpenMesh::MPropHandleT< std::set<EdgeHandle> > getSeamsHandle(TriMesh *mesh);
    ///@}
    
    static bool getHalfedge(const TriMesh *mesh, HalfedgeHandle &heh, int fromNo, int toNo);
    static bool getHalfedge(const TriMesh *mesh, HalfedgeHandle &heh, VertexHandle from, VertexHandle to);
    static bool getEdge(const TriMesh *mesh, EdgeHandle &eh, int v1No, int v2No);
    static bool getEdge(const TriMesh *mesh, EdgeHandle &eh, VertexHandle v1, VertexHandle v2);
    static bool getFace(const TriMesh *mesh, FaceHandle &fh, VertexHandle v1, VertexHandle v2, VertexHandle v3);
    
    /// Calculate the sum of inner (clockwise) angles by iterating from one halfedge to another.
    static double getSumInnerAngle(const TriMesh *mesh, HalfedgeHandle heh1, HalfedgeHandle heh2);
    static double getSumInnerAngle(const TriMesh *mesh, VertexHandle v);
    
    /// Get boundary for a given opened mesh
    static bool getBoundaryOfOpenedMesh(std::vector< std::vector<HalfedgeHandle> > &boundaries, const TriMesh *mesh, bool getInteriorHalfedge);
    
    /** Expand selection by n-ring connectivity **/
    static void expandVertice(const TriMesh *mesh, const VertexHandle centerV, std::set<VertexHandle> &verticesSelection, int n, double maxDistance);
    static void expandVertice(const TriMesh *mesh, std::set<VertexHandle> &verticesSelection, int n);
    static void shrinkVertice(const TriMesh *mesh, std::set<VertexHandle> &verticesSelection, int n);
    
    static double calcArea(OpenMesh::Vec3d p1, OpenMesh::Vec3d p2, OpenMesh::Vec3d p3);
    
    enum InteriorParameterizationMethod { Barycentric, Conformal, MeanValue };
    
    PlushPatternGenerator(TriMesh *mesh, QString meshName);
    ~PlushPatternGenerator();
    
    void loadCurvature();
    
    bool loadSkeleton();
    bool loadBoneWeight();
    bool saveBoneWeight();
    
    bool calcCurvature();
    
    bool calcSelection(std::vector<VertexHandle> &targetVertices, double threshold=0.001);
    
    void calcGeodesic(std::vector<VertexHandle> targetVertices);
    bool loadGeodesic();
    bool saveGeodesic(std::vector<VertexHandle> selectedVertices);
    
    bool calcSeams(std::vector<VertexHandle> selectedVertices, int limitNum = 0, bool elimination = false, bool allPaths = true);
    bool calcCircularSeams(TriMesh *mesh);
    
    bool calcSkeletonWeight();
    
    // Used in selection
    void calc_parameterization_weight_matrix(TriMesh *mesh, Eigen::SparseMatrix<double> &M, InteriorParameterizationMethod method);
    
    /// Flatten 3D loops into 2D loops using LPFB.
    bool calcFlattenedGraph();
    
    std::vector<TriMesh>* getFlattenedMeshes();
    std::set<EdgeHandle>* getSeams();
    
private:
    TriMesh *m_mesh;
    QString m_meshName;
    
    bool isJobCanceled;
    
    bool isIntersected(std::vector<VertexHandle> path1, std::vector<VertexHandle> path2);
    
    bool splitWithBoundary(std::vector<TriMesh> *subMeshes, std::set<EdgeHandle> *seams);
    
    /// @name Flattening
    ///@{
    bool calcLPFB(TriMesh *mesh, std::map<VertexHandle, OpenMesh::Vec3d> *boundaryPosition);
    bool calcInteriorPoints(TriMesh *mesh, std::map<VertexHandle, OpenMesh::Vec3d> *boundaryPosition);
    /// Calculate the distortion between original 3D mesh and flattened graph, the result is stored in distortionHandle
    void calcDistortion(std::vector<TriMesh> *flattenedMeshes);
    /// Organize 2D loops to prevent overlapping while minimize bounding area.
    bool packFlattenedGraph(std::vector<TriMesh> *flattenedMeshes, const int nColumns=4);
    ///@}
    
    void initProperties();
    void uninitProperties();
    
    VertexHandle get_original_handle(TriMesh *mesh, const VertexHandle vh) const;
    EdgeHandle get_original_handle(TriMesh *mesh, const EdgeHandle eh) const;
    HalfedgeHandle get_original_handle(TriMesh *mesh, const HalfedgeHandle eh) const;
    FaceHandle get_original_handle(TriMesh *mesh, const FaceHandle fh) const;
    
    /// Utils
    template<class T>
    void markSelection(const std::vector< std::set<T> > &group_of_handles, TriMesh *mesh) {
        OpenMesh::VPropHandleT<VertexHandle> inverseMapping = getInverseMappingHandle(mesh);
        for (auto group_it = group_of_handles.begin(); group_it != group_of_handles.end(); group_it++) {
            markSelection(*group_it, mesh, inverseMapping);
        }
    }
    
    template<class T>
    inline void markSelection(const std::set<T> &handles, TriMesh *mesh, OpenMesh::VPropHandleT<VertexHandle> &inverseMapping) {
        std::vector<T> tmp(handles.begin(), handles.end());
        markSelection(tmp, mesh, inverseMapping);
    }
    inline void markSelection(const std::vector<FaceHandle> &faces, TriMesh *mesh, OpenMesh::VPropHandleT<VertexHandle> &inverseMapping) {
        std::vector<int> facesId;
        for (auto f_it = faces.begin(); f_it != faces.end(); f_it++) {
            VertexHandle originalV[3];
            int count = 0;
            for (TriMesh::ConstFaceVertexIter cfv_it = mesh->cfv_iter(*f_it); cfv_it; cfv_it++) {
                originalV[count++] = mesh->property(inverseMapping, *cfv_it);
            }
            FaceHandle originalF;
            assert(getFace(m_mesh, originalF, originalV[0], originalV[1], originalV[2]));
            facesId.push_back(originalF.idx());
        }
        MeshSelection::selectFaces(m_mesh, facesId);
    }
    
    inline void markSelection(const std::vector<VertexHandle> &vertices, TriMesh *mesh, OpenMesh::VPropHandleT<VertexHandle> &inverseMapping) {
        std::vector<int> verticesId;
        for (auto v_it = vertices.begin(); v_it != vertices.end(); v_it++) {
            VertexHandle originalV = mesh->property(inverseMapping, *v_it);
            verticesId.push_back(originalV.idx());
        }
        MeshSelection::selectVertices(m_mesh, verticesId);
    }
    
    inline void markSelection(const std::vector<EdgeHandle> &edes, TriMesh *mesh, OpenMesh::VPropHandleT<VertexHandle> &inverseMapping) {
        std::vector<int> edgesId;
        for (auto e_it = edes.begin(); e_it != edes.end(); e_it++) {
            HalfedgeHandle heh = mesh->halfedge_handle(*e_it, 0);
            VertexHandle originalV1 = mesh->property(inverseMapping, mesh->from_vertex_handle(heh));
            VertexHandle originalV2 = mesh->property(inverseMapping, mesh->to_vertex_handle(heh));
            EdgeHandle original_eh;
            assert(getEdge(m_mesh, original_eh, originalV1, originalV2));
            edgesId.push_back(original_eh.idx());
        }
        MeshSelection::selectEdges(m_mesh, edgesId);
    }
    
    public slots:
    /// Cancel time-consuming works. (calcGeodesic, calcCurvature, etc.)
    void cancelJob();
};
#endif //PLUSH_PATTERN_GENERATOR_HH