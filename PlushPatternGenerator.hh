#ifndef PLUSH_PATTERN_GENERATOR_HH
#define PLUSH_PATTERN_GENERATOR_HH

// This one must be included before any other OpenFlipper headers, otherwise it will cause weird runtime errors.
#include <OpenFlipper/BasePlugin/BaseInterface.hh>

#include "SuperDeform/Skeleton.hh"

#include "Common.hh"

#include <Eigen/sparse>

#include <ObjectTypes/PolyLine/PolyLine.hh>
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
    void log(QString msg);
    
    /// Use this signal to inform user current job progress.
    void setJobState(int val);
    
    void updateView();
    
public:
    /// Used for log
    ///@{
    static const int LOGERR = 0;
    static const int LOGINFO = 1;
    ///@}

    /// @name Merging property handle
    ///@{
    static OpenMesh::MPropHandleT<int> merge_iterations_handle;

    static OpenMesh::MPropHandleT< std::map< std::vector<HalfedgeHandle>, double> > joint_boundary_area_handle;
    static OpenMesh::MPropHandleT< std::map< std::vector<HalfedgeHandle>, double> > joint_boundary_distortion_handle;
    
    static OpenMesh::EPropHandleT<int> segment_no_handle;
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
    
    static bool is_different_texture(TriMesh *mesh, EdgeHandle eh);
    static bool is_different_texture(TriMesh *mesh, HalfedgeHandle eh);
    
    enum InteriorParameterizationMethod { Barycentric, Conformal, MeanValue };
    
    PlushPatternGenerator(TriMesh *mesh, QString meshName);
    ~PlushPatternGenerator();
    
    bool load_seams();
    bool save_seams();

    void loadCurvature();
    
    bool loadSkeleton();
    bool loadBoneWeight();
    bool saveBoneWeight();
    
    bool calcCurvature();
    
    bool calcSelection(std::vector<VertexHandle> &targetVertices, double threshold=0.001);
    
    void calcGeodesic(std::vector<VertexHandle> targetVertices);
    bool loadGeodesic();
    bool saveGeodesic(std::vector<VertexHandle> selectedVertices);
    void optimize_patches(double threshold, bool step);

    bool calcSeams(std::vector<VertexHandle> selectedVertices,
                   double developable_threshold = 0.1,
                   int limitNum = 0,
                   bool elimination = false,
                   bool allPaths = true);
    bool calcSkeletonWeight();
    
    // Used in selection
    void calc_parameterization_weight_matrix(TriMesh *mesh, Eigen::SparseMatrix<double> &M, InteriorParameterizationMethod method);
    
    /// Flatten 3D loops into 2D loops using LPFB.
    bool calcFlattenedGraph();
    
    void set_geodesic_coeffifients(double distanceCoefficient,
                                   double textureCoefficient,
                                   double curvatureCoefficient,
                                   double skeletonCoefficient,
                                   double pathCoefficient);
    
    std::vector<TriMesh>* getFlattenedMeshes();
    std::set<EdgeHandle>* getSeams();
    
    VertexHandle get_original_handle(TriMesh *mesh, const VertexHandle vh) const;
    EdgeHandle get_original_handle(TriMesh *mesh, const EdgeHandle eh) const;
    HalfedgeHandle get_original_handle(TriMesh *mesh, const HalfedgeHandle eh) const;
    FaceHandle get_original_handle(TriMesh *mesh, const FaceHandle fh) const;
    
    void get_intersection_points(std::set<EdgeHandle> *seams,
                                 std::set<VertexHandle> &intersection_points);
    void get_segments_from_seams(std::vector< std::vector<EdgeHandle> > &segments);
    void get_segments_from_seams(std::vector< std::vector<HalfedgeHandle> > &segments);
    void get_segments_from_seams(std::vector< std::vector<HalfedgeHandle> > &segments,
                                 std::set<VertexHandle> intersection_points);

    bool calcLocalSeams(TriMesh *mesh, double developable_threshold);
    
    void show_intersection_points();

private:
    TriMesh *m_mesh;
    QString m_meshName;
    
    bool isJobCanceled;
    
    /// @name Geodesic coefficients
    ///@{
    double m_distanceCoefficient = 0.2;
    double m_textureCoefficient = 1;
    double m_curvatureCoefficient = 0.4;
    double m_skeletonCoefficient = 0.7;
    double m_pathCoefficient = 0.2;
    ///@}
    
    /// @name Seams
    ///@{
    bool isIntersected(std::vector<VertexHandle> path1, std::vector<VertexHandle> path2);
    bool calcCircularSeams(TriMesh *mesh, double threshold, bool is_local_seam);
    bool calcStructralSeams(TriMesh *mesh,
                            std::vector<VertexHandle> selectedVertices,
                            int limitNum,
                            bool elimination,
                            bool allPaths);
    ///@}
    
    /// @name Sub mesh
    ///@{
    bool splitWithBoundary(std::vector<TriMesh> *subMeshes, std::set<EdgeHandle> *seams);
    bool extract_mesh_with_boundary(TriMesh *new_mesh, FaceHandle root_face, std::set<EdgeHandle> *seams);
    bool extract_mesh_with_boundary(TriMesh *new_mesh, FaceHandle root_face, std::vector<HalfedgeHandle> *seams);
    ///@}

    /// @name Boundary
    ///@{
    void get_closed_boundaries_of_seams(std::vector< std::vector<HalfedgeHandle> > *closed_seams, std::set<EdgeHandle> *seams);
    
    bool get_adjacent_boundary(HalfedgeHandle start_heh,
                            std::set<EdgeHandle> &seams,
                            std::vector<HalfedgeHandle> &adj_boundary);
    bool get_joint_boundary(std::vector<HalfedgeHandle> &adj_boundary1,
                            std::vector<HalfedgeHandle> &adj_boundary2,
                            std::set<EdgeHandle> &seams,
                            std::vector<HalfedgeHandle> &joint_boundary);
    ///@}
    /// @name Flattening
    ///@{
    bool calcLPFB(TriMesh *mesh, std::map<VertexHandle, OpenMesh::Vec3d> *boundaryPosition);
    bool calcInteriorPoints(TriMesh *mesh, std::map<VertexHandle, OpenMesh::Vec3d> *boundaryPosition);
    /// Calculate the distortion between original 3D mesh and flattened graph, the result is stored in distortionHandle
    void calcDistortion(std::vector<TriMesh> *flattenedMeshes);
    void calcDistortion(TriMesh &flattenedMeshes);
    /// Organize 2D loops to prevent overlapping while minimize bounding area.
    bool packFlattenedGraph(std::vector<TriMesh> *flattenedMeshes, const int nColumns=4);
    ///@}
    
    void initProperties();
    void uninitProperties();
    
    std::vector<HalfedgeHandle> prevBoundary;
    std::vector<HalfedgeHandle> prevSegment;
    
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
            bool found = getFace(m_mesh, originalF, originalV[0], originalV[1], originalV[2]);
            assert(found);
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
            bool found = getEdge(m_mesh, original_eh, originalV1, originalV2);
            assert(found);
            edgesId.push_back(original_eh.idx());
        }
        MeshSelection::selectEdges(m_mesh, edgesId);
    }
    
    public slots:
    /// Cancel time-consuming works. (calcGeodesic, calcCurvature, etc.)
    void cancelJob();
};
#endif //PLUSH_PATTERN_GENERATOR_HH
