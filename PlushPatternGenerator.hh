#ifndef PLUSH_PATTERN_GENERATOR_HH
#define PLUSH_PATTERN_GENERATOR_HH

// This one must be included before any other OpenFlipper headers, otherwise it will cause weird runtime errors.
#include <OpenFlipper/BasePlugin/BaseInterface.hh>

#include "SuperDeform/Skeleton.hh"

#include "FilteredTriMesh.hh"
#include "Common.hh"

#include <Eigen/sparse>

#include <boost/graph/adjacency_list.hpp>

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
    typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
    boost::property<boost::vertex_owner_t, FilteredTriMesh>,
    boost::property<boost::edge_index_t, std::size_t,
    boost::property<boost::edge_owner_t, std::vector<HalfedgeHandle>>>
    > SubMesh_graph;

    /// Used for log
    ///@{
    static const int LOGERR = 0;
    static const int LOGINFO = 1;
    ///@}

    /// @name Merging property handle
    ///@{
    static OpenMesh::MPropHandleT< std::map< std::set<HalfedgeHandle>, double> > joint_boundary_area_handle;
    static OpenMesh::MPropHandleT< std::map< std::set<HalfedgeHandle>, double> > joint_boundary_distortion_handle;

    static OpenMesh::MPropHandleT< std::set<VertexHandle> > intersection_points_handle;

    static OpenMesh::EPropHandleT<int> segment_no_handle;
    
    static OpenMesh::FPropHandleT<int> face_to_patch_idx_handle;

    static OpenMesh::MPropHandleT< std::set<EdgeHandle> > seams_handle;
    
    static OpenMesh::MPropHandleT< std::vector< std::vector<FilteredTriMesh> > > hierarchical_patches_handle;
    ///@}


    /// @name Curvature property handle
    ///@{
    /// Curvature value/direction of each vertex. It's empty before curvature calculation.
    static OpenMesh::VPropHandleT<double> meanCurvatureHandle;
    static OpenMesh::VPropHandleT<double> gaussianCurvatureHandle;
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
    static OpenMesh::MPropHandleT<SubMesh_graph> flattenedSubMeshesHandle;
    /// Distortion indicator for visualization
    static OpenMesh::FPropHandleT<double> distortionFHandle;
    static OpenMesh::VPropHandleT<double> distortionVHandle;
    static OpenMesh::VPropHandleT<TriMesh::Point> flattenedPositionHandle;
    ///@}
    
    /// @name Handle of new sub-mesh
    ///@{
    /// Mapping flattened vertex handle back to original vertex handle
    static OpenMesh::VPropHandleT<VertexHandle> getInverseMappingHandle(TriMesh *mesh);
    ///@}
    
    static bool getHalfedge(const FilteredTriMesh &mesh, HalfedgeHandle &heh, VertexHandle from, VertexHandle to);
    static bool getHalfedge(const TriMesh *mesh, HalfedgeHandle &heh, int fromNo, int toNo);
    static bool getHalfedge(const TriMesh *mesh, HalfedgeHandle &heh, VertexHandle from, VertexHandle to);
    static bool getEdge(const TriMesh *mesh, EdgeHandle &eh, int v1No, int v2No);
    static bool getEdge(const TriMesh *mesh, EdgeHandle &eh, VertexHandle v1, VertexHandle v2);
    static bool getFace(const TriMesh *mesh, FaceHandle &fh, VertexHandle v1, VertexHandle v2, VertexHandle v3);
    
    /// Calculate the sum of inner (clockwise) angles by iterating from one halfedge to another.
    static double getSumInnerAngle(const TriMesh *mesh, VertexHandle v);
    
    /// Get boundary for a given opened mesh
    static std::vector< std::vector<HalfedgeHandle> >  getBoundaryOfOpenedMesh(const FilteredTriMesh &mesh,
                                        bool getInteriorHalfedge);
    
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

    /// Seams
    ///@{
    bool calcSeams(std::vector<VertexHandle> selectedVertices,
                   double developable_threshold = 0.1,
                   int limitNum = 0,
                   bool elimination = false,
                   bool allPaths = true);
    bool calcLocalSeams(FilteredTriMesh &mesh, double developable_threshold);
    bool set_seams_from_selection();
    bool load_seams();
    bool save_seams();
    
    std::set<EdgeHandle> get_seams();
    ///@}

    /// Skeleton
    ///@{
    bool calcSkeletonWeight();
    bool loadSkeleton();
    bool loadBoneWeight();
    bool saveBoneWeight();
    ///@}

    /// Curvature
    ///@{
    bool calcCurvature();
    void loadCurvature();
    ///@}
    
    bool calcSelection(std::vector<VertexHandle> &targetVertices, double threshold=0.001);
    
    /// Geodesic
    ///@{
    void calcGeodesic(std::vector<VertexHandle> targetVertices);
    bool loadGeodesic();
    bool saveGeodesic(std::vector<VertexHandle> selectedVertices);

    void set_geodesic_coeffifients(double distanceCoefficient,
                                   double textureCoefficient,
                                   double curvatureCoefficient,
                                   double skeletonCoefficient,
                                   double pathCoefficient);
    ///@}
    
    /// Merging
    ///@{
    void construct_subsets(double threshold);
    FilteredTriMesh merge_patch(FilteredTriMesh &patch1,
                                FilteredTriMesh &patch2,
                                int seam_segment_idx,
                                std::set<EdgeHandle> &seam_segment);
    void optimize_patches(double threshold, bool step);
    
    std::vector< std::vector<FilteredTriMesh> > get_hierarchical_patches();
    ///@}

    // Used in selection and flattening
    Eigen::SparseMatrix<double> calc_parameterization_weight_matrix(FilteredTriMesh &mesh,
                                                                    InteriorParameterizationMethod method);
    void calc_parameterization_weight_matrix(TriMesh *mesh,
                                             Eigen::SparseMatrix<double> &M,
                                             InteriorParameterizationMethod method);
    
    /// Flattening
    ///@{
    /// Flatten 3D loops into 2D loops using LPFB.
    bool calcFlattenedGraph();
    std::vector<FilteredTriMesh> getFlattenedSubMeshes();
    bool get_triMesh_from_subMesh(TriMesh *result_triMesh, FilteredTriMesh &subMesh, bool use_flattened_position);
    ///@}
    
    VertexHandle get_original_handle(TriMesh *mesh, const VertexHandle vh) const;
    EdgeHandle get_original_handle(TriMesh *mesh, const EdgeHandle eh) const;
    HalfedgeHandle get_original_handle(TriMesh *mesh, const HalfedgeHandle eh) const;
    FaceHandle get_original_handle(TriMesh *mesh, const FaceHandle fh) const;
    
    /// Visualization utils
    ///@{
    std::set<VertexHandle> get_intersection_points();
    std::vector< std::vector<HalfedgeHandle> > get_seams_segments();
    ///@}
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
    
    std::vector< std::vector<EdgeHandle> > get_segments_from_seams(std::set<EdgeHandle> &seams);
    std::vector< std::vector<HalfedgeHandle> > get_halfedge_segments_from_seams(std::set<EdgeHandle> &seams);
    std::vector< std::vector<HalfedgeHandle> > get_halfedge_segments_from_seams(std::set<VertexHandle> intersection_points,
                                                                                std::set<EdgeHandle> &seams);
    std::set<VertexHandle> get_intersection_points(std::set<EdgeHandle> &seams);
    ///@}
    
    /// @name Sub-mesh
    ///@{
    SubMesh_graph get_subMeshes_with_boundary(std::set<EdgeHandle> &seams);
    FilteredTriMesh get_subMesh_with_boundary(FaceHandle root_face, std::set<HalfedgeHandle> &seams);
    FilteredTriMesh get_subMesh_with_boundary(FaceHandle root_face, std::set<EdgeHandle> &seams);
    ///@}

    /// @name Boundary
    ///@{
    void get_closed_boundaries_of_seams(std::vector< std::vector<HalfedgeHandle> > &closed_seams, std::set<EdgeHandle> &seams);
    
    bool get_adjacent_boundary(HalfedgeHandle start_heh,
                            std::set<EdgeHandle> &seams,
                            std::set<HalfedgeHandle> &adj_boundary);
    bool get_joint_boundary(std::set<HalfedgeHandle> &adj_boundary1,
                            std::set<HalfedgeHandle> &adj_boundary2,
                            std::vector<HalfedgeHandle> &joint_seam_segment,
                            std::set<HalfedgeHandle> &joint_boundary,
                            std::set<EdgeHandle> &seams);
    bool get_joint_submesh(std::vector<HalfedgeHandle> &joint_seam_segment,
                           std::set<EdgeHandle> &seams);
    
    ///@}
    /// @name Flattening
    ///@{
    bool calcLPFB(FilteredTriMesh &mesh, std::map<HalfedgeHandle, OpenMesh::Vec3d> &boundaryPosition);
    bool calcInteriorPoints(FilteredTriMesh &mesh, std::map<HalfedgeHandle, OpenMesh::Vec3d> &boundaryPosition);
    /// Calculate the distortion between original 3D mesh and flattened graph, the result is stored in distortionHandle
    void calcDistortion(std::vector<FilteredTriMesh> &flattenedMeshes);
    void calcDistortion(FilteredTriMesh &flattenedMeshes);
    ///@}
    
    void initProperties();
    void uninitProperties();
    
    std::set<HalfedgeHandle> prevBoundary;
    std::vector<HalfedgeHandle> prevSegment;
    
    public slots:
    /// Cancel time-consuming works. (calcGeodesic, calcCurvature, etc.)
    void cancelJob();
};
#endif //PLUSH_PATTERN_GENERATOR_HH
