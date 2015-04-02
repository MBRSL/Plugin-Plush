#ifndef PLUSH_PATTERN_GENERATOR_HH
#define PLUSH_PATTERN_GENERATOR_HH

// This one must be included before any other OpenFlipper headers, otherwise it will cause weird runtime errors.
#include <OpenFlipper/BasePlugin/BaseInterface.hh>
#include <ObjectTypes/PolyLine/PolyLine.hh>

#include "CGAL_Polyhedron_builder.hh"

#include "SuperDeform/Skeleton.hh"

#include "Common.hh"

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
    TriMesh *m_mesh;
    QString m_meshName;
    
    /// Used for log
    ///@{
    static const int LOGERR = 0;
    static const int LOGINFO = 1;
    ///@}

    /// @name Curvature property handle
    ///@{
    /// Curvature value/direction of each vertex. It's empty before curvature calculation.
    static OpenMesh::VPropHandleT<double> minCurvatureHandle;
    static OpenMesh::VPropHandleT<double> maxCurvatureHandle;
    static OpenMesh::VPropHandleT<OpenMesh::Vec3d> minCurvatureDirectionHandle;
    static OpenMesh::VPropHandleT<OpenMesh::Vec3d> maxCurvatureDirectionHandle;
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
    static OpenMesh::MPropHandleT< std::vector< std::vector<HalfedgeHandle> > > boundariesHandle;
    /// This stores flattened meshes.
    static OpenMesh::MPropHandleT< std::vector<TriMesh> > flattenedMeshesHandle;
    /// Distortion indicator for visualization
    static OpenMesh::FPropHandleT<double> distortionFHandle;
    static OpenMesh::VPropHandleT<double> distortionVHandle;
    /// Mapping flattened vertex handle back to original vertex handle
    static OpenMesh::VPropHandleT<VertexHandle> getInverseMappingHandle(TriMesh *mesh);
    static OpenMesh::MPropHandleT< std::vector<EdgeHandle> > getSeamsHandle(TriMesh *mesh);
    ///@}
    
    static bool getHalfedge(TriMesh *mesh, HalfedgeHandle &heh, int fromNo, int toNo);
    static bool getHalfedge(TriMesh *mesh, HalfedgeHandle &heh, VertexHandle from, VertexHandle to);
    static bool getEdge(TriMesh *mesh, EdgeHandle &eh, int v1No, int v2No);
    static bool getEdge(TriMesh *mesh, EdgeHandle &eh, VertexHandle v1, VertexHandle v2);
    static bool getFace(TriMesh *mesh, FaceHandle &fh, VertexHandle v1, VertexHandle v2, VertexHandle v3);
    static OpenMesh::Vec3d getVector(TriMesh *mesh, EdgeHandle &_eh);
    static OpenMesh::Vec3d getVector(TriMesh *mesh, HalfedgeHandle &_heh);
    
    /// Calculate the sum of inner (clockwise) angles by iterating from one halfedge to another.
    static double getSumInnerAngle(TriMesh *mesh, HalfedgeHandle heh1, HalfedgeHandle heh2);
    
    /// Get boundary for a given opened mesh
    static bool getBoundaryOfOpenedMesh(std::vector<HalfedgeHandle> &boundary, TriMesh &mesh);
    
    PlushPatternGenerator(TriMesh *mesh, QString meshName);
    ~PlushPatternGenerator();
    
    void loadCurvature();
    
    bool loadSkeleton();
    bool loadBoneWeight();
    bool saveBoneWeight();
    
    bool calcCurvature();
    
    void calcGeodesic(std::vector<VertexHandle> targetVertices);
    bool loadGeodesic();
    bool saveGeodesic(std::vector<VertexHandle> selectedVertices);
    
    bool calcSeams(std::vector<VertexHandle> selectedVertices, int limitNum = 0, bool elimination = false, bool allPaths = true);

    bool calcSkeletonWeight();
    
    /// Flatten 3D loops into 2D loops using LPFB.
    bool calcFlattenedGraph();
    
    std::vector<TriMesh>* getFlattenedMeshes();
    std::vector<EdgeHandle>* getSeams();
    
private:
    Polyhedron m_polyhedron;
    
    /// Mapping from TriMesh::VertexHandle::idx to boost_vertex_descriptor
    std::map<int, boost_vertex_descriptor> m_verticesMapping;
    
    bool isJobCanceled;
    
    bool isIntersected(std::vector<VertexHandle> path1, std::vector<VertexHandle> path2);
    
    bool splitWithBoundary(std::vector<TriMesh> *subMeshes, std::vector< std::vector<HalfedgeHandle> > *loops);

    void getBoundariesByEdges(std::vector< std::vector<HalfedgeHandle> > *boundaries, std::vector<EdgeHandle> *separator);

    /// @name Flattening
    ///@{
    bool calcLPFB(TriMesh *mesh);
    bool calcInteriorPoints(std::vector<TriMesh> *flattenedMeshes);
    /// Calculate the distortion between original 3D mesh and flattened graph, the result is stored in distortionHandle
    void calcDistortion();
    /// Organize 2D loops to prevent overlapping while minimize bounding area.
    bool packFlattenedGraph(std::vector<TriMesh> *flattenedMeshes, const int nColumns=4);
    ///@}
    
    void initProperties();
    void uninitProperties();
    
    public slots:
    /// Cancel time-consuming works. (calcGeodesic, calcCurvature, etc.)
    void cancelJob();
};
#endif //PLUSH_PATTERN_GENERATOR_HH
