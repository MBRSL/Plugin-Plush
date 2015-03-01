#ifndef PLUSH_PATTERN_GENERATOR_HH
#define PLUSH_PATTERN_GENERATOR_HH

// This one must be included before any other OpenFlipper headers, otherwise it will cause weird runtime errors.
#include <OpenFlipper/BasePlugin/BaseInterface.hh>

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
    static OpenMesh::MPropHandleT< std::map<std::pair<VertexHandle, VertexHandle>, std::vector<int> > > geodesicPathHandle;
    ///@}
    
    /// @name Skeleton property handle
    ///@{
    
    /// The skeleton of this mesh. It's built during mesh loading unless error occurs
    static OpenMesh::MPropHandleT<Skeleton*> skeletonHandle;
    /// The weight of corresponding bones for each vertex
    static OpenMesh::VPropHandleT<double*> bonesWeightHandle;
    ///@}
    
    PlushPatternGenerator(TriMesh *mesh, QString meshName);
    ~PlushPatternGenerator();
    
    void loadCurvature();
    
    bool loadSkeleton();
    bool loadBoneWeight();
    bool saveBoneWeight();
    
    bool calcCurvature();
    void calcGeodesic(VertexHandle sourceHandle, std::vector<int> targetVertices);
    bool calcSpanningTree(std::vector<EdgeHandle> &spanningTree, std::vector<int> selectedVertices, int limitNum, bool allPaths);
    bool calcSkeletonWeight();
    
private:
    Polyhedron m_polyhedron;
    
    /// Mapping from TriMesh::VertexHandle::idx to boost_vertex_descriptor
    std::map<int, boost_vertex_descriptor> m_verticesMapping;
    
    static int LOGERR;
    static int LOGINFO;
    
    bool isJobCanceled;
    
    bool isIntersected(std::vector<int> path1, std::vector<int> path2);
    
    static bool getEdge(TriMesh *mesh, EdgeHandle &_eh, int v1No, int v2No);
    static bool getEdge(TriMesh *mesh, EdgeHandle &_eh, VertexHandle v1, VertexHandle v2);
    
    void initProperties();
    void uninitProperties();
    
    public slots:
    /// Cancel time-consuming works. (calcGeodesic, calcCurvature, etc.)
    void cancelJob();
};
#endif //PLUSH_PATTERN_GENERATOR_HH
