#ifndef PLUSHPLUGIN_HH
#define PLUSHPLUGIN_HH

#include <OpenFlipper/BasePlugin/BaseInterface.hh>
#include <OpenFlipper/BasePlugin/ToolboxInterface.hh>
#include <OpenFlipper/BasePlugin/LoggingInterface.hh>
#include <OpenFlipper/BasePlugin/LoadSaveInterface.hh>
#include <OpenFlipper/BasePlugin/KeyInterface.hh>
#include <OpenFlipper/BasePlugin/PluginFunctions.hh>
#include <OpenFlipper/BasePlugin/RPCInterface.hh>
#include <OpenFlipper/BasePlugin/ProcessInterface.hh>

#include <OpenFlipper/common/Types.hh>

#include <ObjectTypes/TriangleMesh/TriangleMesh.hh>

#include "CGAL_Polyhedron_builder.hh"

#include "SuperDeform/Skeleton.hh"

typedef TriMesh::VertexHandle VertexHandle;
typedef TriMesh::EdgeHandle EdgeHandle;
typedef TriMesh::HalfedgeHandle HalfedgeHandle;
typedef TriMesh::FaceHandle FaceHandle;
typedef TriMesh::VertexIter VertexIter;
typedef TriMesh::EdgeIter EdgeIter;
typedef TriMesh::HalfedgeIter HalfedgeIter;
typedef TriMesh::FaceIter FaceIter;

class PlushPlugin : public QObject,
        BaseInterface,
        ToolboxInterface,
        LoggingInterface,
        LoadSaveInterface,
        KeyInterface,
        RPCInterface,
        ProcessInterface
{
    Q_OBJECT
    Q_INTERFACES(BaseInterface)
    Q_INTERFACES(ToolboxInterface)
    Q_INTERFACES(LoggingInterface)
    Q_INTERFACES(LoadSaveInterface)
    Q_INTERFACES(KeyInterface)
    Q_INTERFACES(RPCInterface)
    Q_INTERFACES(ProcessInterface)
    
signals:
    //BaseInterface
    void updateView();
    void updatedObject(int _identifier, const UpdateType& _type);
    //LoggingInterface
    void log(Logtype _type, QString _message);
    void log(QString _message);
    // ToolboxInterface
    void addToolbox( QString _name, QWidget* _widget );
    // LoadSaveInterface
    void load( QString _filename, DataType _dataType, int &_id);
    void addEmptyObject( DataType _type, int& _id );
    // PickingInterface
    void addPickMode( const std::string& _mode );
    void addHiddenPickMode(const std::string& _mode);
    // KeyInterface
    void registerKey(int _key, Qt::KeyboardModifiers _modifiers, QString _description, bool _multiUse = false);
    // RPCInterface
    void pluginExists( QString _pluginName, bool &_exists );
    void functionExists ( QString _pluginName, QString _functionName, bool &_exists );
    // ProcessInterface
    void startJob (QString _jobId, QString _description, int _min, int _max, bool _blocking=false);
    void finishJob (QString _jobId);
    void setJobState (QString _jobId, int _value);
    void setJobDescription (QString _jobId, QString _text);
    
public:
    PlushPlugin();
    ~PlushPlugin();
    // BaseInterface
    QString name() { return (QString("Simple Plush")); };
    QString description( ) { return (QString("Smooths the active Mesh")); };

    // Curvature
    static OpenMesh::VPropHandleT<double> minCurvatureHandle;
    static OpenMesh::VPropHandleT<double> maxCurvatureHandle;
    static OpenMesh::VPropHandleT<OpenMesh::Vec3d> minCurvatureDirectionHandle;
    static OpenMesh::VPropHandleT<OpenMesh::Vec3d> maxCurvatureDirectionHandle;
    
    // Geodesic
    static OpenMesh::EPropHandleT<double> edgeWeightHandle;
    static OpenMesh::MPropHandleT< std::map<std::pair<VertexHandle, VertexHandle>, double> > geodesicDistanceHandle;
    static OpenMesh::MPropHandleT< std::map<std::pair<VertexHandle, VertexHandle>, IdList> > geodesicPathHandle;
    
    // Skeleton
    static OpenMesh::MPropHandleT<Skeleton*> skeletonHandle;
    // The weight of corresponding bones for each vertex
    static OpenMesh::VPropHandleT<double*> bonesWeightHandle;

private:
    QSpinBox *geodesicEdges;
    QPushButton *geodesicButton;
    QPushButton *geodesicAllButton;
    QPushButton *loadSelectionButton;
    QPushButton *saveSelectionButton;
    QPushButton *clearSelectionButton;
    QPushButton *calcCurvatureButton;
    
    bool showAllPath;

    bool isJobCanceled;
    OpenFlipperThread *thread;
    
    std::vector<char*> *requiredPlugins;
    
    bool isIntersected(TriMesh *mesh, IdList path1, IdList path2);

    bool getEdge(TriMesh *mesh, EdgeHandle &_eh, int v1No, int v2No);
    bool getEdge(TriMesh *mesh, EdgeHandle &_eh, VertexHandle v1, VertexHandle v2);
    
    void loadCurvature(TriMesh *mesh, QString meshName);
    
    int loadSelection(int meshId, QString meshName);
    void saveSelection(int meshId, QString meshName);
    void clearSelection(int meshId);
    
    bool loadSkeleton(int meshId);
    bool loadBoneWeight(int meshId);
    void saveBoneWeight(int meshId);
    
    bool calcCurvature(QString _jobId, int meshId);
    void calcGeodesic(TriMesh *mesh,
                      Polyhedron &P,
                      std::map<int, Polyhedron::Vertex_handle> &verticesMapping,
                      VertexHandle sourceHandle,
                      IdList targetVertices);
    bool calcSpanningTree(QString _jobId, int meshId, std::vector<EdgeHandle> &spanningTree, IdList selectedVertices, int limitNum, bool allPaths);
    
    void initProperties(TriMesh *mesh);
    void uninitProperties(TriMesh *mesh);
    
private slots:
    // BaseInterface
    void initializePlugin();
    void pluginsInitialized();
    
    void calcSkeletonWeightButtonClicked();
    void saveSkeletonWeightButtonClicked();
    void showGeodesicButtonClicked();
    void loadSelectionButtonClicked();
    void saveSelectionButtonClicked();
    void clearSelectionButtonClicked();
    void calcCurvatureButtonClicked();

    // Starter function for thread
    void showGeodesicThread(QString _jobId);
    void calcCurvatureThread(QString _jobId);

    void slotKeyEvent( QKeyEvent* _event );
    
    // LoadSaveInterface
    void fileOpened(int _id);
    void objectDeleted(int _id);
    
    // ProcessInterface
    void canceledJob(QString _job);
    void finishedJob(QString _job);

public slots:
    QString version() { return QString("1.0"); };
};
#endif //PLUSHPLUGIN_HH
