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

#include "GeodesicDistance/geodesic_mesh.hh"


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

private:
    bool isJobCanceled;
    OpenFlipperThread *thread;
    
    std::vector<char*> *requiredPlugins;

    std::map<std::pair<VertexHandle, VertexHandle>, double> *geodesicDistance;
    std::map<std::pair<VertexHandle, VertexHandle>, IdList> *geodesicPath;
    
    //std::map<VertexHandle, double> *curvatureK2;
    OpenMesh::VPropHandleT<double> minCurvatureHandle;
    OpenMesh::VPropHandleT<double> maxCurvatureHandle;
    OpenMesh::VPropHandleT<OpenMesh::Vec3d> minCurvatureDirectionHandle;
    OpenMesh::VPropHandleT<OpenMesh::Vec3d> maxCurvatureDirectionHandle;
    
    bool getEdge(TriMesh *mesh, EdgeHandle &_eh, int v1No, int v2No);
    bool getEdge(TriMesh *mesh, EdgeHandle &_eh, VertexHandle v1, VertexHandle v2);
    double getEdgeWeight(TriMesh *mesh, int v1No, int v2No);
    double getEdgeWeight(TriMesh *mesh, VertexHandle v1, VertexHandle v2);
    double getEdgeWeight(TriMesh *mesh, EdgeHandle eh);
    
    void calcGeodesic(TriMesh *mesh, VertexHandle sourceHandle);
//    void findPath(TriMesh *mesh, std::set<EdgeHandle> &spanningTree, std::vector<VertexHandle> &path, VertexHandle sourceHandle, VertexHandle destHandle);
    void calcSpanningTree(TriMesh *mesh, int meshId, std::set<EdgeHandle> &result, IdList selectedVertices);
    void getOrderedSelectedVertices(TriMesh *mesh, int meshId, IdList *selectedVertices);
    
private slots:
    // BaseInterface
    void initializePlugin();
    void pluginsInitialized();

    void showGeodesic();
    void showRidge();
    
    void loadCurvatureButtonClicked();
    void calcCurvatureButtonClicked();
    void showCurvatureButtonClicked();
    void calcCurvature(QString _jobId);

    void slotKeyEvent( QKeyEvent* _event );
    
    // ProcessInterface
    void canceledJob(QString _job);
    void finishedJob(QString _job);

public slots:
    QString version() { return QString("1.0"); };
};
#endif //PLUSHPLUGIN_HH