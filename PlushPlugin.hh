#ifndef PLUSHPLUGIN_HH
#define PLUSHPLUGIN_HH
#include <OpenFlipper/BasePlugin/BaseInterface.hh>
#include <OpenFlipper/BasePlugin/ToolboxInterface.hh>
#include <OpenFlipper/BasePlugin/LoggingInterface.hh>
#include <OpenFlipper/BasePlugin/LoadSaveInterface.hh>
#include <OpenFlipper/BasePlugin/KeyInterface.hh>
#include <OpenFlipper/common/Types.hh>
#include <ObjectTypes/TriangleMesh/TriangleMesh.hh>
#include "GeodesicDistance/geodesic_mesh.hh"

class PlushPlugin : public QObject,
        BaseInterface,
        ToolboxInterface,
        LoggingInterface,
        LoadSaveInterface,
        KeyInterface
{
    Q_OBJECT
    Q_INTERFACES(BaseInterface)
    Q_INTERFACES(ToolboxInterface)
    Q_INTERFACES(LoggingInterface)
    Q_INTERFACES(LoadSaveInterface)
    Q_INTERFACES(KeyInterface)
    
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

public:
    PlushPlugin();
    // BaseInterface
    QString name() { return (QString("Simple Plush")); };
    QString description( ) { return (QString("Smooths the active Mesh")); };

private:
    std::vector<char*> requiredPlugins;
    std::map<int, geodesic::Mesh*> meshes_geodesic;
    
    bool pickEdge(TriMesh *mesh, OpenMesh::EdgeHandle &_eh, OpenMesh::Vec3d p1, OpenMesh::Vec3d p2);
    bool translate_openMesh_to_geodesic_mesh(TriMesh *mesh, std::vector<double> &points, std::vector<unsigned> &faces);
private slots:
    // BaseInterface
    void initializePlugin();
    void pluginsInitialized();

    void calcGeodesic(PluginFunctions::ObjectIterator o_it);
    void showGeodesic();
    
    void slotKeyEvent( QKeyEvent* _event );

public slots:
    QString version() { return QString("1.0"); };
};
#endif //PLUSHPLUGIN_HH