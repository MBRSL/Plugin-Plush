#ifndef PLUSHPLUGIN_HH
#define PLUSHPLUGIN_HH

#include <OpenFlipper/BasePlugin/BaseInterface.hh>
#include <OpenFlipper/BasePlugin/ToolboxInterface.hh>
#include <OpenFlipper/BasePlugin/LoggingInterface.hh>
#include <OpenFlipper/BasePlugin/LoadSaveInterface.hh>
#include <OpenFlipper/BasePlugin/KeyInterface.hh>
#include <OpenFlipper/BasePlugin/PluginFunctions.hh>
#include <OpenFlipper/BasePlugin/ProcessInterface.hh>
#include <OpenFlipper/BasePlugin/RPCInterface.hh>

#include "PlushPatternGenerator.hh"

#include "Common.hh"

/**
 * @brief The main class for this plugin.
 *
 * Including basic initializePlugin(), pluginsInitialized(), and other UI-related functions.
 */
class PlushPlugin : public QObject,
        BaseInterface,
        ToolboxInterface,
        LoggingInterface,
        LoadSaveInterface,
        KeyInterface,
        ProcessInterface,
        RPCInterface
{
    Q_OBJECT
    Q_INTERFACES(BaseInterface)
    Q_INTERFACES(ToolboxInterface)
    Q_INTERFACES(LoggingInterface)
    Q_INTERFACES(LoadSaveInterface)
    Q_INTERFACES(KeyInterface)
    Q_INTERFACES(ProcessInterface)
    Q_INTERFACES(RPCInterface)
    
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
    void deleteObject( int _id );
    // PickingInterface
    void addPickMode( const std::string& _mode );
    void addHiddenPickMode(const std::string& _mode);
    // KeyInterface
    void registerKey(int _key, Qt::KeyboardModifiers _modifiers, QString _description, bool _multiUse = false);
    // ProcessInterface
    void startJob (QString _jobId, QString _description, int _min, int _max, bool _blocking=false);
    void finishJob (QString _jobId);
    void setJobState (QString _jobId, int _value);
    void setJobDescription (QString _jobId, QString _text);
    
    // Custom canceling signal
    void cancelingJob();
    
public:
    PlushPlugin();
    ~PlushPlugin();
    // BaseInterface
    QString name() { return (QString("Simple Plush")); };
    QString description( ) { return (QString("Smooths the active Mesh")); };

private:
    PlushPatternGenerator *m_patternGenerator;
    /// This obj contains the mesh used by pattern generator. It self contains info such as id and filename in OpenFlipper
    TriMeshObject *m_triMeshObj;
    
    // seams
    QCheckBox *seamElimination;
    QSpinBox *seamNumPaths;
    QPushButton *seamShowSingleButton;
    QPushButton *seamShowAllButton;
    QPushButton *seamLocalButton;
    QPushButton *seamMergeButton;
    QCheckBox *seamMergeIsStep;
    QDoubleSpinBox *seamMergeThreshold;
    
    // geodesic
    QDoubleSpinBox *geodesic_distance_coef;
    QDoubleSpinBox *geodesic_skeleton_coef;
    QDoubleSpinBox *geodesic_curvature_coef;
    QDoubleSpinBox *geodesic_smoothness_coef;
    
    // subset
    QDoubleSpinBox *subset_merge_threshold;

    bool showAllPath;
    /// Flag for thread stopping
    bool isJobStopped;
    QString m_currentJobId;
    
    std::vector<char*> *requiredPlugins;
    
    int loadSelection(TriMesh *mesh, QString meshName);
    void saveSelection(TriMesh *mesh, QString meshName);
    void clearSelection(TriMesh *mesh);
    void calcSelection(TriMesh *mesh);

private slots:
    // BaseInterface
    void initializePlugin();
    void pluginsInitialized();
    
    void seamShowButtonClicked();
    void seamSaveButtonClicked();
    void seamLoadButtonClicked();
    void seamMergeButtonClicked();
    
    void subset_calc_button_clicked();
    void subset_show_button_clicked();
    
    void calcGeodesicButtonClicked();
    
    void loadSelectionButtonClicked();
    void saveSelectionButtonClicked();
    void clearSelectionButtonClicked();
    void calcSelectionButtonClicked();
    
    void showFlattenedGrpahButtonClicked();
    void calcFlattenedGraphButtonClicked();
    
    void calcCurvatureButtonClicked();
    
    void calcSkeletonWeightButtonClicked();
    
    void vis_intersection_points_button_clicked();
    void vis_seam_segments_button_clicked();
    void vis_skeleton_button_clicked();

    // Starter function for thread
    void calcGeodesicThread();
    void calcFlattenedGraphThread();
    void calcMergeSegmentThread();
    void calcSubsetThread();
    
    void slotKeyEvent( QKeyEvent* _event );
    
    // LoadSaveInterface
    void fileOpened(int _id);
    void objectDeleted(int _id);
    
    // ProcessInterface
    void canceledJob(QString _job);
    // Custom handler
    void finishedJobHandler();
    
    bool checkIfGeneratorExist();
    
    /// Receive job state from m_patternGenerator and then re-emit with jobId
    void receiveJobState(int);
    
    /// Receive message from m_patternGenerator and then re-emit
    void receiveLog(int, QString);
    
    void receiveUpdate();

public slots:
    QString version() { return QString("1.0"); };
};
#endif //PLUSHPLUGIN_HH
