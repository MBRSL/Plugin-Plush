#include <ObjectTypes/PolyMesh/PolyMesh.hh>
#include <ObjectTypes/TriangleMesh/TriangleMesh.hh>

#include <OpenFlipper/BasePlugin/RPCWrappers.hh>
#include <OpenFlipper/common/Types.hh>

#include "PlushPlugin.hh"
#include "SuperDeform/Weight.hh"

OpenMesh::VPropHandleT<double> PlushPlugin::minCurvatureHandle;
OpenMesh::VPropHandleT<double> PlushPlugin::maxCurvatureHandle;
OpenMesh::VPropHandleT<OpenMesh::Vec3d> PlushPlugin::minCurvatureDirectionHandle;
OpenMesh::VPropHandleT<OpenMesh::Vec3d> PlushPlugin::maxCurvatureDirectionHandle;

OpenMesh::EPropHandleT<double> PlushPlugin::edgeWeightHandle;
OpenMesh::MPropHandleT< std::map<std::pair<VertexHandle, VertexHandle>, double> > PlushPlugin::geodesicDistanceHandle;
OpenMesh::MPropHandleT< std::map<std::pair<VertexHandle, VertexHandle>, IdList> > PlushPlugin::geodesicPathHandle;

OpenMesh::MPropHandleT<Skeleton*> PlushPlugin::skeletonHandle;
OpenMesh::VPropHandleT<double*> PlushPlugin::bonesWeightHandle;

PlushPlugin::PlushPlugin()
{
    requiredPlugins = new std::vector<char*>();
    
    thread = NULL;
}

void PlushPlugin::initializePlugin()
{
    requiredPlugins->push_back("meshobjectselection");
    
    // Register keys
    emit registerKey(Qt::Key_X,     Qt::NoModifier, "Normal mode");
    emit registerKey(Qt::Key_Z,     Qt::NoModifier, "Plush picking mode");
}

PlushPlugin::~PlushPlugin() {
    delete requiredPlugins;
}
void PlushPlugin::pluginsInitialized() {
    bool isPluginsExist;
    for (std::vector<char*>::iterator it = requiredPlugins->begin(); it != requiredPlugins->end(); it++) {
        pluginExists(*it, isPluginsExist);
        if (!isPluginsExist) {
            char str[100];
            sprintf(str, "Plugin %s not found.", *it);
            emit log(LOGERR, str);
            return;
        }
    }
    
    // Create the Toolbox Widget
    QWidget* toolBox = new QWidget();
    QGridLayout* layout = new QGridLayout(toolBox);
    
    // Create button that can be toggled
    // to (de)activate plugin's picking mode
    QGroupBox *skeletonWeightGroup = new QGroupBox(tr("Skeleton weight"));
    QPushButton *calcSkeletonWeightButton = new QPushButton(tr("Calculate"));
    QPushButton *saveSkeletonWeightButton = new QPushButton(tr("Save"));
    QHBoxLayout *skeletonWeightLayout = new QHBoxLayout;
    skeletonWeightLayout->addWidget(calcSkeletonWeightButton);
    skeletonWeightLayout->addWidget(saveSkeletonWeightButton);
    skeletonWeightGroup->setLayout(skeletonWeightLayout);
    
    QGroupBox *geodesicGroup = new QGroupBox(tr("Geodesic paths"));
    QLabel *geodesicNumberLabel = new QLabel(tr("#"));
    geodesicEdges = new QSpinBox();
    geodesicEdges->setMinimum(0);
    geodesicButton = new QPushButton(tr("Show"));
    geodesicAllButton = new QPushButton(tr("Show all"));
    QHBoxLayout *geodesicLayout = new QHBoxLayout;
    geodesicLayout->addWidget(geodesicNumberLabel);
    geodesicLayout->addWidget(geodesicEdges);
    geodesicLayout->addWidget(geodesicButton);
    geodesicLayout->addWidget(geodesicAllButton);
    geodesicGroup->setLayout(geodesicLayout);
    
    QGroupBox *selectionGroup = new QGroupBox(tr("Selection"));
    loadSelectionButton = new QPushButton(tr("Load"));
    saveSelectionButton = new QPushButton(tr("Save"));
    clearSelectionButton = new QPushButton(tr("Clear"));
    QHBoxLayout *selectionLayout = new QHBoxLayout;
    selectionLayout->addWidget(loadSelectionButton);
    selectionLayout->addWidget(saveSelectionButton);
    selectionLayout->addWidget(clearSelectionButton);
    selectionGroup->setLayout(selectionLayout);
    
    QGroupBox *curvatureGroup = new QGroupBox(tr("Curvature"));
    calcCurvatureButton = new QPushButton(tr("Calculate"));
    QHBoxLayout *curvatureLayout = new QHBoxLayout;
    curvatureLayout->addWidget(calcCurvatureButton);
    curvatureGroup->setLayout(curvatureLayout);
    
    layout->addWidget(skeletonWeightGroup);
    layout->addWidget(geodesicGroup);
    layout->addWidget(selectionGroup);
    layout->addWidget(curvatureGroup);
    
    connect(calcSkeletonWeightButton, SIGNAL(clicked()), this, SLOT(calcSkeletonWeightButtonClicked()));
    connect(saveSkeletonWeightButton, SIGNAL(clicked()), this, SLOT(saveSkeletonWeightButtonClicked()));
    connect(geodesicButton, SIGNAL(clicked()), this, SLOT(showGeodesicButtonClicked()));
    connect(geodesicAllButton, SIGNAL(clicked()), this, SLOT(showGeodesicButtonClicked()));
    connect(loadSelectionButton, SIGNAL(clicked()), this, SLOT(loadSelectionButtonClicked()));
    connect(saveSelectionButton, SIGNAL(clicked()), this, SLOT(saveSelectionButtonClicked()));
    connect(clearSelectionButton, SIGNAL(clicked()), this, SLOT(clearSelectionButtonClicked()));
    connect(calcCurvatureButton, SIGNAL(clicked()), this, SLOT(calcCurvatureButtonClicked()));
    
    emit addToolbox(tr("Plush"), toolBox);
}

void PlushPlugin::fileOpened(int _id) {
    // There are some bug loading/saving property in OpenFlipper. we have to handle them by ourselves.
    BaseObjectData *obj;
    PluginFunctions::getObject(_id, obj);
    QString meshName = QFileInfo(obj->name()).baseName();
    
    if (obj->dataType(DATA_TRIANGLE_MESH)) {
        TriMesh *mesh;
        mesh = PluginFunctions::triMesh(obj);
        
        // setup properties
        initProperties(mesh);
        
        // Load curvature
        loadCurvature(mesh, meshName);
        
        // Load selection
        int selectedCount = loadSelection(_id, meshName);
        geodesicEdges->setMaximum(selectedCount*(selectedCount-1)/2);
        
        // Load skeleton
        loadSkeleton(_id);

        // Load bone weight with skeleton
        loadBoneWeight(_id);
    }
}

void PlushPlugin::objectDeleted(int _id) {
    BaseObjectData *obj;
    PluginFunctions::getObject(_id, obj);
    QString meshName = QFileInfo(obj->name()).baseName();
    
    if (obj->dataType(DATA_TRIANGLE_MESH)) {
        TriMesh *mesh;
        mesh = PluginFunctions::triMesh(obj);
        uninitProperties(mesh);
    }
}

void PlushPlugin::initProperties(TriMesh *mesh) {
    mesh->add_property(minCurvatureHandle, "Min Curvature");
    mesh->add_property(maxCurvatureHandle, "Max Curvature");
    mesh->add_property(minCurvatureDirectionHandle, "Min curvature direction");
    mesh->add_property(maxCurvatureDirectionHandle, "Max curvature direction");
    
    mesh->add_property(edgeWeightHandle, "Edge weight");
    mesh->add_property(geodesicDistanceHandle, "Geodesic distance between vertices pair");
    mesh->add_property(geodesicPathHandle, "Geodesic path between vertices pair");
    
    mesh->add_property(skeletonHandle, "Skeleton");
    mesh->add_property(bonesWeightHandle, "Bone weights for each vertex");
    
//    for (VertexIter v_it = mesh->vertices_begin(); v_it != mesh->vertices_end(); v_it++) {
//    }
    for (EdgeIter e_it = mesh->edges_begin(); e_it != mesh->edges_end(); e_it++) {
        mesh->property(edgeWeightHandle, *e_it) = -1;
    }
}

void PlushPlugin::uninitProperties(TriMesh *mesh) {
    if (mesh->property(skeletonHandle) != NULL) {
        delete mesh->property(skeletonHandle);
    }
    for (VertexIter v_it = mesh->vertices_begin(); v_it != mesh->vertices_end(); v_it++) {
        if (mesh->property(bonesWeightHandle, *v_it) != NULL) {
            delete[] mesh->property(bonesWeightHandle, *v_it);
        }
    }
}

void PlushPlugin::calcSkeletonWeightButtonClicked() {
    for (PluginFunctions::ObjectIterator o_it(PluginFunctions::TARGET_OBJECTS); o_it != PluginFunctions::objectsEnd();
         ++o_it) {
        int meshId = o_it->id();
        
        if (o_it->dataType(DATA_TRIANGLE_MESH)) {
            TriMesh *mesh = PluginFunctions::triMesh(*o_it);
            QString meshName = QFileInfo(o_it->name()).baseName();
            
            if (mesh->property(skeletonHandle) == NULL) {
                if (!loadSkeleton(meshId)) {
                    return;
                }
            }
            
            Skeleton *skeleton = mesh->property(skeletonHandle);
            Weight weightGenerator;
            mesh->request_face_normals();
            weightGenerator.computeBoneWeight(mesh, mesh->property(skeletonHandle));
            mesh->release_face_normals();
            
            Eigen::MatrixXd weight = weightGenerator.getWeight();
            for (VertexIter v_it = mesh->vertices_begin(); v_it != mesh->vertices_end(); v_it++) {
                double *w = new double[skeleton->bones.size()];
                for (size_t i = 0; i < skeleton->bones.size(); i++) {
                    w[i] = weight(i, v_it->idx());
                }
                if (mesh->property(bonesWeightHandle, *v_it) != NULL) {
                    delete[] mesh->property(bonesWeightHandle, *v_it);
                }
                mesh->property(bonesWeightHandle, *v_it) = w;
            }
        }
    }
}

void PlushPlugin::saveSkeletonWeightButtonClicked() {
    for (PluginFunctions::ObjectIterator o_it(PluginFunctions::TARGET_OBJECTS); o_it != PluginFunctions::objectsEnd();
         ++o_it) {
        int meshId = o_it->id();
        
        if (o_it->dataType(DATA_TRIANGLE_MESH)) {
            saveBoneWeight(meshId);
        }
    }
}

void PlushPlugin::calcCurvatureButtonClicked() {
    if (thread == NULL) {
        thread = new OpenFlipperThread("calcCurvature");
        connect(thread, SIGNAL(finished(QString)), this, SLOT(finishedJob(QString)));
        connect(thread, SIGNAL(function(QString)), this, SLOT(calcCurvatureThread(QString)), Qt::DirectConnection);
        emit startJob("calcCurvature", "calculate curvature of whole mesh", 0, 100, true);
        thread->start();
        thread->startProcessing();
    }
}

void PlushPlugin::showGeodesicButtonClicked() {
    if (thread == NULL) {
        if (sender() == geodesicButton) {
            showAllPath = false;
        } else {
            showAllPath = true;
        }
        thread = new OpenFlipperThread("calcGeodesic");
        connect(thread, SIGNAL(finished(QString)), this, SLOT(finishedJob(QString)));
        connect(thread, SIGNAL(function(QString)), this, SLOT(showGeodesicThread(QString)), Qt::DirectConnection);
        emit startJob("calcGeodesic", "calculate geodesic of whole mesh", 0, 100, true);
        thread->start();
        thread->startProcessing();
    }
}

void PlushPlugin::saveSelectionButtonClicked() {
    for (PluginFunctions::ObjectIterator o_it(PluginFunctions::TARGET_OBJECTS); o_it != PluginFunctions::objectsEnd();
         ++o_it) {
        int meshId = o_it->id();
        
        if (o_it->dataType(DATA_TRIANGLE_MESH)) {
            TriMesh *mesh = PluginFunctions::triMesh(*o_it);
            
            QString meshName = QFileInfo(o_it->name()).baseName();
            saveSelection(meshId, meshName);
        }
    }
}

void PlushPlugin::loadSelectionButtonClicked() {
    for (PluginFunctions::ObjectIterator o_it(PluginFunctions::TARGET_OBJECTS); o_it != PluginFunctions::objectsEnd();
         ++o_it) {
        int meshId = o_it->id();
        
        if (o_it->dataType(DATA_TRIANGLE_MESH)) {
            TriMesh *mesh = PluginFunctions::triMesh(*o_it);
            
            QString meshName = QFileInfo(o_it->name()).baseName();
            int selectedCount = loadSelection(meshId, meshName);
            
            geodesicEdges->setMaximum(selectedCount*(selectedCount-1)/2);
        }
    }
}

void PlushPlugin::clearSelectionButtonClicked() {
    for (PluginFunctions::ObjectIterator o_it(PluginFunctions::TARGET_OBJECTS); o_it != PluginFunctions::objectsEnd();
         ++o_it) {
        int meshId = o_it->id();
        
        if (o_it->dataType(DATA_TRIANGLE_MESH)) {
            clearSelection(meshId);
        }
    }
}

void PlushPlugin::calcCurvatureThread(QString _jobId) {
    if (PluginFunctions::objectCount() == 0)
    {
        emit log(LOGERR, "Load a model first.");
        return;
    }
    
    for (PluginFunctions::ObjectIterator o_it(PluginFunctions::TARGET_OBJECTS); o_it != PluginFunctions::objectsEnd();
         ++o_it) {
        if (o_it->dataType(DATA_TRIANGLE_MESH)) {
            int meshId = o_it->id();
            calcCurvature(_jobId, meshId);
        }
    }
    emit log(LOGINFO, "Curvature calculation done.");
}

void PlushPlugin::showGeodesicThread(QString _jobId) {
    if (PluginFunctions::objectCount() == 0)
    {
        emit log(LOGERR, "Load a model first.");
        return;
    }

    for (PluginFunctions::ObjectIterator o_it(PluginFunctions::TARGET_OBJECTS); o_it != PluginFunctions::objectsEnd();
         ++o_it) {
        if (o_it->dataType(DATA_TRIANGLE_MESH)) {
            int meshId = o_it->id();
            QString meshName = QFileInfo(o_it->name()).baseName();
            TriMesh *mesh = PluginFunctions::triMesh(*o_it);
            
            IdList selectedVertices;
            selectedVertices = RPC::callFunctionValue<IdList> ("meshobjectselection", "getVertexSelection", meshId);
            geodesicEdges->setMaximum(selectedVertices.size()*(selectedVertices.size()-1)/2);

            std::vector<std::pair<IdList, double> > spanningTree;
            if (calcSpanningTree(_jobId, meshId, spanningTree, selectedVertices)) {
                struct Comparator {
                    bool operator() (std::pair<IdList, double> a,
                                     std::pair<IdList, double> b) {
                        return a.second < b.second;
                    }
                } comparator;
                
                std::sort(spanningTree.begin(), spanningTree.end(), comparator);
                
                IdList edgeList;
                int beginNo = 0;
                int endNo = geodesicEdges->value();
                if (!showAllPath && endNo != 0) {
                    beginNo = endNo - 1;
                }
                for (int i = beginNo; i < endNo; i++) {
                    IdList path = spanningTree[i].first;
                    double weight = spanningTree[i].second;
                    int prevIdx = *path.begin();
                    IdList::iterator vIdx_it = (path.begin()+1);
                    for (; vIdx_it != path.end(); vIdx_it++) {
                        EdgeHandle eh;
                        assert(getEdge(mesh, eh, prevIdx, *vIdx_it));
                        edgeList.push_back(eh.idx());
                        prevIdx = *vIdx_it;
                    }
                    emit log(LOGINFO, QString("Total weight of path #%1: %2").arg(QString::number(i+1), QString::number(weight)));
                }

                RPC::callFunction<int>("meshobjectselection", "clearEdgeSelection", meshId);
                RPC::callFunction<int, IdList>("meshobjectselection", "selectEdges", meshId, edgeList);
            }
        }
    }
    emit log(LOGINFO, "Geodesic calculation done.");
}

bool PlushPlugin::getEdge(TriMesh *mesh, EdgeHandle &_eh, VertexHandle v1, VertexHandle v2) {
    return getEdge(mesh, _eh, v1.idx(), v2.idx());
}

bool PlushPlugin::getEdge(TriMesh *mesh, EdgeHandle &_eh, int v1No, int v2No) {
    VertexHandle v1 = mesh->vertex_handle(v1No);
    VertexHandle v2 = mesh->vertex_handle(v2No);
    for (TriMesh::VertexEdgeIter ve_it = mesh->ve_iter(v1); ve_it; ve_it++)
    {
        HalfedgeHandle heh = mesh->halfedge_handle(*ve_it, 0);
        if (mesh->from_vertex_handle(heh) == v2
        ||  mesh->to_vertex_handle(heh) == v2)
        {
            _eh = *ve_it;
            return true;
        }
    }
    return false;
}

void PlushPlugin::canceledJob(QString _job) {
    // do nothing
    isJobCanceled = true;
    thread = NULL;
}

void PlushPlugin::finishedJob(QString _job) {
    // do nothing
    thread = NULL;
    
    emit finishJob(_job);
}

void PlushPlugin::slotKeyEvent( QKeyEvent* _event ) {
    // Switch pressed keys
//    switch (_event->key())
//    {
//        case Qt::Key_X:
//            if (!pickButton_->isChecked()) {
//                pickButton_->click();
//            }
//            break;
//        case Qt::Key_Z:
//            if (pickButton_->isChecked()) {
//                pickButton_->click();
//            }
//            break;
//        default:
//            break;
//    }
}

Q_EXPORT_PLUGIN2( plushPlugin , PlushPlugin );