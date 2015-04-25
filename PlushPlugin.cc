#include "PlushPlugin.hh"

#include <MeshTools/MeshSelectionT.hh>
#include <ObjectTypes/PolyLine/PolyLine.hh>

PlushPlugin::PlushPlugin()
{
    requiredPlugins = new std::vector<char*>();
    
    m_currentJobId = "";
    
    m_patternGenerator = NULL;
    m_triMeshObj = NULL;
}

void PlushPlugin::initializePlugin()
{
    // Create the Toolbox Widget
    QWidget* toolBox = new QWidget();
    QGridLayout* layout = new QGridLayout(toolBox);
    
    // Create button that can be toggled
    // to (de)activate plugin's picking mode
    QGroupBox *geodesicGroup = new QGroupBox(tr("Geodesic paths"));
    QLabel *geodesicNumberLabel = new QLabel(tr("#"));
    geodesicNumPaths = new QSpinBox();
    geodesicNumPaths->setMinimum(0);
    geodesicNumPaths->setValue(999);
    geodesicElimination = new QCheckBox(tr("Eliminate crossover paths"));
    geodesicElimination->setChecked(true);
    geodesicShowSingleButton = new QPushButton(tr("Show single path"));
    geodesicShowAllButton = new QPushButton(tr("Show all path"));
    geodesicCalcButton = new QPushButton(tr("Calculate"));
    QVBoxLayout *geodesicLayout = new QVBoxLayout;
    QHBoxLayout *geodesicRow1Layout = new QHBoxLayout;
    QHBoxLayout *geodesicRow2Layout = new QHBoxLayout;
    QHBoxLayout *geodesicRow3Layout = new QHBoxLayout;
    geodesicRow1Layout->addWidget(geodesicNumberLabel);
    geodesicRow1Layout->addWidget(geodesicNumPaths);
    geodesicRow1Layout->addWidget(geodesicElimination);
    geodesicRow2Layout->addWidget(geodesicShowSingleButton);
    geodesicRow2Layout->addWidget(geodesicShowAllButton);
    geodesicRow3Layout->addWidget(geodesicCalcButton);
    geodesicLayout->addLayout(geodesicRow1Layout);
    geodesicLayout->addLayout(geodesicRow2Layout);
    geodesicLayout->addLayout(geodesicRow3Layout);
    geodesicGroup->setLayout(geodesicLayout);
    
    QGroupBox *flatteningGroup = new QGroupBox(tr("Flattening"));
    QPushButton *calcFlattenButton = new QPushButton(tr("Calculate flattened graph"));
    QPushButton *showFlattenButton = new QPushButton(tr("Show graph"));
    QHBoxLayout *flatteningLayout = new QHBoxLayout;
    flatteningLayout->addWidget(calcFlattenButton);
    flatteningLayout->addWidget(showFlattenButton);
    flatteningGroup->setLayout(flatteningLayout);
    
    QGroupBox *selectionGroup = new QGroupBox(tr("Selection"));
    QPushButton *loadSelectionButton = new QPushButton(tr("Load"));
    QPushButton *saveSelectionButton = new QPushButton(tr("Save"));
    QPushButton *clearSelectionButton = new QPushButton(tr("Clear"));
    QPushButton *calcSelectionButton = new QPushButton(tr("Calculate selction"));
    QVBoxLayout *selectionLayout = new QVBoxLayout;
    QHBoxLayout *selectionRow1Layout = new QHBoxLayout;
    QHBoxLayout *selectionRow2Layout = new QHBoxLayout;
    selectionRow1Layout->addWidget(loadSelectionButton);
    selectionRow1Layout->addWidget(saveSelectionButton);
    selectionRow1Layout->addWidget(clearSelectionButton);
    selectionRow2Layout->addWidget(calcSelectionButton);
    selectionLayout->addLayout(selectionRow1Layout);
    selectionLayout->addLayout(selectionRow2Layout);
    selectionGroup->setLayout(selectionLayout);
    
    QGroupBox *curvatureGroup = new QGroupBox(tr("Curvature"));
    QPushButton *calcCurvatureButton = new QPushButton(tr("Calculate"));
    QHBoxLayout *curvatureLayout = new QHBoxLayout;
    curvatureLayout->addWidget(calcCurvatureButton);
    curvatureGroup->setLayout(curvatureLayout);
    
    QGroupBox *skeletonWeightGroup = new QGroupBox(tr("Skeleton weight"));
    QPushButton *calcSkeletonWeightButton = new QPushButton(tr("Calculate"));
    QHBoxLayout *skeletonWeightLayout = new QHBoxLayout;
    skeletonWeightLayout->addWidget(calcSkeletonWeightButton);
    skeletonWeightGroup->setLayout(skeletonWeightLayout);
    
    layout->addWidget(geodesicGroup);
    layout->addWidget(selectionGroup);
    layout->addWidget(flatteningGroup);
    layout->addWidget(curvatureGroup);
    layout->addWidget(skeletonWeightGroup);
    
    connect(calcSkeletonWeightButton, SIGNAL(clicked()), this, SLOT(calcSkeletonWeightButtonClicked()));
    connect(geodesicShowSingleButton, SIGNAL(clicked()), this, SLOT(showGeodesicButtonClicked()));
    connect(geodesicShowAllButton, SIGNAL(clicked()), this, SLOT(showGeodesicButtonClicked()));
    connect(geodesicCalcButton, SIGNAL(clicked()), this, SLOT(calcGeodesicButtonClicked()));
    connect(loadSelectionButton, SIGNAL(clicked()), this, SLOT(loadSelectionButtonClicked()));
    connect(saveSelectionButton, SIGNAL(clicked()), this, SLOT(saveSelectionButtonClicked()));
    connect(clearSelectionButton, SIGNAL(clicked()), this, SLOT(clearSelectionButtonClicked()));
    connect(calcSelectionButton, SIGNAL(clicked()), this, SLOT(calcSelectionButtonClicked()));
    connect(calcFlattenButton, SIGNAL(clicked()), this, SLOT(calcFlattenedGraphButtonClicked()));
    connect(showFlattenButton, SIGNAL(clicked()), this, SLOT(showFlattenedGrpahButtonClicked()));
    connect(calcCurvatureButton, SIGNAL(clicked()), this, SLOT(calcCurvatureButtonClicked()));
    
    emit addToolbox(tr("Plush"), toolBox);

    // Register keys
    emit registerKey(Qt::Key_X,     Qt::NoModifier, "Next path");
    emit registerKey(Qt::Key_Z,     Qt::NoModifier, "Previous path");
}

PlushPlugin::~PlushPlugin() {
    delete requiredPlugins;
}
void PlushPlugin::pluginsInitialized() {
}

void PlushPlugin::fileOpened(int _id) {
    // There are some bug loading/saving property in OpenFlipper. we have to handle them by ourselves.
    BaseObjectData *obj;
    PluginFunctions::getObject(_id, obj);
    
    if (obj->dataType(DATA_TRIANGLE_MESH)) {
        if (!m_patternGenerator) {
            m_triMeshObj = PluginFunctions::triMeshObject(obj);
            TriMesh *mesh = PluginFunctions::triMesh(obj);
            QString meshName = QFileInfo(obj->name()).baseName();
            m_patternGenerator = new PlushPatternGenerator(mesh, meshName);
            
            // Load selection
            int selectedCount = loadSelection(mesh, meshName);
            geodesicNumPaths->setMaximum(selectedCount*(selectedCount-1)/2);
            
            // When user click on "cancel" button, cancel this job
            connect(this, SIGNAL(cancelingJob()), m_patternGenerator, SLOT(cancelJob()));
            // Get progress info from m_patternGenerator, then re-emit it to OpenFlipper system inside receiveJobState
            connect(m_patternGenerator, SIGNAL(setJobState(int)), this, SLOT(receiveJobState(int)));
            // Propagate messages from m_patternGenerator to OpenFlipper
            connect(m_patternGenerator, SIGNAL(log(int, QString)), this, SLOT(receiveLog(int, QString)));
        }
    }
}

void PlushPlugin::objectDeleted(int _id) {
    if (m_triMeshObj && _id == m_triMeshObj->id()) {
        disconnect(this, SIGNAL(cancelingJob()), m_patternGenerator, SLOT(cancelJob()));
        disconnect(m_patternGenerator, SIGNAL(setJobState(int)), this, SLOT(receiveJobState(int)));
        disconnect(m_patternGenerator, SIGNAL(log(int, QString)), this, SLOT(receiveLog(int, QString)));

        delete m_patternGenerator;
        
        m_patternGenerator = NULL;
        m_triMeshObj = NULL;
    }
}

int PlushPlugin::loadSelection(TriMesh *mesh, QString meshName) {
    QFile file(meshName + "_selection.txt");
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
        emit log(LOGERR, QString("Unable to read file %1").arg(meshName + "_selection.txt"));
        return 0;
    }
    
    QTextStream in(&file);
    
    std::vector<int> selectedVertices;
    QString v;
    while(!in.atEnd()) {
        v = in.readLine();
        bool ok;
        int vId = v.toInt(&ok);
        if (ok) {
            selectedVertices.push_back(vId);
        }
    }
    
    MeshSelection::selectVertices(mesh, selectedVertices);
    emit updatedObject(m_triMeshObj->id(), UPDATE_SELECTION);
    return selectedVertices.size();
}

void PlushPlugin::saveSelection(TriMesh *mesh, QString meshName) {
    std::vector<int> selectedVertices = MeshSelection::getVertexSelection(mesh);
    
    // Prepare file for saving data
    QFile file(meshName + "_selection.txt");
    file.open(QIODevice::WriteOnly | QIODevice::Text);
    QTextStream out(&file);
    
    for (std::vector<int>::iterator v_it = selectedVertices.begin(); v_it != selectedVertices.end(); v_it++) {
        out << *v_it << "\n";
    }
    file.close();
}

void PlushPlugin::clearSelection(TriMesh *mesh) {
    MeshSelection::clearVertexSelection(mesh);
    emit updatedObject(m_triMeshObj->id(), UPDATE_SELECTION);
}

void PlushPlugin::calcSelection(TriMesh *mesh) {
    std::vector<VertexHandle> targetVertices;
    m_patternGenerator->calcSelection(targetVertices);
    
    std::vector<int> targetVerticesId;
    for (auto v : targetVertices) {
        targetVerticesId.push_back(v.idx());
    }
    MeshSelection::selectVertices(mesh, targetVerticesId);
    emit updatedObject(m_triMeshObj->id(), UPDATE_SELECTION);
}

void PlushPlugin::calcSkeletonWeightButtonClicked() {
    if (!checkIfGeneratorExist()) {
        return;
    }

    m_patternGenerator->calcSkeletonWeight();
    m_patternGenerator->saveBoneWeight();
}

void PlushPlugin::calcCurvatureButtonClicked() {
    if (!checkIfGeneratorExist()) {
        return;
    }
    
    m_currentJobId = "calcCurvature";
    OpenFlipperThread *thread = new OpenFlipperThread(m_currentJobId);
    connect(thread, SIGNAL(finished(QString)), this, SIGNAL(finishJob(QString)));
    connect(thread, SIGNAL(function(QString)), this, SLOT(calcCurvatureThread()), Qt::DirectConnection);
    
    // Custom handler to set m_currentJobId to "" after finishing job. Also updating screen.
    connect(thread, SIGNAL(finished(QString)), this, SLOT(finishedJobHandler()));
    
    emit startJob(m_currentJobId, "calculate curvature of whole mesh", 0, 100, true);
    thread->start();
    thread->startProcessing();
}

void PlushPlugin::calcGeodesicButtonClicked() {
    if (!checkIfGeneratorExist()) {
        return;
    }

    m_currentJobId = "calcGeodesic";
    OpenFlipperThread *thread = new OpenFlipperThread(m_currentJobId);
    connect(thread, SIGNAL(finished(QString)), this, SIGNAL(finishJob(QString)));
    connect(thread, SIGNAL(function(QString)), this, SLOT(calcGeodesicThread()), Qt::DirectConnection);
    
    // Custom handler to set m_currentJobId to "" after finishing job.
    connect(thread, SIGNAL(finished(QString)), this, SLOT(finishedJobHandler()));
    
    emit startJob(m_currentJobId, "calculate geodesic of whole mesh", 0, 100, true);
    thread->start();
    thread->startProcessing();
}

void PlushPlugin::calcFlattenedGraphButtonClicked() {
    std::set<EdgeHandle> *seams = m_patternGenerator->getSeams();
    if (!seams) {
        emit log(LOGERR, "No seams. Use show geodesic path.");
        return;
    }
    
    m_currentJobId = "calcFlattenedGraph";
    OpenFlipperThread *thread = new OpenFlipperThread(m_currentJobId);
    connect(thread, SIGNAL(finished(QString)), this, SIGNAL(finishJob(QString)));
    connect(thread, SIGNAL(function(QString)), this, SLOT(calcFlattenedGraphThread()), Qt::DirectConnection);
    
    // Custom handler to set m_currentJobId to "" after finishing job.
    connect(thread, SIGNAL(finished(QString)), this, SLOT(finishedJobHandler()));
    
    emit startJob(m_currentJobId, "Calculate flattened graph using spanning tree", 0, 100, true);
    thread->start();
    thread->startProcessing();
}

void PlushPlugin::showGeodesicButtonClicked() {
    if (!checkIfGeneratorExist()) {
        return;
    }
    
    if (sender() == geodesicShowSingleButton) {
        showAllPath = false;
    } else {
        showAllPath = true;
    }
    
    TriMesh *mesh = m_triMeshObj->mesh();
    
    std::vector<int> selectedVerticesId;
    selectedVerticesId = MeshSelection::getVertexSelection(mesh);
    geodesicNumPaths->setMaximum(selectedVerticesId.size()*(selectedVerticesId.size()-1)/2);
    
    std::vector<VertexHandle> selectedVertices;
    for (size_t i = 0; i < selectedVerticesId.size(); i++) {
        selectedVertices.push_back(mesh->vertex_handle(selectedVerticesId[i]));
    }
    
    MeshSelection::clearEdgeSelection(mesh);
    m_patternGenerator->calcSeams(selectedVertices, geodesicNumPaths->value(), geodesicElimination->isChecked(), showAllPath);
    m_patternGenerator->calcCircularSeams(mesh);
    std::vector<TriMesh> *subMeshes = m_patternGenerator->getFlattenedMeshes();
    if (subMeshes) {
        for (size_t i = 0; i < subMeshes->size(); i++) {
            m_patternGenerator->calcCircularSeams(&subMeshes->at(i));
        }
    }
    std::set<EdgeHandle> *seams = m_patternGenerator->getSeams();
    if (seams) {
        std::vector<int> edgeList;
        for (auto e_it = seams->begin(); e_it != seams->end(); e_it++) {
            edgeList.push_back(e_it->idx());
        }
        MeshSelection::selectEdges(mesh, edgeList);
        emit updatedObject(m_triMeshObj->id(), UPDATE_SELECTION);
    }
}

void PlushPlugin::showFlattenedGrpahButtonClicked() {
    std::vector<TriMesh> *flattenedMeshes = m_patternGenerator->getFlattenedMeshes();
    if (!flattenedMeshes) {
        emit log(LOGERR, "Calculate flattened meshes first.");
        return;
    }
    
    // 1 means right viewer
    int viewerId = 1;
    // We don't want to show both models and polylines at same the time.
    // So we place polylines to a very far place so that it will not be shown on original viewer
    OpenMesh::Vec3d eye(100, 0, 1);
    OpenMesh::Vec3d far(100, 0, 0);
    OpenMesh::Vec3d right(0, 1, 0);
    
    RPC::callFunction ("core", "multiViewMode", 1);
    PluginFunctions::setActiveExaminer(viewerId);
    PluginFunctions::setFixedView(PluginFunctions::VIEW_FRONT, viewerId);
    PluginFunctions::allowRotation(false, viewerId);
    PluginFunctions::orthographicProjection(viewerId);
    // Set this viewer to a far place
    PluginFunctions::lookAt(eye, far, right);


    // Add new object and copy m_flattenedGraph into it
    int n = flattenedMeshes->size();
    int *newTriMeshIds = new int[n];
    
    std::vector<int> seamsId;
    for (int i = 0; i < n; i++) {
        emit addEmptyObject(DATA_TRIANGLE_MESH, newTriMeshIds[i]);

        // Get the newly created object
        TriMeshObject *object = 0;
        PluginFunctions::getObject(newTriMeshIds[i], object);
        TriMesh *subMesh = object->mesh();
        OpenMesh::VPropHandleT<VertexHandle> inverseMappingHandle = PlushPatternGenerator::getInverseMappingHandle(subMesh);
        OpenMesh::MPropHandleT< std::set<EdgeHandle> > seamsHandle = PlushPatternGenerator::getSeamsHandle(subMesh);
        
        TriMesh &oldSubMesh = flattenedMeshes->at(i);
        OpenMesh::VPropHandleT<VertexHandle> oldInverseMappingHandle = PlushPatternGenerator::getInverseMappingHandle(&oldSubMesh);
        OpenMesh::MPropHandleT< std::set<EdgeHandle> > oldSeamsHandle = PlushPatternGenerator::getSeamsHandle(&oldSubMesh);
        std::map<VertexHandle, VertexHandle> oldToNewMapping;

        for (VertexIter v_it = oldSubMesh.vertices_begin(); v_it != oldSubMesh.vertices_end(); v_it++) {
            TriMesh::Point p = oldSubMesh.point(*v_it);
            VertexHandle newV = subMesh->add_vertex(p + far);
            oldToNewMapping.emplace(*v_it, newV);
            
            // Update inverse map
            subMesh->property(inverseMappingHandle, newV) = oldSubMesh.property(oldInverseMappingHandle, *v_it);
        }

        for (FaceIter f_it = oldSubMesh.faces_begin(); f_it != oldSubMesh.faces_end(); f_it++) {
            HalfedgeHandle heh = oldSubMesh.halfedge_handle(*f_it);
            VertexHandle v1 = oldSubMesh.from_vertex_handle(heh);
            VertexHandle v2 = oldSubMesh.to_vertex_handle(heh);
            VertexHandle v3 = oldSubMesh.to_vertex_handle(oldSubMesh.next_halfedge_handle(heh));
            
            std::vector<VertexHandle> face_vhs;
            face_vhs.push_back(oldToNewMapping[v1]);
            face_vhs.push_back(oldToNewMapping[v2]);
            face_vhs.push_back(oldToNewMapping[v3]);
            FaceHandle newF = subMesh->add_face(face_vhs);
            
            // Set properties
            subMesh->set_color(newF, oldSubMesh.color(*f_it));
        }
        
        // Update seams
        std::vector<int> subSeamsId;

        std::set<EdgeHandle> &seams = subMesh->property(seamsHandle);
        std::set<EdgeHandle> &oldSeams = oldSubMesh.property(oldSeamsHandle);
        for (auto e_it = oldSeams.begin(); e_it != oldSeams.end(); e_it++) {
            HalfedgeHandle old_heh = oldSubMesh.halfedge_handle(*e_it, 0);
            VertexHandle oldV1 = oldSubMesh.from_vertex_handle(old_heh);
            VertexHandle oldV2 = oldSubMesh.to_vertex_handle(old_heh);
            EdgeHandle he;
            bool edgeExist = PlushPatternGenerator::getEdge(subMesh, he, oldToNewMapping[oldV1], oldToNewMapping[oldV2]);
            assert(edgeExist);
            
            // Update inverse map
            seams.insert(he);
            subSeamsId.push_back(he.idx());
            
            VertexHandle oldOriginalV1 = oldSubMesh.property(oldInverseMappingHandle, oldV1);
            VertexHandle oldOriginalV2 = oldSubMesh.property(oldInverseMappingHandle, oldV2);
            EdgeHandle original_he;
            edgeExist = PlushPatternGenerator::getEdge(m_triMeshObj->mesh(), original_he, oldOriginalV1, oldOriginalV2);
            assert(edgeExist);
            seamsId.push_back(original_he.idx());
        }
        MeshSelection::selectEdges(subMesh, subSeamsId);

        MeshSelection::selectBoundaryEdges(subMesh);

        subMesh->request_face_normals();
        subMesh->update_face_normals();
    
        // Replace old sub mesh with new sub mesh to m_flattenedGraph
        flattenedMeshes->at(i) = *subMesh;
    }
    MeshSelection::selectEdges(m_triMeshObj->mesh(), seamsId);
    emit updatedObject(m_triMeshObj->id(), UPDATE_SELECTION_VERTICES);
}

void PlushPlugin::saveSelectionButtonClicked() {
    if (!checkIfGeneratorExist()) {
        return;
    }

    QString meshName = QFileInfo(m_triMeshObj->name()).baseName();
    saveSelection(m_triMeshObj->mesh(), meshName);
}

void PlushPlugin::loadSelectionButtonClicked() {
    if (!checkIfGeneratorExist()) {
        return;
    }
    
    QString meshName = QFileInfo(m_triMeshObj->name()).baseName();
    int selectedCount = loadSelection(m_triMeshObj->mesh(), meshName);

    geodesicNumPaths->setMaximum(selectedCount*(selectedCount-1)/2);
}

void PlushPlugin::clearSelectionButtonClicked() {
    if (!checkIfGeneratorExist()) {
        return;
    }
    
    clearSelection(m_triMeshObj->mesh());
}

void PlushPlugin::calcSelectionButtonClicked() {
    if (!checkIfGeneratorExist()) {
        return;
    }
    
    calcSelection(m_triMeshObj->mesh());
}

void PlushPlugin::calcCurvatureThread() {
    if (!checkIfGeneratorExist()) {
        return;
    }
    
    m_patternGenerator->calcCurvature();
    emit log(LOGINFO, "Curvature calculation done.");
}

void PlushPlugin::calcGeodesicThread() {
    if (!checkIfGeneratorExist()) {
        return;
    }
    
    TriMesh *mesh = m_triMeshObj->mesh();
    
    std::vector<int> selectedVerticesId;
    selectedVerticesId = MeshSelection::getVertexSelection(mesh);
    geodesicNumPaths->setMaximum(selectedVerticesId.size()*(selectedVerticesId.size()-1)/2);
    
    std::vector<VertexHandle> selectedVertices;
    for (size_t i = 0; i < selectedVerticesId.size(); i++) {
        selectedVertices.push_back(mesh->vertex_handle(selectedVerticesId[i]));
    }
    m_patternGenerator->calcGeodesic(selectedVertices);
    m_patternGenerator->saveGeodesic(selectedVertices);

    emit log(LOGINFO, "Geodesic calculation done.");
}

void PlushPlugin::calcFlattenedGraphThread() {
    m_patternGenerator->calcFlattenedGraph();
}

void PlushPlugin::canceledJob(QString _job) {
    emit cancelingJob();
}

void PlushPlugin::finishedJobHandler() {
    m_currentJobId = "";
    emit updatedObject(m_triMeshObj->id(), UPDATE_SELECTION);
}
void PlushPlugin::receiveJobState(int state) {
    if (m_currentJobId != "") {
        emit setJobState(m_currentJobId, state);
    }
}

void PlushPlugin::receiveLog(int type, QString msg) {
    switch (type) {
        case PlushPatternGenerator::LOGERR:
            emit log(LOGERR, msg);
            break;
        case PlushPatternGenerator::LOGINFO:
            emit log(LOGINFO, msg);
            break;
    }
}

void PlushPlugin::slotKeyEvent( QKeyEvent* _event ) {
    // Switch pressed keys
    int currentEdgeNo = geodesicNumPaths->value();
    switch (_event->key())
    {
        case Qt::Key_X:
            if (currentEdgeNo < geodesicNumPaths->maximum()) {
                geodesicNumPaths->setValue(currentEdgeNo+1);
            }
            geodesicShowSingleButton->click();
            break;
        case Qt::Key_Z:
            if (currentEdgeNo > 0) {
                geodesicNumPaths->setValue(currentEdgeNo-1);
            }
            geodesicShowSingleButton->click();
            break;
        default:
            break;
    }
}

bool PlushPlugin::checkIfGeneratorExist() {
    if (!m_patternGenerator || !m_triMeshObj)
    {
        emit log(LOGERR, "Load a model first.");
        return false;
    }
    return true;
}

Q_EXPORT_PLUGIN2( plushPlugin , PlushPlugin );
