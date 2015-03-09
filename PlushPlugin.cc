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
    geodesicElimination = new QCheckBox(tr("Eliminate crossover paths"));
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
    QPushButton *flattenButton = new QPushButton(tr("Add plane"));
    QHBoxLayout *flatteningLayout = new QHBoxLayout;
    flatteningLayout->addWidget(flattenButton);
    flatteningGroup->setLayout(flatteningLayout);
    
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
    connect(flattenButton, SIGNAL(clicked()), this, SLOT(flattenButtonClicked()));
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
    emit updatedObject(m_triMeshObj->id(), UPDATE_SELECTION_EDGES);
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
    emit updatedObject(m_triMeshObj->id(), UPDATE_SELECTION_EDGES);
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
    
    // Custom handler to set m_currentJobId to "" after finishing job. Also updating screen.
    connect(thread, SIGNAL(finished(QString)), this, SLOT(finishedJobHandler()));
    
    emit startJob(m_currentJobId, "calculate geodesic of whole mesh", 0, 100, true);
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
    
    if (m_patternGenerator->calcSpanningTree(selectedVertices, geodesicNumPaths->value(), geodesicElimination->isChecked(), showAllPath)) {
        std::vector<int> edgeList;
        for (size_t i = 0; i < m_patternGenerator->m_spanningTree.size(); i++) {
            edgeList.push_back(m_patternGenerator->m_spanningTree[i].idx());
        }
        
        MeshSelection::clearEdgeSelection(mesh);
        MeshSelection::selectEdges(mesh, edgeList);
        emit updatedObject(m_triMeshObj->id(), UPDATE_SELECTION_EDGES);
    }
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

void PlushPlugin::flattenButtonClicked() {
    if (m_patternGenerator->m_spanningTree.size() == 0) {
        emit log(LOGERR, "No spanning tree calculated. Use show geodesic path.");
        return;
    }
    
    // Remove previous polyline
    for (PluginFunctions::ObjectIterator o_it(PluginFunctions::ALL_OBJECTS); o_it != PluginFunctions::objectsEnd();
         ++o_it) {
        int objId = o_it->id();
        
        if (o_it->dataType(DATA_POLY_LINE)) {
            emit deleteObject(objId);
        }
    }
    
    TriMesh *mesh = m_patternGenerator->m_mesh;

    std::vector<int> selectedVerticesId;
    selectedVerticesId = MeshSelection::getVertexSelection(mesh);
    
    std::vector<VertexHandle> selectedVertices;
    for (size_t i = 0; i < selectedVerticesId.size(); i++) {
        selectedVertices.push_back(mesh->vertex_handle(selectedVerticesId[i]));
    }

    std::vector< std::vector<HalfedgeHandle> > loops;
    m_patternGenerator->getLoops(loops, selectedVertices);
    
    for (size_t i = 0; i < loops.size(); i++) {
        int objId;
        emit addEmptyObject(DATA_POLY_LINE, objId);
        
        // Get the newly created object
        PolyLineObject *object = 0;
        PluginFunctions::getObject(objId, object);
        PolyLine *polyLine = object->line();

        polyLine->add_point(mesh->point(mesh->from_vertex_handle(loops[i][0])));
        for (size_t j = 0; j < loops[i].size(); j++) {
            HalfedgeHandle heh = loops[i][j];
            polyLine->add_point(mesh->point(mesh->to_vertex_handle(heh)));
        }
    }
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

void PlushPlugin::canceledJob(QString _job) {
    emit cancelingJob();
}

void PlushPlugin::finishedJobHandler() {
    m_currentJobId = "";
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
