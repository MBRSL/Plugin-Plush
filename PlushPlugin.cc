#include "PlushPlugin.hh"

#include <MeshTools/MeshSelectionT.hh>

PlushPlugin::PlushPlugin()
{
    requiredPlugins = new std::vector<char*>();
    
    m_currentJobId = "";
    
    m_patternGenerator = NULL;
    m_triMeshObj = NULL;
}

void PlushPlugin::initializePlugin()
{
    // Register keys
    emit registerKey(Qt::Key_X,     Qt::NoModifier, "Normal mode");
    emit registerKey(Qt::Key_Z,     Qt::NoModifier, "Plush picking mode");
}

PlushPlugin::~PlushPlugin() {
    delete requiredPlugins;
}
void PlushPlugin::pluginsInitialized() {
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
    
    if (obj->dataType(DATA_TRIANGLE_MESH)) {
        if (!m_patternGenerator) {
            m_triMeshObj = PluginFunctions::triMeshObject(obj);
            TriMesh *mesh = PluginFunctions::triMesh(obj);
            QString meshName = QFileInfo(obj->name()).baseName();
            m_patternGenerator = new PlushPatternGenerator(mesh, meshName);
            
            // Load selection
            int selectedCount = loadSelection(mesh, meshName);
            geodesicEdges->setMaximum(selectedCount*(selectedCount-1)/2);
            
            // When user click on "cancel" button, cancel this job
            connect(this, SIGNAL(cancelingJob()), m_patternGenerator, SLOT(cancelJob()));
            // Get progress info from m_patternGenerator, then re-emit it to OpenFlipper system inside receiveJobState
            connect(m_patternGenerator, SIGNAL(setJobState(int)), this, SLOT(receiveJobState(int)));
        }
    }
}

void PlushPlugin::objectDeleted(int _id) {
    if (!m_triMeshObj && _id == m_triMeshObj->id()) {
        disconnect(this, SIGNAL(cancelingJob()), m_patternGenerator, SLOT(cancelJob()));
        disconnect(m_patternGenerator, SIGNAL(setJobState(int)), this, SLOT(receiveJobState(int)));

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
}

void PlushPlugin::saveSkeletonWeightButtonClicked() {
    if (!checkIfGeneratorExist()) {
        return;
    }

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

void PlushPlugin::showGeodesicButtonClicked() {
    if (!checkIfGeneratorExist()) {
        return;
    }

    if (sender() == geodesicButton) {
        showAllPath = false;
    } else {
        showAllPath = true;
    }

    m_currentJobId = "calcGeodesic";
    OpenFlipperThread *thread = new OpenFlipperThread(m_currentJobId);
    connect(thread, SIGNAL(finished(QString)), this, SIGNAL(finishJob(QString)));
    connect(thread, SIGNAL(function(QString)), this, SLOT(showGeodesicThread()), Qt::DirectConnection);
    
    // Custom handler to set m_currentJobId to "" after finishing job. Also updating screen.
    connect(thread, SIGNAL(finished(QString)), this, SLOT(finishedJobHandler()));
    
    emit startJob(m_currentJobId, "calculate geodesic of whole mesh", 0, 100, true);
    thread->start();
    thread->startProcessing();
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

    geodesicEdges->setMaximum(selectedCount*(selectedCount-1)/2);
}

void PlushPlugin::clearSelectionButtonClicked() {
    if (!checkIfGeneratorExist()) {
        return;
    }
    
    clearSelection(m_triMeshObj->mesh());
}

void PlushPlugin::calcCurvatureThread() {
    if (!checkIfGeneratorExist()) {
        return;
    }
    
    m_patternGenerator->calcCurvature();
    emit log(LOGINFO, "Curvature calculation done.");
}

void PlushPlugin::showGeodesicThread() {
    if (!checkIfGeneratorExist()) {
        return;
    }

    TriMesh *mesh = m_triMeshObj->mesh();

    std::vector<int> selectedVertices;
    selectedVertices = MeshSelection::getVertexSelection(mesh);
    geodesicEdges->setMaximum(selectedVertices.size()*(selectedVertices.size()-1)/2);

    std::vector<EdgeHandle> spanningTree;
    if (m_patternGenerator->calcSpanningTree(spanningTree, selectedVertices, geodesicEdges->value(), showAllPath)) {
        std::vector<int> edgeList;
        for (size_t i = 0; i < spanningTree.size(); i++) {
            edgeList.push_back(spanningTree[i].idx());
        }

        MeshSelection::clearEdgeSelection(mesh);
        MeshSelection::selectEdges(mesh, edgeList);
    }
    emit log(LOGINFO, "Geodesic calculation done.");
}

void PlushPlugin::canceledJob(QString _job) {
    emit cancelingJob();
}

void PlushPlugin::finishedJobHandler() {
    m_currentJobId = "";
    
    // Update visual changes outside thread, otherwise it will crash
    emit updatedObject(m_triMeshObj->id(), UPDATE_SELECTION_EDGES);
}
void PlushPlugin::receiveJobState(int state) {
    if (m_currentJobId != "") {
        emit setJobState(m_currentJobId, state);
    }
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

bool PlushPlugin::checkIfGeneratorExist() {
    if (!m_patternGenerator || !m_triMeshObj)
    {
        emit log(LOGERR, "Load a model first.");
        return false;
    }
    return true;
}

Q_EXPORT_PLUGIN2( plushPlugin , PlushPlugin );
