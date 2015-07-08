#include "PlushPlugin.hh"

#include <ACG/Scenegraph/TextNode.hh>
#include <ACG/Utils/ColorCoder.hh>

#include <MeshTools/MeshSelectionT.hh>
#include <MeshTools/MeshInfoT.hh>

#include <ObjectTypes/PolyLine/PolyLine.hh>

#include <sys/time.h>
double get_wall_time(){
    struct timeval time;
    if (gettimeofday(&time,NULL)){
        //  Handle error
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}
double get_cpu_time(){
    return (double)clock() / CLOCKS_PER_SEC;
}

double cpu0, cpu1;

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
    QGroupBox *seamGroup = new QGroupBox(tr("Seams"));
    seamNumPaths = new QSpinBox();
    seamNumPaths->setMinimum(0);
    seamElimination = new QCheckBox(tr("Eliminate crossover paths"));
//    geodesicElimination->setChecked(true);
    seamShowSingleButton = new QPushButton(tr("Single seam"));
    seamShowAllButton = new QPushButton(tr("All seams"));
    seamLocalButton = new QPushButton(tr("Local seams"));
    QPushButton *seam_from_selection_button = new QPushButton(tr("Seam from selection"));
    QPushButton *seamLoadButton = new QPushButton(tr("Load"));
    QPushButton *seamSaveButton = new QPushButton(tr("Save"));
    seamMergeIsStep = new QCheckBox();
    seamMergeThreshold = new QDoubleSpinBox();
    seamMergeThreshold->setDecimals(5);
    seamMergeThreshold->setMinimum(0.00001);
    seamMergeThreshold->setMaximum(1);
    seamMergeThreshold->setValue(0.05);
    seamMergeThreshold->setSingleStep(0.005);
    seamMergeButton = new QPushButton(tr("Merge"));
    QVBoxLayout *seamLayout = new QVBoxLayout;
    QHBoxLayout *seamRow1Layout = new QHBoxLayout;
    QHBoxLayout *seamRow2Layout = new QHBoxLayout;
    QHBoxLayout *seamRow3Layout = new QHBoxLayout;
    QHBoxLayout *seamRow4Layout = new QHBoxLayout;
    seamRow1Layout->addWidget(new QLabel(tr("#")));
    seamRow1Layout->addWidget(seamNumPaths);
    seamRow1Layout->addWidget(seamElimination);
    seamRow2Layout->addWidget(seamShowSingleButton);
    seamRow2Layout->addWidget(seamShowAllButton);
    seamRow2Layout->addWidget(seamLocalButton);
    seamRow3Layout->addWidget(seam_from_selection_button);
    seamRow3Layout->addWidget(seamLoadButton);
    seamRow3Layout->addWidget(seamSaveButton);
    seamRow4Layout->addWidget(new QLabel(tr("Is step")));
    seamRow4Layout->addWidget(seamMergeIsStep);
    seamRow4Layout->addWidget(seamMergeThreshold);
    seamRow4Layout->addWidget(seamMergeButton);
    seamLayout->addLayout(seamRow1Layout);
    seamLayout->addLayout(seamRow2Layout);
    seamLayout->addLayout(seamRow3Layout);
    seamLayout->addLayout(seamRow4Layout);
    seamGroup->setLayout(seamLayout);
    
    QGroupBox *subsetGroup = new QGroupBox(tr("Subset"));
    QHBoxLayout *subsetLayout = new QHBoxLayout;
    subset_merge_threshold = new QDoubleSpinBox();
    subset_merge_threshold->setDecimals(5);
    subset_merge_threshold->setMinimum(0.00001);
    subset_merge_threshold->setMaximum(1);
    subset_merge_threshold->setValue(0.05);
    subset_merge_threshold->setSingleStep(0.005);
    QPushButton *subset_calc_button = new QPushButton(tr("Calc subset"));
    QPushButton *subset_show_button = new QPushButton(tr("Show subset"));
    subsetLayout->addWidget(subset_merge_threshold);
    subsetLayout->addWidget(subset_calc_button);
    subsetLayout->addWidget(subset_show_button);
    subsetGroup->setLayout(subsetLayout);
    
    QGroupBox *geodesicGroup = new QGroupBox(tr("Geodesic"));
    geodesic_distance_coef = new QDoubleSpinBox();
    geodesic_curvature_coef = new QDoubleSpinBox();
    geodesic_skeleton_coef = new QDoubleSpinBox();
    geodesic_smoothness_coef = new QDoubleSpinBox();
    geodesic_distance_coef->setValue(1.0);
    geodesic_curvature_coef->setValue(1.2);
    geodesic_skeleton_coef->setValue(0.2);
    geodesic_smoothness_coef->setValue(0.0);
    geodesic_distance_coef->setSingleStep(0.1);
    geodesic_curvature_coef->setSingleStep(0.1);
    geodesic_skeleton_coef->setSingleStep(0.1);
    geodesic_smoothness_coef->setSingleStep(0.1);
    QPushButton *geodesicCalcButton = new QPushButton(tr("Calculate"));
    QGridLayout *geodesicLayout = new QGridLayout;
    geodesicLayout->addWidget(new QLabel(tr("distance")), 0, 0);
    geodesicLayout->addWidget(new QLabel(tr("skeleton")), 0, 2);
    geodesicLayout->addWidget(new QLabel(tr("curvature")), 1, 0);
    geodesicLayout->addWidget(new QLabel(tr("smoothness")), 1, 2);
    geodesicLayout->addWidget(geodesic_distance_coef, 0, 1);
    geodesicLayout->addWidget(geodesic_skeleton_coef, 0, 3);
    geodesicLayout->addWidget(geodesic_curvature_coef, 1, 1);
    geodesicLayout->addWidget(geodesic_smoothness_coef, 1, 3);
    geodesicLayout->addWidget(geodesicCalcButton, 2, 0, 1, 4);
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
    
    QGroupBox *visualization_group = new QGroupBox(tr("Visualization"));
    QPushButton *vis_intersection_points_button = new QPushButton(tr("Intersection points"));
    vis_seam_segments_button = new QPushButton(tr("Segments"));
    vis_seam_segments_importance_button = new QPushButton(tr("Segments importance"));
    QPushButton *vis_skeleton_button = new QPushButton(tr("Skeleton"));
    QPushButton *vis_save_home_view_button = new QPushButton(tr("save view"));
    QPushButton *vis_load_home_view_button = new QPushButton(tr("load view"));
    QVBoxLayout *vis_layout = new QVBoxLayout;
    QHBoxLayout *vis_layout_row1 = new QHBoxLayout;
    QHBoxLayout *vis_layout_row2 = new QHBoxLayout;
    vis_layout_row1->addWidget(vis_intersection_points_button);
    vis_layout_row1->addWidget(vis_seam_segments_button);
    vis_layout_row1->addWidget(vis_seam_segments_importance_button);
    vis_layout_row1->addWidget(vis_skeleton_button);
    vis_layout_row2->addWidget(vis_save_home_view_button);
    vis_layout_row2->addWidget(vis_load_home_view_button);
    vis_layout->addLayout(vis_layout_row1);
    vis_layout->addLayout(vis_layout_row2);
    visualization_group->setLayout(vis_layout);

    layout->addWidget(seamGroup);
    layout->addWidget(subsetGroup);
    layout->addWidget(geodesicGroup);
    layout->addWidget(selectionGroup);
    layout->addWidget(flatteningGroup);
    layout->addWidget(curvatureGroup);
    layout->addWidget(skeletonWeightGroup);
    layout->addWidget(visualization_group);
    
    connect(seamShowSingleButton, SIGNAL(clicked()), this, SLOT(seamShowButtonClicked()));
    connect(seamShowAllButton, SIGNAL(clicked()), this, SLOT(seamShowButtonClicked()));
    connect(seamLocalButton, SIGNAL(clicked()), this, SLOT(seamShowButtonClicked()));
    connect(seam_from_selection_button, SIGNAL(clicked()), this, SLOT(seam_from_selection_button_clicked()));
    connect(seamLoadButton, SIGNAL(clicked()), this, SLOT(seamLoadButtonClicked()));
    connect(seamSaveButton, SIGNAL(clicked()), this, SLOT(seamSaveButtonClicked()));
    connect(seamMergeButton, SIGNAL(clicked()), this, SLOT(seamMergeButtonClicked()));
    
    connect(subset_calc_button, SIGNAL(clicked()), this, SLOT(subset_calc_button_clicked()));
    connect(subset_show_button, SIGNAL(clicked()), this, SLOT(subset_show_button_clicked()));

    connect(geodesicCalcButton, SIGNAL(clicked()), this, SLOT(calcGeodesicButtonClicked()));

    connect(loadSelectionButton, SIGNAL(clicked()), this, SLOT(loadSelectionButtonClicked()));
    connect(saveSelectionButton, SIGNAL(clicked()), this, SLOT(saveSelectionButtonClicked()));
    connect(clearSelectionButton, SIGNAL(clicked()), this, SLOT(clearSelectionButtonClicked()));
    connect(calcSelectionButton, SIGNAL(clicked()), this, SLOT(calcSelectionButtonClicked()));

    connect(calcFlattenButton, SIGNAL(clicked()), this, SLOT(calcFlattenedGraphButtonClicked()));
    connect(showFlattenButton, SIGNAL(clicked()), this, SLOT(showFlattenedGrpahButtonClicked()));

    connect(calcCurvatureButton, SIGNAL(clicked()), this, SLOT(calcCurvatureButtonClicked()));
    
    connect(calcSkeletonWeightButton, SIGNAL(clicked()), this, SLOT(calcSkeletonWeightButtonClicked()));

    connect(vis_intersection_points_button, SIGNAL(clicked()), this, SLOT(vis_intersection_points_button_clicked()));
    connect(vis_seam_segments_button, SIGNAL(clicked()), this, SLOT(vis_seam_segments_button_clicked()));
    connect(vis_seam_segments_importance_button, SIGNAL(clicked()), this, SLOT(vis_seam_segments_button_clicked()));
    connect(vis_skeleton_button, SIGNAL(clicked()), this, SLOT(vis_skeleton_button_clicked()));
    connect(vis_save_home_view_button, SIGNAL(clicked()), this, SLOT(vis_save_home_view_button_clicked()));
    connect(vis_load_home_view_button, SIGNAL(clicked()), this, SLOT(vis_load_home_view_button_clicked()));
    
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
            seamNumPaths->setMaximum(selectedCount*(selectedCount-1)/2);
            seamNumPaths->setValue(selectedCount*(selectedCount-1)/2);

            // When user click on "cancel" button, cancel this job
            connect(this, SIGNAL(cancelingJob()), m_patternGenerator, SLOT(cancelJob()));
            // Get progress info from m_patternGenerator, then re-emit it to OpenFlipper system inside receiveJobState
            connect(m_patternGenerator, SIGNAL(setJobState(int)), this, SLOT(receiveJobState(int)));
            // Propagate messages from m_patternGenerator to OpenFlipper
            connect(m_patternGenerator, SIGNAL(log(int, QString)), this, SLOT(receiveLog(int, QString)));
            
            connect(m_patternGenerator, SIGNAL(updateView()), this, SLOT(receiveUpdate()));
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
    MeshSelection::clearEdgeSelection(mesh);
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
    
    m_patternGenerator->calcCurvature();
}

void PlushPlugin::calcGeodesicButtonClicked() {
    if (!checkIfGeneratorExist()) {
        return;
    }
    
    m_patternGenerator->set_geodesic_coeffifients(geodesic_distance_coef->value(),
                                                  0,
                                                  geodesic_curvature_coef->value(),
                                                  geodesic_skeleton_coef->value(),
                                                  geodesic_smoothness_coef->value());

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

void PlushPlugin::seam_from_selection_button_clicked() {
    if (!checkIfGeneratorExist()) {
        return;
    }
    m_patternGenerator->set_seams_from_selection();
}

void PlushPlugin::seamLoadButtonClicked() {
    if (!checkIfGeneratorExist()) {
        return;
    }
    m_patternGenerator->load_seams();
}

void PlushPlugin::seamSaveButtonClicked() {
    if (!checkIfGeneratorExist()) {
        return;
    }
    
    m_patternGenerator->save_seams();
}

void PlushPlugin::subset_calc_button_clicked() {
    if (!checkIfGeneratorExist()) {
        return;
    }
    
    m_triMeshObj->setObjectDrawMode(ACG::SceneGraph::DrawModes::EDGES_COLORED | ACG::SceneGraph::DrawModes::SOLID_FLAT_SHADED);
    m_triMeshObj->materialNode()->enable_alpha_test(0.1);
    if (seamMergeIsStep->isChecked()) {
        calcMergeSegmentThread();
    } else {
        m_currentJobId = "calc subset";
        OpenFlipperThread *thread = new OpenFlipperThread(m_currentJobId);
        connect(thread, SIGNAL(finished(QString)), this, SIGNAL(finishJob(QString)));
        connect(thread, SIGNAL(function(QString)), this, SLOT(calcSubsetThread()), Qt::DirectConnection);
        
        // Custom handler to set m_currentJobId to "" after finishing job.
        connect(thread, SIGNAL(finished(QString)), this, SLOT(finishedJobHandler()));
        
        emit startJob(m_currentJobId, "calc subset", 0, 100, false);
        thread->start();
        thread->startProcessing();
    }
}

void PlushPlugin::seamMergeButtonClicked() {
    if (!checkIfGeneratorExist()) {
        return;
    }
    
    m_triMeshObj->setObjectDrawMode(ACG::SceneGraph::DrawModes::EDGES_COLORED | ACG::SceneGraph::DrawModes::SOLID_FLAT_SHADED);
    m_triMeshObj->materialNode()->enable_alpha_test(0.1);
    if (seamMergeIsStep->isChecked()) {
        calcMergeSegmentThread();
    } else {
        m_currentJobId = "merge";
        OpenFlipperThread *thread = new OpenFlipperThread(m_currentJobId);
        connect(thread, SIGNAL(finished(QString)), this, SIGNAL(finishJob(QString)));
        connect(thread, SIGNAL(function(QString)), this, SLOT(calcMergeSegmentThread()), Qt::DirectConnection);
        
        // Custom handler to set m_currentJobId to "" after finishing job.
        connect(thread, SIGNAL(finished(QString)), this, SLOT(finishedJobHandler()));
        
        emit startJob(m_currentJobId, "merging sub meshes", 0, 100, false);
        thread->start();
        thread->startProcessing();
    }
}

void PlushPlugin::calcFlattenedGraphButtonClicked() {
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

void PlushPlugin::seamShowButtonClicked() {
    if (!checkIfGeneratorExist()) {
        return;
    }
    
    TriMesh *mesh = m_triMeshObj->mesh();

    if (sender() == seamShowSingleButton || sender() == seamShowAllButton) {
        std::vector<int> selectedVerticesId;
        selectedVerticesId = MeshSelection::getVertexSelection(mesh);
        seamNumPaths->setMaximum(selectedVerticesId.size()*(selectedVerticesId.size()-1)/2);
        
        std::vector<VertexHandle> selectedVertices;
        for (size_t i = 0; i < selectedVerticesId.size(); i++) {
            selectedVertices.push_back(mesh->vertex_handle(selectedVerticesId[i]));
        }
        
        if (sender() == seamShowSingleButton) {
            showAllPath = false;
        } else {
            showAllPath = true;
        }
        MeshSelection::clearEdgeSelection(mesh);
        cpu0 = get_cpu_time();
        m_patternGenerator->calcSeams(selectedVertices,
                                      seamMergeThreshold->value(),
                                      seamNumPaths->value(),
                                      seamElimination->isChecked(),
                                      showAllPath);
        cpu1 = get_cpu_time();
    } else if (sender() == seamLocalButton) {
//        cpu0 = get_cpu_time();
//        m_patternGenerator->calcLocalSeams(mesh, seamMergeThreshold->value());
//        cpu1 = get_cpu_time();
//        emit log(LOGINFO, QString("Time: %1").arg(cpu1-cpu0));
    } else {
        return;
    }
    
    std::set<EdgeHandle> seams = m_patternGenerator->get_seams();
    if (!seams.empty()) {
        for (EdgeHandle eh : seams) {
            mesh->status(eh).set_selected(true);
        }
        emit updatedObject(m_triMeshObj->id(), UPDATE_SELECTION);
    }
}

void PlushPlugin::vis_intersection_points_button_clicked() {
    if (!checkIfGeneratorExist()) {
        return;
    }
    
    std::set<VertexHandle> intersection_points = m_patternGenerator->get_intersection_points();
    TriMesh *mesh = m_triMeshObj->mesh();
    for (VertexHandle v : intersection_points) {
        mesh->status(v).set_feature(true);
    }
    emit updatedObject(m_triMeshObj->id(), UPDATE_ALL);
}

void PlushPlugin::vis_seam_segments_button_clicked() {
    if (!checkIfGeneratorExist()) {
        return;
    }
    
    TriMesh *mesh = m_triMeshObj->mesh();
    m_triMeshObj->materialNode()->enable_alpha_test(0.1);
    m_triMeshObj->materialNode()->set_line_smooth(true);
    m_triMeshObj->materialNode()->set_line_width(3);
    
    std::vector< std::vector<HalfedgeHandle> > heh_segments = m_patternGenerator->get_seams_segments();
    
    MeshSelection::clearVertexSelection(mesh);
    MeshSelection::clearEdgeSelection(mesh);

    int num = heh_segments.size();
    for (EdgeHandle eh : mesh->edges()) {
        mesh->set_color(eh, TriMesh::Color(1,1,1,0));
    }
    
    if (sender() == vis_seam_segments_button) {
        // generate unique random sequence
        int *sequence = new int[num];
        for (int i = 0; i < num; i++) {
            sequence[i] = i;
        }
        // randomly swap
        for (int i = 0; i < num; i++) {
            int r = rand() % num;
            int tmp = sequence[i];
            sequence[i] = sequence[r];
            sequence[r] = tmp;
        }
        int counter = 0;
        ACG::ColorCoder color_coder(0, num, false);
        for (auto group : heh_segments) {
            TriMesh::Color color = color_coder.color_float4(sequence[counter++] % num);
            for (HalfedgeHandle heh : group) {
                mesh->set_color(mesh->edge_handle(heh), color);
            }
        }
    } else {
        double max = -1;
        double *group_weight_sum = new double[heh_segments.size()];
        int counter = 0;
        for (auto group : heh_segments) {
            for (HalfedgeHandle heh : group) {
                if (!PlushPatternGenerator::is_different_texture(mesh, mesh->edge_handle(heh))) {
                        group_weight_sum[counter] += mesh->property(PlushPatternGenerator::edgeWeightHandle, mesh->edge_handle(heh));
                }
            }
            if (group_weight_sum[counter] > max) {
                max = group_weight_sum[counter];
            }
            counter++;
        }
        ACG::ColorCoder color_coder(0, log2(max), false);
        counter = 0;
        for (auto group : heh_segments) {
            TriMesh::Color color = color_coder.color_float4(log2(group_weight_sum[counter++]));
            for (HalfedgeHandle heh : group) {
                mesh->set_color(mesh->edge_handle(heh), color);
            }
        }
        delete[] group_weight_sum;
    }
    PluginFunctions::setDrawMode(ACG::SceneGraph::DrawModes::EDGES_COLORED | ACG::SceneGraph::DrawModes::SOLID_FLAT_SHADED);
    emit updatedObject(m_triMeshObj->id(), UPDATE_COLOR | UPDATE_SELECTION);
}

void PlushPlugin::vis_skeleton_button_clicked() {
    // Visualize skeleton
    Skeleton *skeleton = m_triMeshObj->mesh()->property(PlushPatternGenerator::skeletonHandle);
    for (Bone bone : skeleton->bones) {
        int polyLineObjId;
        emit addEmptyObject(DATA_POLY_LINE, polyLineObjId);
        PolyLineObject *object = 0;
        PluginFunctions::getObject(polyLineObjId, object);
        PolyLine *polyLine = object->line();
        polyLine->add_point(bone.getA());
        polyLine->add_point(bone.getB());
    }
    PluginFunctions::setDrawMode(ACG::SceneGraph::DrawModes::WIREFRAME);
    emit updatedObject(m_triMeshObj->id(), UPDATE_ALL);
}

void PlushPlugin::vis_save_home_view_button_clicked() {
    QString view = PluginFunctions::getEncodedExaminerView();

    QString meshName = QFileInfo(m_triMeshObj->name()).baseName();
    // Prepare file for saving data
    QFile file(meshName + "_view.txt");
    file.open(QIODevice::WriteOnly | QIODevice::Text);
    QTextStream out(&file);
    out << view;
    file.close();
}

void PlushPlugin::vis_load_home_view_button_clicked() {
    QString meshName = QFileInfo(m_triMeshObj->name()).baseName();
    QFile file(meshName + "_view.txt");
    file.open(QIODevice::ReadOnly | QIODevice::Text);
    QTextStream in(&file);
    QString view = in.readAll();
    file.close();

    PluginFunctions::setEncodedExaminerView(view);
    emit updateView();
}

void PlushPlugin::showFlattenedGrpahButtonClicked() {
    std::vector<FilteredTriMesh> flattenedSubMeshes = m_patternGenerator->getFlattenedSubMeshes();
    if (flattenedSubMeshes.empty()) {
        emit log(LOGERR, "Calculate flattened meshes first.");
        return;
    }
    
    // 1 means right viewer
    int viewerId = 1;
    // We don't want to show both models at same the time.
    // So we place the flattened model to a very far place so that it will not be shown on original viewer
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

    int n = flattenedSubMeshes.size();
    int *newTriMeshIds = new int[n];

    // Used for displaying sub-mesh number
    QFont font("Arial Black", QFont::Bold);
    font.setPointSize(30);

//    ACG::Vec3d mesh_bbMin, mesh_bbMax;
//    MeshInfo::getBoundingBox(m_triMeshObj->mesh(), mesh_bbMin, mesh_bbMax);
//    double mesh_length = (mesh_bbMax-mesh_bbMin).norm();

//    std::vector<int> seamsId;
    const int nColumns = 4;
    std::vector<TriMesh::Point> minP(n);
    std::vector<TriMesh::Point> maxP(n);
    std::vector<double> columnWidth(nColumns, 0);
    std::vector<double> rowHeight(ceil((double)n/nColumns), 0);
    double width = 0, height = 0, depth = 0;
    for (int i = 0; i < n; i++) {
        // Add an empty object for creating new sub-mesh from FilteredTriMesh
        emit addEmptyObject(DATA_TRIANGLE_MESH, newTriMeshIds[i]);
        TriMeshObject *object = 0;
        PluginFunctions::getObject(newTriMeshIds[i], object);
        object->materialNode()->set_line_smooth(true);
        object->materialNode()->enable_alpha_test(0.1);
        object->manipulatorNode()->translate(far);

        TriMesh *mesh = object->mesh();
        mesh->request_edge_colors();
        
        m_patternGenerator->get_triMesh_from_subMesh(mesh, flattenedSubMeshes.at(i), true);

        MeshInfo::getBoundingBox(mesh, minP[i], maxP[i]);
        columnWidth[i % nColumns] = max(columnWidth[i % nColumns], maxP[i][0] - minP[i][0]);
        rowHeight[i / nColumns] = max(rowHeight[i / nColumns], maxP[i][1] - minP[i][1]);
        width = max(width, maxP[i][0] - minP[i][0]);
        height = max(height, maxP[i][1] - minP[i][1]);
        depth = max(depth, maxP[i][2] - minP[i][2]);

        // numbers on sub meshes
//        ACG::SceneGraph::TransformNode *sub_mesh_text_transform_node;
//        std::stringstream sub_mesh_text_transform_node_name;
//        sub_mesh_text_transform_node_name << "sub_mesh_text_transform_" << i;
//        if (!object->getAdditionalNode(sub_mesh_text_transform_node, "Plush", sub_mesh_text_transform_node_name.str().c_str())) {
//            sub_mesh_text_transform_node = new ACG::SceneGraph::TransformNode(object->manipulatorNode(),
//                                                                   sub_mesh_text_transform_node_name.str());
//            object->addAdditionalNode(sub_mesh_text_transform_node, "Plush", sub_mesh_text_transform_node_name.str().c_str());
//        }
//        ACG::Vec3d sub_mesh_min, sub_mesh_max;
//        MeshInfo::getBoundingBox(subMesh, sub_mesh_min, sub_mesh_max);
//        sub_mesh_text_transform_node->loadIdentity();
//        sub_mesh_text_transform_node->translate(sub_mesh_min);
//        
//        ACG::SceneGraph::TextNode *sub_mesh_text_node;
//        std::stringstream sub_mesh_text_node_name;
//        sub_mesh_text_node_name << "sub_mesh_patch_number_" << i;
//        if (!object->getAdditionalNode(sub_mesh_text_node, "Plush", sub_mesh_text_node_name.str().c_str())) {
//            sub_mesh_text_node = new ACG::SceneGraph::TextNode(sub_mesh_text_transform_node,
//                                                    "<TextNode>",
//                                                    ACG::SceneGraph::TextNode::SCREEN_ALIGNED,
//                                                    false);
//            object->addAdditionalNode(sub_mesh_text_node, "Plush", sub_mesh_text_node_name.str().c_str());
//        }
//        
//        sub_mesh_text_node->setFont(font);
//        sub_mesh_text_node->setText(std::to_string(i+1));
//        sub_mesh_text_node->setSize(0.5);
//
//        // numbers on original meshes
//        ACG::SceneGraph::TransformNode *text_transform_node;
//        std::stringstream text_transform_node_name;
//        text_transform_node_name << "text_transform_" << i;
//        if (!m_triMeshObj->getAdditionalNode(text_transform_node, "Plush", text_transform_node_name.str().c_str())) {
//            text_transform_node = new ACG::SceneGraph::TransformNode(m_triMeshObj->manipulatorNode(),
//                                                                   text_transform_node_name.str());
//            m_triMeshObj->addAdditionalNode(text_transform_node, "Plush" , text_transform_node_name.str().c_str());
//        }
//        ACG::Vec3d avg_point(0,0,0);
//        ACG::Vec3d avg_normal(0,0,0);
//        for (VertexHandle v : oldSubMesh.vertices()) {
//            VertexHandle original_v = m_patternGenerator->get_original_handle(&oldSubMesh, v);
//            avg_point += m_triMeshObj->mesh()->point(original_v);
//            avg_normal += m_triMeshObj->mesh()->normal(original_v);
//        }
//        avg_point /= subMesh->n_vertices();
////        avg_normal[2] = 0;
//        avg_normal.normalize();
//
//        double min_distance = std::numeric_limits<double>::max();
//        ACG::Vec3d closed_point(0,0,0);
//        for (VertexHandle v : oldSubMesh.vertices()) {
//            VertexHandle original_v = m_patternGenerator->get_original_handle(&oldSubMesh, v);
//            ACG::Vec3d p = m_triMeshObj->mesh()->point(original_v);
//            if ((p-avg_point).norm() < min_distance) {
//                min_distance = (p-avg_point).norm();
//                closed_point = p;
//            }
//        }
//        
//        text_transform_node->loadIdentity();
//        if (avg_normal[0] < 0 || avg_normal[1] < 0) {
//            text_transform_node->translate(closed_point + avg_normal * mesh_length * 0.01 * (abs(avg_normal[0])+abs(avg_normal[1])));
//        } else {
//            text_transform_node->translate(closed_point);
//        }
//        
//        ACG::SceneGraph::TextNode *text_node;
//        std::stringstream text_node_name;
//        text_node_name << "patch_number_" << i;
//        if (!m_triMeshObj->getAdditionalNode(text_node, "Plush", text_node_name.str().c_str())) {
//            text_node = new ACG::SceneGraph::TextNode(text_transform_node,
//                                                     "<TextNode>",
//                                                     ACG::SceneGraph::TextNode::SCREEN_ALIGNED_STATIC_SIZE,
//                                                     false);
//            m_triMeshObj->addAdditionalNode(text_node, "Plush", text_node_name.str().c_str());
//        }
//        
//        text_node->setFont(font);
//        text_node->setText(std::to_string(i+1));
//        text_node->setPixelSize(20);
    }
    
    // Packing
    double spacing = max(height, width)/5;
    // The anchor is at bottom-left corner of a mesh
    OpenMesh::Vec3d anchor(0,0,0);
    
    for (int i = 0; i < n; i++) {
        TriMeshObject *object = 0;
        PluginFunctions::getObject(newTriMeshIds[i], object);
        TriMesh *mesh = object->mesh();
        
        // For a new row, calculate the bottom-left anchor position
        if (i % nColumns == 0) {
            if (i > 0) {
                // Add row height of previous row
                anchor[1] += rowHeight[i / nColumns - 1] + spacing;
            }
            anchor[0] = 0;
        } else {
            anchor[0] += columnWidth[(i-1) % nColumns] + spacing;
        }
        
        // Move each mesh to the grid center it belongs to
        OpenMesh::Vec3d additional_offset(columnWidth[i % nColumns]/2 - (maxP[i][0] + minP[i][0])/2,
                                          rowHeight[i / nColumns]/2 - (maxP[i][1] + minP[i][1])/2,
                                          -depth/2);
        object->manipulatorNode()->translate(anchor + additional_offset);
    }
//    MeshSelection::selectEdges(m_triMeshObj->mesh(), seamsId);
    PluginFunctions::setDrawMode(ACG::SceneGraph::DrawModes::EDGES_COLORED);
    emit updatedObject(m_triMeshObj->id(), UPDATE_SELECTION);
}

void PlushPlugin::subset_show_button_clicked() {
    auto hierarchical_patches = m_patternGenerator->get_hierarchical_patches();
    if (hierarchical_patches.size() > 0) {
        auto hierarchical_patches_layer = hierarchical_patches[hierarchical_patches.size()-1];
        int counter = 0;
        for (FilteredTriMesh subMesh : hierarchical_patches_layer) {
            int id;
            emit addEmptyObject(DATA_TRIANGLE_MESH, id);
            TriMeshObject *object = 0;
            PluginFunctions::getObject(id, object);
            
            TriMesh *mesh = object->mesh();
            
            m_patternGenerator->get_triMesh_from_subMesh(mesh, subMesh, false);
            MeshSelection::selectBoundaryEdges(mesh);
            counter++;
            if (counter > 5) {
                break;
            }
        }
    }
    
    emit updatedObject(m_triMeshObj->id(), UPDATE_SELECTION);
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

    seamNumPaths->setMaximum(selectedCount*(selectedCount-1)/2);
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

void PlushPlugin::calcGeodesicThread() {
    if (!checkIfGeneratorExist()) {
        return;
    }
    
    TriMesh *mesh = m_triMeshObj->mesh();
    
    std::vector<int> selectedVerticesId;
    selectedVerticesId = MeshSelection::getVertexSelection(mesh);
    seamNumPaths->setMaximum(selectedVerticesId.size()*(selectedVerticesId.size()-1)/2);
    
    std::vector<VertexHandle> selectedVertices;
    for (size_t i = 0; i < selectedVerticesId.size(); i++) {
        selectedVertices.push_back(mesh->vertex_handle(selectedVerticesId[i]));
    }
    cpu0 = get_cpu_time();
    m_patternGenerator->calcGeodesic(selectedVertices);
    cpu1 = get_cpu_time();
    emit log(LOGINFO, QString("Time: %1").arg(cpu1-cpu0));

    m_patternGenerator->saveGeodesic(selectedVertices);

    emit log(LOGINFO, "Geodesic calculation done.");
}

void PlushPlugin::calcFlattenedGraphThread() {
    m_patternGenerator->calcFlattenedGraph();
}

void PlushPlugin::calcMergeSegmentThread() {
//    MeshSelection::clearEdgeSelection(m_triMeshObj->mesh());
    cpu0 = get_cpu_time();
    m_patternGenerator->optimize_patches(seamMergeThreshold->value(), seamMergeIsStep->isChecked());
    cpu1 = get_cpu_time();
    emit log(LOGINFO, QString("Time: %1").arg(cpu1-cpu0));
}


void PlushPlugin::calcSubsetThread() {
    m_patternGenerator->construct_subsets(seamMergeThreshold->value());
//    emit updatedObject(m_triMeshObj->id(), UPDATE_SELECTION);
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

void PlushPlugin::receiveUpdate() {
    emit updatedObject(m_triMeshObj->id(), UPDATE_SELECTION | UPDATE_COLOR);
}

void PlushPlugin::slotKeyEvent( QKeyEvent* _event ) {
    // Switch pressed keys
    int currentEdgeNo = seamNumPaths->value();
    switch (_event->key())
    {
        case Qt::Key_X:
            if (currentEdgeNo < seamNumPaths->maximum()) {
                seamNumPaths->setValue(currentEdgeNo+1);
            }
            seamShowSingleButton->click();
            break;
        case Qt::Key_Z:
            if (currentEdgeNo > 0) {
                seamNumPaths->setValue(currentEdgeNo-1);
            }
            seamShowSingleButton->click();
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
