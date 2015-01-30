#include <ObjectTypes/Coordsys/Coordsys.hh>
#include <ObjectTypes/PolyMesh/PolyMesh.hh>
#include <ObjectTypes/PolyLine/PolyLine.hh>
#include <ObjectTypes/TriangleMesh/TriangleMesh.hh>

#include <OpenFlipper/BasePlugin/RPCWrappers.hh>
#include <OpenFlipper/common/Types.hh>

#include "PlushPlugin.hh"
#include "GeodesicDistance/geodesic_algorithm_dijkstra.hh"

PlushPlugin::PlushPlugin()
{
    requiredPlugins = new std::vector<char*>();
    
    geodesicDistance = new std::map<std::pair<VertexHandle, VertexHandle>, double>();
    geodesicPath = new std::map<std::pair<VertexHandle, VertexHandle>, IdList>();
    
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
    
    delete geodesicDistance;
    delete geodesicPath;
    
    delete geodesicEdges;
    delete geodesicButton;
    delete geodesicEdges;
    delete ridgeButton;
    delete loadCurvatureButton;
    delete showCurvatureButton;
    delete calcCurvatureButton;

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
    geodesicEdges = new QSpinBox();
    geodesicButton = new QPushButton(tr("Show path"));
    ridgeButton = new QPushButton(tr("Show ridge"));
    loadCurvatureButton = new QPushButton(tr("Load curvature from file"));
    showCurvatureButton = new QPushButton(tr("Show curvature"));
    calcCurvatureButton = new QPushButton(tr("Calculate curvature"));
    layout->addWidget(geodesicButton, 0, 0);
    layout->addWidget(geodesicEdges, 0, 1);
    layout->addWidget(geodesicButton, 0, 0);
    layout->addWidget(ridgeButton, 1, 0);
    layout->addWidget(loadCurvatureButton, 2, 0);
    layout->addWidget(showCurvatureButton, 3, 0);
    layout->addWidget(calcCurvatureButton, 4, 0);
    connect(geodesicButton, SIGNAL(clicked()), this, SLOT(showGeodesic()));
    connect(ridgeButton, SIGNAL(clicked()), this, SLOT(showRidge()));
    connect(loadCurvatureButton, SIGNAL(clicked()), this, SLOT(loadCurvatureButtonClicked()));
    connect(showCurvatureButton, SIGNAL(clicked()), this, SLOT(showCurvatureButtonClicked()));
    connect(calcCurvatureButton, SIGNAL(clicked()), this, SLOT(calcCurvatureButtonClicked()));
    
    emit addToolbox(tr("Plush"), toolBox);
}

void PlushPlugin::showRidge() {
    // parse file
    QFile input("/Users/MBRSL/Dropbox/Curricular/Plush project/OpenFlipper/Plugin-Plush/ridge.txt");
    if (input.open(QIODevice::ReadOnly)) {
        QTextStream in(&input);
        // each line is one polyline
        while (!in.atEnd()) {
            QString line = in.readLine();
            
            QTextStream textStream(&line);
            int lineType;
            double strength, sharpness;
            textStream >> lineType >> strength >> sharpness;
            
            if (strength >= 5) {
                int objId;
                emit addEmptyObject(DATA_POLY_LINE, objId);
                PolyLineObject *object = 0;
                PluginFunctions::getObject(objId, object);
                PolyLine *polyLine = object->line();

                double x,y,z;
                while (!(textStream >> x >> y >> z).atEnd()) {
                    polyLine->add_point(PolyLine::Point(x,y,z));
                }
            }
        }
        input.close();
    }
}

void PlushPlugin::calcCurvatureButtonClicked() {
    if (thread == NULL) {
        thread = new OpenFlipperThread("calcCurvature");
        connect(thread, SIGNAL(finished(QString)), this, SLOT(finishedJob(QString)));
        connect(thread, SIGNAL(function(QString)), this, SLOT(calcCurvature(QString)), Qt::DirectConnection);
        emit startJob("calcCurvature", "calculate curvature of whole mesh", 0, 100, true);
        thread->start();
        thread->startProcessing();
    }
}

void PlushPlugin::loadCurvatureButtonClicked() {
    for (PluginFunctions::ObjectIterator o_it(PluginFunctions::TARGET_OBJECTS); o_it != PluginFunctions::objectsEnd();
         ++o_it) {
        int meshId = o_it->id();
        
        if (o_it->dataType(DATA_TRIANGLE_MESH)) {
            TriMesh *mesh = PluginFunctions::triMesh(*o_it);

            QFile file(o_it->name()+"_curvature.txt");
            if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
                emit log(LOGERR, QString("Unable to read file %1").arg(o_it->name()+"_curvature.txt"));
                continue;
            }
            
            QTextStream in(&file);

            // Add property so that we can store it into TriMesh
            mesh->add_property(maxCurvatureHandle, "Max Curvature");
            mesh->add_property(minCurvatureHandle, "Min Curvature");
            mesh->add_property(maxCurvatureDirectionHandle, "Max curvature direction");
            mesh->add_property(minCurvatureDirectionHandle, "Min curvature direction");

            double minK1 = 1e9, minK2 = 1e9;
            double d1x, d1y, d1z, d2x, d2y, d2z, k1, k2;
            for (VertexIter v_it = mesh->vertices_begin(); v_it != mesh->vertices_end(); v_it++) {
                in >> d1x >> d1y >> d1z >> d2x >> d2y >> d2z >> k1 >> k2;
                mesh->property(maxCurvatureDirectionHandle, *v_it) = OpenMesh::Vec3d(d1x, d1y, d1z);
                mesh->property(minCurvatureDirectionHandle, *v_it) = OpenMesh::Vec3d(d2x, d2y, d2z);
                mesh->property(maxCurvatureHandle, *v_it) = k1;
                mesh->property(minCurvatureHandle, *v_it) = k2;
            }
        }
    }
    emit log(LOGINFO, "Load curvature complete.");
}

void PlushPlugin::showCurvatureButtonClicked() {
    for (PluginFunctions::ObjectIterator o_it(PluginFunctions::TARGET_OBJECTS); o_it != PluginFunctions::objectsEnd();
         ++o_it) {
        int meshId = o_it->id();
        
        if (o_it->dataType(DATA_TRIANGLE_MESH)) {
            TriMesh *mesh = PluginFunctions::triMesh(*o_it);
            
            QString meshName = o_it->name();
            
            ACG::Vec3d boundingBoxMin, boundingBoxMax;
            o_it->getBoundingBox(boundingBoxMin, boundingBoxMax);
            double meshSize = (boundingBoxMax-boundingBoxMin).length();

            for (VertexIter v_it = mesh->vertices_begin(); v_it != mesh->vertices_end(); v_it++) {
                // show curvature direction
            }
        }
    }
}

void PlushPlugin::showGeodesic() {
    if (PluginFunctions::objectCount() == 0)
    {
        emit log(LOGERR, "Load a model first.");
        return;
    }

    for (PluginFunctions::ObjectIterator o_it(PluginFunctions::TARGET_OBJECTS); o_it != PluginFunctions::objectsEnd();
         ++o_it) {
        if (o_it->dataType(DATA_TRIANGLE_MESH)) {
            int meshId = o_it->id();
            QString meshName = o_it->name();
            TriMesh *mesh = PluginFunctions::triMesh(*o_it);
            
            IdList selectedVertices;
            selectedVertices = RPC::callFunctionValue<IdList> ("meshobjectselection", "getVertexSelection", meshId);
            // if the selection is not empty, save it as default selection for next time.
//            if (selectedVertices.size() > 0) {
//                RPC::callFunction<int, QString>("meshobjectselection", "saveFlipperModelingSelection", meshId, meshName + "_selection.txt");
//            } else {
//                RPC::callFunction<int, QString>("meshobjectselection", "loadFlipperModelingSelection", meshId, meshName + "_selection.txt");
//                selectedVertices = RPC::callFunctionValue<IdList> ("meshobjectselection", "getVertexSelection", meshId);
//            }

            std::set<EdgeHandle> spanningTree;
            calcSpanningTree(mesh, meshId, spanningTree, selectedVertices);
            
            IdList edgeList;
            for (std::set<EdgeHandle>::iterator e_it = spanningTree.begin(); e_it != spanningTree.end(); e_it++) {
                edgeList.push_back(e_it->idx());
            }
            RPC::callFunction<int, IdList>("meshobjectselection", "selectEdges", meshId, edgeList);
        }
    }
}

void PlushPlugin::getOrderedSelectedVertices(TriMesh *mesh, int meshId, IdList *selectedVertices) {
    *selectedVertices = RPC::callFunctionValue<IdList> ("meshobjectselection", "getVertexSelection", meshId);
    
    if (selectedVertices->size() < 2) {
        QString str = QString("Object %1 contains less than 2 vertices. no path showed.").arg(meshId);
        emit log(LOGINFO, str);
        return;
    }
    
    // select vertex with min Z coordinate
    double minZ = 1e9;
    size_t minZ_idx = 0;
    for (size_t sortI = 0; sortI < selectedVertices->size(); sortI++)
    {
        TriMesh::Point p1 = mesh->point(mesh->vertex_handle((*selectedVertices)[sortI]));
        if (p1[2] < minZ)
        {
            minZ = p1[2];
            minZ_idx = sortI;
        }
    }
    // select shortest path greedly
    for (size_t sortI = 0; sortI < selectedVertices->size()-1; sortI++)
    {
        TriMesh::Point p1 = mesh->point(mesh->vertex_handle((*selectedVertices)[sortI]));
        double minDist = 1e9;
        size_t min_idx = sortI+1;
        for (size_t sortJ = sortI+1; sortJ < selectedVertices->size(); sortJ++)
        {
            TriMesh::Point p2 = mesh->point(mesh->vertex_handle((*selectedVertices)[sortJ]));
            if ((p1-p2).length() < minDist)
            {
                minDist = (p1-p2).length();
                min_idx = sortJ;
            }
        }
        
        // exchange position
        size_t tmpId = (*selectedVertices)[sortI+1];
        (*selectedVertices)[sortI+1] = (*selectedVertices)[min_idx];
        (*selectedVertices)[min_idx] = tmpId;
    }
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