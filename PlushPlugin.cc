#include "PlushPlugin.hh"
#include "GeodesicDistance/geodesic_algorithm_dijkstra.hh"
#include <ObjectTypes/TriangleMesh/TriangleMesh.hh>
#include <ObjectTypes/PolyMesh/PolyMesh.hh>
#include <OpenFlipper/BasePlugin/PluginFunctions.hh>
#include <OpenFlipper/BasePlugin/RPCWrappers.hh>

PlushPlugin::PlushPlugin()
{
}
void PlushPlugin::initializePlugin()
{
    filename = "";
    sourceIdx = -1;
    destIdx = -1;
    
    requiredPlugins.push_back("meshobjectselection");
    
    
    // Register keys
    emit registerKey(Qt::Key_X,     Qt::NoModifier, "Normal mode");
    emit registerKey(Qt::Key_Z,     Qt::NoModifier, "Plush picking mode");
}

void PlushPlugin::pluginsInitialized() {
    
    
    // Create the Toolbox Widget
    QWidget* toolBox = new QWidget();
    QGridLayout* layout = new QGridLayout(toolBox);
    
    // Create button that can be toggled
    // to (de)activate plugin's picking mode
    QPushButton *geodesicButton = new QPushButton(tr("calculate"));
    QLabel *label = new QLabel("Geodesic distance");
    layout->addWidget(label, 0, 0);
    layout->addWidget(geodesicButton, 0, 1);
    connect(geodesicButton, SIGNAL(clicked()), this, SLOT(calcGeodesic()));
    
    emit addToolbox(tr("Plush"), toolBox);

//    emit load("../../../box.obj", DATA_TRIANGLE_MESH, objId);
    
//    // Add context menu entry/Users/MBRSL/Downloads/geodesic_cpp_03_02_2008/readme.txt
//    // Create submenu
//    contextMenuEntry_ = new QMenu("Plush plugin");
//    QAction* lastAction;
//    // Add action to recently created menu
//    lastAction = contextMenuEntry_->addAction( "Hide object" );
//    lastAction->setToolTip("Hide selected object");
//    lastAction->setStatusTip( lastAction->toolTip() );
//    // Finally insert created context submenu to OpenFlipper's context menu
//    // We want our action to be visible for triangle and polygon meshes
//    emit addContextMenuItem(contextMenuEntry_->menuAction() , DATA_TRIANGLE_MESH , CONTEXTOBJECTMENU );
//    //    emit addContextMenuItem(contextMenuEntry_->menuAction() , DATA_POLY_MESH , CONTEXTOBJECTMENU );
//    // Connect the created context menu entry to local function contextMenuItemSelected(QAction*)
//    connect(contextMenuEntry_, SIGNAL(triggered(QAction*)), this, SLOT(contextMenuItemSelected(QAction*)));
}
void PlushPlugin::showGeodesic() {
    if (PluginFunctions::objectCount() == 0)
    {
        emit log(LOGERR, "Load a model first.");
        return;
    }
    
    for (PluginFunctions::ObjectIterator o_it(PluginFunctions::TARGET_OBJECTS); o_it != PluginFunctions::objectsEnd();
         ++o_it) {
        // pick two vertices
        TriMesh *mesh = NULL;

        QString filename;
        int meshId = -1;
        int sourceIdx = -1;
        int destIdx = -1;

        if (o_it->dataType(DATA_TRIANGLE_MESH)) {
            meshId = o_it->id();
            mesh = PluginFunctions::triMesh(*o_it);
            filename = o_it->path() + OpenFlipper::Options::dirSeparator() + o_it->filename();
            
            IdList selectedVertices = RPC::callFunctionValue<IdList> ("meshobjectselection", "getVertexSelection", o_it->id());
            if (selectedVertices.size() == 2) {
                sourceIdx = selectedVertices[0];
                destIdx = selectedVertices[1];
                TriMesh::Point p1 = mesh->point(mesh->vertex_handle(sourceIdx));
                TriMesh::Point p2 = mesh->point(mesh->vertex_handle(destIdx));
                char strS[100], strD[100];
                sprintf(strS, "source: (%lf, %lf, %lf)", p1[0], p1[1], p1[2]);
                sprintf(strD, "dest: (%lf, %lf, %lf)", p2[0], p2[1], p2[2]);
                emit log(LOGINFO, strS);
                emit log(LOGINFO, strD);
            }
            break;
        }

        if (sourceIdx == -1 || destIdx == -1) {
            emit log(LOGERR, "Select exact two points first.");
            return;
        }
    }
}
void PlushPlugin::calcGeodesic() {
    if (PluginFunctions::objectCount() == 0)
    {
        emit log(LOGERR, "Load a model first.");
        return;
    }
    
    TriMesh *mesh = NULL;
    int meshId = -1;

    for (PluginFunctions::ObjectIterator o_it(PluginFunctions::TARGET_OBJECTS); o_it != PluginFunctions::objectsEnd();
         ++o_it) {
        if (o_it->dataType(DATA_TRIANGLE_MESH)) {
            meshId = o_it->id();
            mesh = PluginFunctions::triMesh(*o_it);

            // write each connected component to .obj
            std::ifstream off("tmp.obj", std::ios::out);
            for (VertexIter v_it = mesh->vertices_begin(); v_it != mesh->vertices_end(); v_it++)
            {
                TriMesh::Point p = mesh->point(*v_it);
                off << "v " << p[0] << " " << p[1] << " " << p[2];
            }
            
            std::vector<double> points;
            std::vector<unsigned> faces;
            
    bool success = geodesic::read_mesh_from_file(filename,points,faces);
    if(!success)
    {
        emit log(LOGERR, "Error loading file in calcGeodesic.");
        return;
    }
    
    geodesic::Mesh gMesh;
    gMesh.initialize_mesh_data(points, faces);		//create internal mesh data structure including edges
    
//    geodesic::GeodesicAlgorithmExact algorithm(&gMesh);	//create exact algorithm for the mesh
    geodesic::GeodesicAlgorithmDijkstra algorithm(&gMesh);	//create exact algorithm for the mesh
    
    geodesic::SurfacePoint source(&gMesh.vertices()[sourceIdx]);		//create source
    std::vector<geodesic::SurfacePoint> all_sources(1,source);					//in general, there could be multiple sources, but now we have only one
    
    if(destIdx != -1)	//target vertex specified, compute single path
    {
        geodesic::SurfacePoint target(&gMesh.vertices()[destIdx]);		//create source
        
        std::vector<geodesic::SurfacePoint> path;	//geodesic path is a sequence of SurfacePoints
        
        bool const lazy_people_flag = false;		//there are two ways to do exactly the same
        if(lazy_people_flag)
        {
            algorithm.geodesic(source, target, path); //find a single source-target path
        }
        else		//doing the same thing explicitly for educational reasons
        {
            double const distance_limit = geodesic::GEODESIC_INF;			// no limit for propagation
            std::vector<geodesic::SurfacePoint> stop_points(1, target);	//stop propagation when the target is covered
            algorithm.propagate(all_sources, distance_limit, &stop_points);	//"propagate(all_sources)" is also fine, but take more time because covers the whole mesh
            
            algorithm.trace_back(target, path);		//trace back a single path
        }
        
        geodesic::print_info_about_path(path);
        for(unsigned i = 0; i < path.size()-1; ++i)
        {
            geodesic::SurfacePoint& v1 = path[i];
            geodesic::SurfacePoint& v2 = path[i+1];
            OpenMesh::Vec3d p1 = OpenMesh::Vec3d(v1.x(), v1.y(), v1.z());
            OpenMesh::Vec3d p2 = OpenMesh::Vec3d(v2.x(), v2.y(), v2.z());
            OpenMesh::Vec3d pMiddle = OpenMesh::Vec3d((v1.x()+v2.x())/2, (v1.y()+v2.y())/2, (v1.z()+v2.z())/2);
            char str[100];
            
            OpenMesh::EdgeHandle _eh;
            if (!pickEdge(mesh, _eh, p1, p2))
            {
                sprintf(str, "edge center at (%lf %lf %lf) not found", pMiddle[0], pMiddle[1], pMiddle[2]);
                emit log(LOGERR, str);
            }
            else
            {
                IdList list;
                list.push_back(_eh.idx());
                RPC::callFunction<int, IdList>("meshobjectselection", "selectEdges", meshId, list);
                sprintf(str, "edge at (%lf %lf %lf)", pMiddle[0], pMiddle[1], pMiddle[2]);
                emit log(LOGINFO, str);
            }
        }
    }
    /*
    else		//target vertex is not specified, print distances to all vertices
    {
        algorithm.propagate(all_sources);	//cover the whole mesh
        
        for(unsigned i=0; i<gMesh.vertices().size(); ++i)
        {
            geodesic::SurfacePoint p(&gMesh.vertices()[i]);
            
            double distance;
            unsigned best_source = algorithm.best_source(p,distance);		//for a given surface point, find closets source and distance to this source
            
            char str[100];
            sprintf(str, "%lf ", distance);
            emit log(LOGINFO, str);
        }
    }
     */
    
    return;
}

bool PlushPlugin::pickEdge(TriMesh *mesh, OpenMesh::EdgeHandle &_eh, OpenMesh::Vec3d p1, OpenMesh::Vec3d p2) {
    double minDistanceSqr = -1;
    OpenMesh::VertexHandle v1, v2;
    for (TriMesh::VertexIter v_it = mesh->vertices_begin(); v_it != mesh->vertices_end(); v_it++)
    {
        OpenMesh::Vec3d p = mesh->point(*v_it);
        if ((p-p1).sqrnorm() < minDistanceSqr || minDistanceSqr < 0) {
            minDistanceSqr = (p-p1).sqrnorm();
            // we found p1, find p1-p2 from neighbors
            v1 = *v_it;
            
            double minDistanceSqr2 = -1;
            for (TriMesh::VertexEdgeIter ve_it = mesh->ve_iter(v1); ve_it; ve_it++)
            {
                TriMesh::HalfedgeHandle heh = mesh->halfedge_handle(*ve_it, 0);
                if (mesh->from_vertex_handle(heh) == v1)
                {
                    v2 = mesh->to_vertex_handle(heh);
                }
                else
                {
                    v2 = mesh->from_vertex_handle(heh);
                }
                
                OpenMesh::Vec3d pp = mesh->point(v2);
                if ((pp-p2).sqrnorm() < minDistanceSqr2 || minDistanceSqr2 < 0) {
                    minDistanceSqr2 = (pp-p2).sqrnorm();
                    // we found (possible) p2
                    _eh = *ve_it;
                }
            }
        }
    }
    if (minDistanceSqr < 0)
    {
        return false;
    }
    else
    {
        return true;
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

Q_EXPORT_PLUGIN2( plushPlugin , PlushPlugin );