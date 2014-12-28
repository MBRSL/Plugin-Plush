#include "PlushPlugin.hh"
#include "GeodesicDistance/geodesic_algorithm_dijkstra.hh"
#include <ObjectTypes/Coordsys/Coordsys.hh>
#include <ObjectTypes/TriangleMesh/TriangleMesh.hh>
#include <ObjectTypes/PolyMesh/PolyMesh.hh>
#include <ObjectTypes/PolyLine/PolyLine.hh>
#include <OpenFlipper/BasePlugin/PluginFunctions.hh>
#include <OpenFlipper/BasePlugin/RPCWrappers.hh>

// Curvature
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Monge_via_jet_fitting.h>
#include <CGAL/Eigen_svd.h>
typedef double                   DFT;
typedef CGAL::Simple_cartesian<DFT>     Data_Kernel;
typedef CGAL::Simple_cartesian<DFT>     Local_Kernel;
typedef CGAL::Eigen_svd SvdTraits;
typedef Data_Kernel::Point_3     DPoint;
typedef CGAL::Monge_via_jet_fitting<Data_Kernel,Local_Kernel,SvdTraits> My_Monge_via_jet_fitting;
typedef My_Monge_via_jet_fitting::Monge_form     My_Monge_form;


PlushPlugin::PlushPlugin()
{
}

void PlushPlugin::initializePlugin()
{
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
    QPushButton *geodesicButton = new QPushButton(tr("Show path"));
    QPushButton *ridgeButton = new QPushButton(tr("Show ridge"));
    QPushButton *curvatureButton = new QPushButton(tr("Show curvature"));
    QPushButton *allCurvatureButton = new QPushButton(tr("Show all curvature"));
    layout->addWidget(geodesicButton, 0, 0);
    layout->addWidget(ridgeButton, 0, 1);
    layout->addWidget(curvatureButton, 1, 0);
    layout->addWidget(allCurvatureButton, 1, 1);
    connect(geodesicButton, SIGNAL(clicked()), this, SLOT(showGeodesic()));
    connect(ridgeButton, SIGNAL(clicked()), this, SLOT(showRidge()));
    connect(curvatureButton, SIGNAL(clicked()), this, SLOT(showSelectedCurvature()));
    connect(allCurvatureButton, SIGNAL(clicked()), this, SLOT(showAllCurvature()));
    
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

void PlushPlugin::showSelectedCurvature() {
    for (PluginFunctions::ObjectIterator o_it(PluginFunctions::TARGET_OBJECTS); o_it != PluginFunctions::objectsEnd();
         ++o_it) {
        int meshId = o_it->id();
        
        // clear previous coordinate gizmos
        if (o_it->dataType(DATA_COORDSYS)) {
            emit deleteObject(meshId);
        }
        if (o_it->dataType(DATA_TRIANGLE_MESH)) {
            TriMesh *mesh = PluginFunctions::triMesh(*o_it);
            
            QString meshName = o_it->name();
            
            ACG::Vec3d boundingBoxMin, boundingBoxMax;
            o_it->getBoundingBox(boundingBoxMin, boundingBoxMax);
            double meshSize = (boundingBoxMax-boundingBoxMin).length();
            
            IdList selectedVertices = RPC::callFunctionValue<IdList> ("meshobjectselection", "getVertexSelection", meshId);
            
            showCurvature(mesh, meshName, meshId,
                          selectedVertices, meshSize);
        }
    }
}
void PlushPlugin::showAllCurvature() {
    // 0.5 = 50% of total points
    double sampleAmount = 0.01;
    for (PluginFunctions::ObjectIterator o_it(PluginFunctions::TARGET_OBJECTS); o_it != PluginFunctions::objectsEnd();
         ++o_it) {
        int meshId = o_it->id();
        
        // clear previous coordinate gizmos
        if (o_it->dataType(DATA_COORDSYS)) {
            emit deleteObject(meshId);
        }
        if (o_it->dataType(DATA_TRIANGLE_MESH)) {
            TriMesh *mesh = PluginFunctions::triMesh(*o_it);
            
            QString meshName = o_it->name();
            
            ACG::Vec3d boundingBoxMin, boundingBoxMax;
            o_it->getBoundingBox(boundingBoxMin, boundingBoxMax);
            double meshSize = (boundingBoxMax-boundingBoxMin).length();
            
            IdList selectedVertices;
            int count = 0;
            for (TriMesh::VertexIter v_it = mesh->vertices_begin(); v_it != mesh->vertices_end(); v_it++, count++) {
                if (count % (int)(1/sampleAmount) == 0) {
                    selectedVertices.push_back(v_it->idx());
                }
            }

            showCurvature(mesh, meshName, meshId,
                          selectedVertices, meshSize);
        }
    }
}

void PlushPlugin::showCurvature(TriMesh *mesh, QString meshName, int meshId,
                                IdList selectedVertices, double meshSize) {
    size_t d_fitting = 3;
    size_t d_monge = 3;
    double indicatorSize = meshSize/30;
    
    // store original selction
    IdList originalVerticesSelection = RPC::callFunctionValue<IdList> ("meshobjectselection", "getVertexSelection", meshId);
    
    // pick vertices
    for (size_t i = 0; i < selectedVertices.size(); i++) {
        IdList currentVertex = IdList();
        currentVertex.push_back(selectedVertices[i]);

        // select only one vertex as target vertex
        RPC::callFunctionValue<IdList> ("meshobjectselection", "clearVertexSelection", meshId);
        RPC::callFunctionValue<IdList> ("meshobjectselection", "selectVertices", meshId, currentVertex);
        
        // get 3-ring vertices around target vertex
        RPC::callFunctionValue<IdList> ("meshobjectselection", "growVertexSelection", meshId);
        RPC::callFunctionValue<IdList> ("meshobjectselection", "growVertexSelection", meshId);
        RPC::callFunctionValue<IdList> ("meshobjectselection", "growVertexSelection", meshId);
        IdList twoRingVertices = RPC::callFunctionValue<IdList> ("meshobjectselection", "getVertexSelection", meshId);

        
        if (twoRingVertices.size() < (d_fitting+1)*(d_fitting+2)/2) {
            QString str = QString("Not enough points for fitting for object %1: %2").arg(QString::number(meshId), meshName);
            emit log(LOGERR, str);
            continue;
        }

        std::vector<DPoint> in_points;
        for (size_t i = 0; i < twoRingVertices.size(); i++)
        {
            TriMesh::VertexHandle vh = mesh->vertex_handle(twoRingVertices[i]);
            for (TriMesh::VertexVertexIter vv_it = mesh->vv_iter(vh) ; vv_it; ++vv_it)
            {
                TriMesh::Point p = mesh->point(*vv_it);
                in_points.push_back(DPoint(p[0],p[1],p[2]));
            }
        }
        
        My_Monge_form monge_form;
        My_Monge_via_jet_fitting monge_fit;
        monge_form = monge_fit(in_points.begin(), in_points.end(), d_fitting, d_monge);
        
        int objId;
        emit addEmptyObject(DATA_COORDSYS, objId);
        CoordsysObject *object = 0;
        PluginFunctions::getObject(objId, object);
        CoordsysNode *coordSys = object->coordsysNode();
        ACG::Vec3d min_curvature(monge_form.minimal_principal_direction()[0], monge_form.minimal_principal_direction()[1], monge_form.minimal_principal_direction()[2]);
        ACG::Vec3d max_curvature(monge_form.maximal_principal_direction()[0], monge_form.maximal_principal_direction()[1], monge_form.maximal_principal_direction()[2]);
        ACG::Vec3d normal_curvature(monge_form.normal_direction()[0], monge_form.normal_direction()[1], monge_form.normal_direction()[2]);
        
        TriMesh::Point center = mesh->point(mesh->vertex_handle(currentVertex[0]));
        coordSys->position(center);
        
        ACG::Vec3d normal_mesh = mesh->normal(mesh->vertex_handle(currentVertex[0]));
        if (OpenMesh::dot(normal_mesh, normal_curvature) < 0) {
            ACG::Vec3d tmp = min_curvature;
            min_curvature = -max_curvature;
            max_curvature = -tmp;
            normal_curvature = -normal_curvature;
        }
        
        QMatrix4x4 QrotationM;
        QrotationM.lookAt(QVector3D(0,0,0),
                          QVector3D(max_curvature[0], max_curvature[1], max_curvature[2]),
                          QVector3D(min_curvature[0], min_curvature[1], min_curvature[2]));
        Matrix4x4 rotationM((TriMesh::Scalar*)QrotationM.data());
        rotationM.transpose();
        coordSys->rotation(rotationM);
        coordSys->size(indicatorSize);
    }
    
    // restore selection
    RPC::callFunctionValue<IdList> ("meshobjectselection", "clearVertexSelection", meshId);
    RPC::callFunctionValue<IdList> ("meshobjectselection", "selectVertices", meshId, originalVerticesSelection);
}

bool PlushPlugin::translate_openMesh_to_geodesic_mesh(TriMesh *mesh, std::vector<double> &points, std::vector<unsigned> &faces) {
    for (TriMesh::VertexIter v_it = mesh->vertices_begin(); v_it != mesh->vertices_end(); v_it++)
    {
        TriMesh::Point p = mesh->point(*v_it);
        points.push_back(p[0]);
        points.push_back(p[1]);
        points.push_back(p[2]);
    }
    for (TriMesh::FaceIter f_it = mesh->faces_begin(); f_it != mesh->faces_end(); f_it++)
    {
        for (TriMesh::FaceVertexIter fv_it = mesh->fv_iter(*f_it); fv_it; fv_it++)
        {
            faces.push_back(fv_it->idx());
        }
    }
    if (points.size() < 3 || faces.size() < 1)
    {
        return false;
    }
    return true;
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
            TriMesh *mesh = PluginFunctions::triMesh(*o_it);

            // pick vertices
            IdList selectedVertices = RPC::callFunctionValue<IdList> ("meshobjectselection", "getVertexSelection", o_it->id());

            if (selectedVertices.size() < 2) {
                QString str = QString("Object %1 contains less than 2 vertices. no path showed.").arg(o_it->name());
                emit log(LOGINFO, str);
                continue;
            }

            // select vertex with min Z coordinate
            double minZ = 1e9;
            size_t minZ_idx = 0;
            for (size_t sortI = 0; sortI < selectedVertices.size(); sortI++)
            {
                TriMesh::Point p1 = mesh->point(mesh->vertex_handle(selectedVertices[sortI]));
                if (p1[2] < minZ)
                {
                    minZ = p1[2];
                    minZ_idx = sortI;
                }
            }
            // select shortest path greedly
            for (size_t sortI = 0; sortI < selectedVertices.size()-1; sortI++)
            {
                TriMesh::Point p1 = mesh->point(mesh->vertex_handle(selectedVertices[sortI]));
                double minDist = 1e9;
                size_t min_idx = sortI+1;
                for (size_t sortJ = sortI+1; sortJ < selectedVertices.size(); sortJ++)
                {
                    TriMesh::Point p2 = mesh->point(mesh->vertex_handle(selectedVertices[sortJ]));
                    if ((p1-p2).length() < minDist)
                    {
                        minDist = (p1-p2).length();
                        min_idx = sortJ;
                    }
                }
                
                // exchange position
                size_t tmpId = selectedVertices[sortI+1];
                selectedVertices[sortI+1] = selectedVertices[min_idx];
                selectedVertices[min_idx] = tmpId;
            }
            
            for (size_t i = 1; i < selectedVertices.size(); i++)
            {
                int sourceIdx = selectedVertices[i-1];
                int destIdx = selectedVertices[i];

                // Init model if it's first time
                if (meshes_geodesic.find(meshId) == meshes_geodesic.end()) {
                    calcGeodesic(o_it);
                }
                geodesic::Mesh *mesh_geodesic = meshes_geodesic.at(meshId);
                //    geodesic::GeodesicAlgorithmExact algorithm(&gMesh);	//create exact algorithm for the mesh
                geodesic::GeodesicAlgorithmDijkstra algorithm(mesh_geodesic);	//create exact algorithm for the mesh
                
                geodesic::SurfacePoint source(&(mesh_geodesic->vertices()[sourceIdx]));		//create source
                std::vector<geodesic::SurfacePoint> all_sources(1,source);					//in general, there could be multiple sources, but now we have only one
                
                if(destIdx != -1)	//target vertex specified, compute single path
                {
                    geodesic::SurfacePoint target(&(mesh_geodesic->vertices()[destIdx]));		//create source
                    
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
            }
        }
    }
}
void PlushPlugin::calcGeodesic(PluginFunctions::ObjectIterator o_it) {
    if (o_it->dataType(DATA_TRIANGLE_MESH)) {
        int meshId = o_it->id();
        TriMesh *mesh = PluginFunctions::triMesh(*o_it);

        std::vector<double> points;
        std::vector<unsigned> faces;
        geodesic::Mesh *mesh_geodesic = new geodesic::Mesh();
        //geodesic::read_mesh_from_file(filename, points, faces);
        if (translate_openMesh_to_geodesic_mesh(mesh, points, faces))
        {
            mesh_geodesic->initialize_mesh_data(points, faces);		//create internal mesh data structure including edges
            meshes_geodesic.insert(std::pair<int, geodesic::Mesh*>(meshId, mesh_geodesic));
        }
    }
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