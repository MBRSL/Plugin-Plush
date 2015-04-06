//
//  Geodesic.cc
//  OpenFlipper
//
//  Created by 饒亢 on 2015/1/13.
//
//

#include "PlushPatternGenerator.hh"
#include "WeightFunctor.hh"

#include <boost/graph/dijkstra_shortest_paths.hpp>
//#include <boost/property_map/function_property_map.hpp>

#include <QTextStream>

typedef boost::property_map<Polyhedron, boost::vertex_external_index_t>::type VertexIdPropertyMap;
typedef boost::property_map<Polyhedron, boost::edge_external_index_t>::type EdgeIdPropertyMap;

void PlushPatternGenerator::calcGeodesic(std::vector<VertexHandle> targetVertices)
{
    // Prepare property maps for shortest path
    VertexIdPropertyMap vertex_index_pmap = get(boost::vertex_external_index, m_polyhedron);
    EdgeIdPropertyMap edge_index_pmap = get(boost::edge_external_index, m_polyhedron);

    isJobCanceled = false;
    
    std::map<std::pair<VertexHandle, VertexHandle>, double> &geodesicDistance = m_mesh->property(geodesicDistanceHandle);
    std::map<std::pair<VertexHandle, VertexHandle>, std::vector<VertexHandle> > &geodesicPath = m_mesh->property(geodesicPathHandle);
    geodesicDistance.clear();
    geodesicPath.clear();

    for (size_t i = 0; i < targetVertices.size(); i++) {
        if (isJobCanceled) {
            emit log(LOGINFO, "Geodesic calculation canceled.");
            return;
        }

        VertexHandle sourceHandle = targetVertices[i];
        
        std::vector<boost_vertex_descriptor> predecessor(boost::num_vertices(m_polyhedron));
        boost::iterator_property_map<std::vector<boost_vertex_descriptor>::iterator, VertexIdPropertyMap> predecessor_pmap(predecessor.begin(), vertex_index_pmap);

        std::vector<double> distance(boost::num_vertices(m_polyhedron));
        boost::iterator_property_map<std::vector<double>::iterator, VertexIdPropertyMap> distance_pmap(distance.begin(), vertex_index_pmap);
        
        std::vector<double> weight(boost::num_edges(m_polyhedron));
        boost::iterator_property_map<std::vector<double>::iterator, EdgeIdPropertyMap> weight_pmap(weight.begin(), edge_index_pmap);
        
        boost_vertex_descriptor source = m_verticesMapping[targetVertices[i].idx()];
        
        // HACK: This is a local variable "shared" by both (functional)weightmap & visitor
        //       so that we can first assign it in visitor and then access it in weightmap
        boost_vertex_descriptor currentV = source;
        
        WeightFunctor weightFunctor(m_mesh,
                                    m_polyhedron,
                                    &currentV,
                                    predecessor_pmap,
                                    edge_index_pmap);
//        auto weight_pmap = boost::make_function_property_map< boost_edge_descriptor,
//                                                            double,
//                                                            WeightFunctor > (weightFunctor);

        // If smooth term is not used, pre-calculate edge weights can accelerate calculation
        boost_edge_iterator eit, eite;
        for (boost::tie(eit, eite) = boost::edges(m_polyhedron); eit != eite; eit++) {
            weight_pmap[*eit] = weightFunctor(*eit);
        }
//        Dijkstra_visitor visitor(&currentV);
        
        // Calculate geodesic using Dijkstra
        boost::dijkstra_shortest_paths(m_polyhedron, source,
                                       vertex_index_map(vertex_index_pmap)
                                       .weight_map(weight_pmap)
                                       .distance_map(distance_pmap)
                                       .predecessor_map(predecessor_pmap)
//                                       .visitor(visitor)
                                       );
        
        for (size_t j = 0; j < targetVertices.size(); j++) {
        VertexHandle destHandle = targetVertices[j];
    //    for (VertexIter v_it = m_mesh->vertices_begin(); v_it != m_mesh->vertices_end(); v_it++) {
    //        VertexHandle destHandle = *v_it;
    //        int id = v_it->idx();
            if (destHandle == sourceHandle) {
                continue;
            }

            // Directions. from source to destination, or from destination to source
            std::pair<VertexHandle, VertexHandle> edgeSD = std::make_pair(sourceHandle, destHandle);
            std::pair<VertexHandle, VertexHandle> edgeDS = std::make_pair(destHandle, sourceHandle);
            
            // Back tracking
            std::vector<VertexHandle> path;
            for(boost_vertex_descriptor predecessor = m_verticesMapping[destHandle.idx()]; predecessor_pmap[predecessor] != predecessor;) {
                path.push_back(m_mesh->vertex_handle(predecessor->id()));
                predecessor = predecessor_pmap[predecessor];
            }
            
            if (path.size() <= 0) {
                // Unreachable
                return;
            } else {
                // Add source to path
                path.push_back(sourceHandle);
            }
            // It's reversed using back tracking, flip it back.
            std::reverse(path.begin(), path.end());

            // DEBUGGING
            // The cost is sometimes non symmetric, this should not happen.
    //        if (geodesicDistance.find(edgeSD) != geodesicDistance.end()) {
    //            double prevDist = geodesicDistance[edgeSD];
    //            if (abs(prevDist - distance_pmap[m_verticesMapping[id]]) > 1e-5) {
    //                double newDist = distance_pmap[m_verticesMapping[id]];
    //            }
    //        }
            
            // TODO:
            // The true distance from A to B should be identical with distance from B to A
            // Currently we only choose the direction with smaller cost
            double cost = distance_pmap[m_verticesMapping[destHandle.idx()]];
            
            // First, check if opposite direction is already calculated
            std::map<std::pair<VertexHandle, VertexHandle>, double>::iterator cost_reverse_it = geodesicDistance.find(edgeDS);
            if (cost_reverse_it != geodesicDistance.end()) {
                double cost_reverse = cost_reverse_it->second;

                if (cost_reverse < cost) {
                    // Copy cost and path from reverse direction
                    geodesicDistance.insert(std::make_pair(edgeSD, cost_reverse));
                    
                    std::map<std::pair<VertexHandle, VertexHandle>, std::vector<VertexHandle> >::iterator path_reverse_it = geodesicPath.find(edgeDS);
                    // Path should exist because they are calculated together
                    assert(path_reverse_it != geodesicPath.end());
                    geodesicPath.insert(std::make_pair(edgeSD, path_reverse_it->second));
                } else {
                    // Replace reverse direction with current cost and path
                    geodesicDistance.erase(edgeDS);
                    geodesicPath.erase(edgeDS);
                    
                    geodesicDistance.insert(std::make_pair(edgeSD, cost));
                    geodesicDistance.insert(std::make_pair(edgeDS, cost));
                    geodesicPath.insert(std::make_pair(edgeSD, path));
                    std::reverse(path.begin(), path.end());
                    geodesicPath.insert(std::make_pair(edgeDS, path));
                }
            }
            // If not, just save it into property
            else {
                geodesicDistance.insert(std::make_pair(edgeSD, cost));
                geodesicPath.insert(std::make_pair(edgeSD, path));
            }
        }
        // most of the time is spent on geodesic calculation
        int status = (double)(i+1)/targetVertices.size() * 100;
        emit setJobState(status);
    }
}

bool PlushPatternGenerator::saveGeodesic(std::vector<VertexHandle> selectedVertices) {
    QString geodesicFilename = m_meshName + ".geodesic";
    QFile file(geodesicFilename);
    if (!file.open(QIODevice::WriteOnly | QIODevice::Text)) {
        emit log(LOGERR, QString("Error opening file: %1").arg(geodesicFilename));
        return false;
    }
    
    std::map<std::pair<VertexHandle, VertexHandle>, double> &geodesicDistance = m_mesh->property(geodesicDistanceHandle);
    std::map<std::pair<VertexHandle, VertexHandle>, std::vector<VertexHandle> > &geodesicPath = m_mesh->property(geodesicPathHandle);

    QTextStream out(&file);
    for (size_t i = 0; i < selectedVertices.size(); i++) {
        for (size_t j = i+1; j < selectedVertices.size(); j++) {
            VertexHandle sourceHandle = selectedVertices[i];
            VertexHandle destHandle = selectedVertices[j];

            // Directions. from source to destination, or from destination to source
            std::pair<VertexHandle, VertexHandle> edgeSD = std::make_pair(sourceHandle, destHandle);
            std::pair<VertexHandle, VertexHandle> edgeDS = std::make_pair(destHandle, sourceHandle);
            
            // source / dest / cost
            out << sourceHandle.idx() << " " << destHandle.idx() << " " << geodesicDistance[edgeSD] << endl;
            
            // num of Vertices on path/ v1 v2 ...
            std::vector<VertexHandle> path = geodesicPath[edgeSD];
            out << path.size();
            for (size_t pi = 0; pi < path.size(); pi++) {
                out << " " << path[pi].idx();
            }
            out << endl;
        }
    }
    file.close();
    return true;
}

bool PlushPatternGenerator::loadGeodesic() {
    QString geodesicFilename = m_meshName + ".geodesic";
    
    QFile file(geodesicFilename);
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
        emit log(LOGERR, QString("Error opening file: %1").arg(geodesicFilename));
        return false;
    }
    
    std::map<std::pair<VertexHandle, VertexHandle>, double> &geodesicDistance = m_mesh->property(geodesicDistanceHandle);
    std::map<std::pair<VertexHandle, VertexHandle>, std::vector<VertexHandle> > &geodesicPath = m_mesh->property(geodesicPathHandle);
    
    QTextStream fin(&file);
    while (!fin.atEnd()) {
        int sourceIdx, destIdx;
        double cost;
        fin >> sourceIdx >> destIdx >> cost;
        
        VertexHandle sourceHandle = m_mesh->vertex_handle(sourceIdx);
        VertexHandle destHandle = m_mesh->vertex_handle(destIdx);

        // Directions. from source to destination, or from destination to source
        std::pair<VertexHandle, VertexHandle> edgeSD = std::make_pair(sourceHandle, destHandle);
        std::pair<VertexHandle, VertexHandle> edgeDS = std::make_pair(destHandle, sourceHandle);
        
        geodesicDistance.erase(edgeSD);
        geodesicDistance.erase(edgeDS);
        geodesicDistance.emplace(edgeSD, cost);
        geodesicDistance.emplace(edgeDS, cost);

        // consume endl
        fin.readLine();
        QString pathStr = fin.readLine();
        QStringList vertices = pathStr.split(" ");
        // The first one is number of vertices on path
        int num = vertices[0].toInt();
        
        std::vector<VertexHandle> path;
        for (int i = 1; i <= num; i++) {
            path.push_back(m_mesh->vertex_handle(vertices[i].toInt()));
        }
        
        geodesicPath.erase(edgeSD);
        geodesicPath.erase(edgeDS);
        geodesicPath.emplace(edgeSD, path);
        reverse(path.begin(), path.end());
        geodesicPath.emplace(edgeDS, path);
    }
    file.close();
    return true;
}
