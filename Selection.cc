//
//  Selection.cpp
//  OpenFlipper
//
//  Created by 饒亢 on 2015/1/31.
//
//

#include "PlushPlugin.hh"

int PlushPlugin::loadSelection(int meshId, QString meshName) {
    QFile file(meshName+"_selection.txt");
    if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
        emit log(LOGERR, QString("Unable to read file %1").arg(meshName+"_selection.txt"));
        return;
    }
    
    QTextStream in(&file);
    
    IdList selectedVertices;
    QString v;
    while(!in.atEnd()) {
        v = in.readLine();
        bool ok;
        int vId = v.toInt(&ok);
        if (ok) {
            selectedVertices.push_back(vId);
        }
    }
    RPC::callFunction<int, IdList> ("meshobjectselection", "selectVertices", meshId, selectedVertices);
    return selectedVertices.size();
}

void PlushPlugin::saveSelection(int meshId, QString meshName) {
    IdList selectedVertices = RPC::callFunctionValue<IdList> ("meshobjectselection", "getVertexSelection", meshId);
    
    // Prepare file for saving data
    QFile file(meshName+"_selection.txt");
    file.open(QIODevice::WriteOnly | QIODevice::Text);
    QTextStream out(&file);
    
    for (IdList::iterator v_it = selectedVertices.begin(); v_it != selectedVertices.end(); v_it++) {
        out << *v_it << "\n";
    }
    file.close();
}

void PlushPlugin::clearSelection(int meshId) {
    RPC::callFunction<int> ("meshobjectselection", "clearVertexSelection", meshId);
}