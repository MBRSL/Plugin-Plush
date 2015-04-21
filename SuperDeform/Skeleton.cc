#include "Skeleton.hh"
#include <iostream>
#include <fstream>

Bone::Bone(OpenMesh::Vec3d _a, OpenMesh::Vec3d _b, int _idxa, int _idxb) :
	  a(_a), b(_b), idxa(_idxa), idxb(_idxb)
{
}

Bone::~Bone()
{
}

//add chi
void Bone::updateA(double x ,double y , double z) {
	a = OpenMesh::Vec3d(x,y,z);
}
void Bone::updateB(double x ,double y , double z) {
	b = OpenMesh::Vec3d(x,y,z);
}

OpenMesh::Vec3d Bone::getA() const
{
	return a;
}

void Bone::setA(OpenMesh::Vec3d _a)
{
	a = _a;
}

OpenMesh::Vec3d Bone::getB() const
{
	return b;
}

void Bone::setB(OpenMesh::Vec3d _b)
{
	b = _b;
}

int Bone::getIdxA()
{
	return idxa;
}

int Bone::getIdxB()
{
	return idxb;
}

Skeleton::~Skeleton()
{
}

bool Skeleton::save(const std::string &filename)
{
	ofstream outputFile(filename.c_str());
	if (!outputFile.is_open())
	{
		cout << "Error opening file " << filename.c_str() << endl;
		return false;
	}
	if (verts.empty())
	{
		return false;
	}

	outputFile << "0 " << verts[0][0] << " " << verts[0][1] <<
		" " << verts[0][2] << " -1" << endl;
	for(size_t i = 0; i < bones.size(); i++)
	{
		outputFile << bones[i].getIdxB() << " " << bones[i].getB()[0] << " " << bones[i].getB()[1] <<
			" " << bones[i].getB()[2] << " " << bones[i].getIdxA() << endl;
	}
	outputFile.close();
	return true;
}

bool Skeleton::build(const std::string &filename)
{
    ifstream fin(filename.c_str());
    
    if(!fin.is_open())
    {
        cout << "Error opening file " << filename.c_str() << endl;
        return false;
    }
    
    std::vector<int> ids_tmp;
    std::vector<OpenMesh::Vec3d> verts_tmp;
    std::vector<int> parentids_tmp;
    while(!fin.eof())
    {
        int id = -1;
        int parentid = -1;
        double d[3];
        fin >> id >> d[0] >> d[1] >> d[2] >> parentid;
        if (id == -1)
        {
            break;
        }
        ids_tmp.push_back(id);
        OpenMesh::Vec3d p(d[0], d[1], d[2]);
        verts_tmp.push_back(p);
        parentids_tmp.push_back(parentid);
    }
    
    verts.resize(ids_tmp.size());
    parentids.resize(ids_tmp.size());
    for (size_t i = 0; i < ids_tmp.size(); i++)
    {
        size_t id = ids_tmp[i];
        assert(0 <= id && id < ids_tmp.size() && "Bone index should be [0, #bones)");
        verts[id] = verts_tmp[i];
        parentids[id] = parentids_tmp[id];
        assert(parentids[id] < (int)ids_tmp.size() &&  "Parent bone index should be (-inf, #bones)");
    }
    for (size_t i = 0; i < ids_tmp.size(); i++)
    {
        int id = ids_tmp[i];
        if(parentids[id] >= 0)
        {
            Bone bone(verts[parentids[id]], verts[id], parentids[id], id);
            assert((bone.getA() - bone.getB()).norm() > 0 && "A bone contains no length. Do you have identical joints?");
            bones.push_back(bone);
        }
    }
    
    return true;
}