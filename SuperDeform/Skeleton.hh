#ifndef SKELETON_H
#define SKELETON_H

#include <vector>
#include <OpenFlipper/common/Types.hh>

using namespace std;

class Bone
{
public:
    Bone(OpenMesh::Vec3d _a, OpenMesh::Vec3d _b, int _idxa, int _idxb);
	~Bone();
	OpenMesh::Vec3d getA() const;
	void setA(OpenMesh::Vec3d _a);
	OpenMesh::Vec3d getB() const;
	void setB(OpenMesh::Vec3d _b);
	int getIdxA();
	int getIdxB();

	//add chi
	void updateA(double x , double y , double z);
	void updateB(double x , double y , double z);

private:
	OpenMesh::Vec3d a;
	OpenMesh::Vec3d b;
	int idxa;
	int idxb;
};

class Skeleton
{
public:
	Skeleton() {}
	~Skeleton();
	bool build(const std::string &filename);
	bool save(const std::string &filename);

	vector<Bone> bones;
	vector<OpenMesh::Vec3d> verts;
    vector<int> parentids;
};

#endif //SKELETON_H