#include "Skeleton.hh"
#include "PlushPlugin.hh"

#include <Eigen/Sparse>
//#include "includelib.h"
//#include "vector3d.h"
//#include "taucs_interface.h"

class Weight
{
public:
	Weight() 
	{
//		InitTaucsInterface();
	}
	~Weight();
//	void computeBoneWeight(Polyhedron *P, const vector<Vector> &facetNormal, Skeleton *skeleton);
    void computeBoneWeight(TriMesh *mesh, Skeleton *skeleton);

	Eigen::MatrixXd getWeight();

private:
//	void setMatrix(Polyhedron *P, const std::vector<Vector> &facetNormal, const vector<Bone *> &bones,
    void setMatrix(TriMesh *mesh, const vector<Bone> &bones,
                   Eigen::SparseMatrix<double> &idDA, Eigen::VectorXd &H, std::vector< std::vector<int> > &nearestBoneIdx);
	// collisionDetection
//	bool collision(Polyhedron *P, const std::vector<Vector> &facetNormal, Point3D p1, Point3D p2);
    bool collision(TriMesh *mesh, const std::vector<OpenMesh::Vec3d> &facetNormal, OpenMesh::Vec3d p1, OpenMesh::Vec3d p2);
    void solveWeight(int vsize, int bsize, Eigen::SparseMatrix<double> &idDA, Eigen::VectorXd &H, std::vector<std::vector<int> > &nearestBoneIdx);
	void normalizeWeight(int vsize);

    Eigen::MatrixXd boneweight;
};