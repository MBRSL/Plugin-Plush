#include "Skeleton.hh"
#include "PlushPlugin.hh"

#include <Eigen/Sparse>

class Weight
{
public:
	Weight() 
	{
	}
	~Weight();
    void computeBoneWeight(TriMesh *mesh, Skeleton *skeleton);

	Eigen::MatrixXd getWeight();

private:
    void setMatrix(TriMesh *mesh, const vector<Bone> &bones,
                   Eigen::SparseMatrix<double> &idDA, Eigen::VectorXd &H, std::vector< std::vector<int> > &nearestBoneIdx);
	// collisionDetection
    bool collision(TriMesh *mesh, const std::vector<OpenMesh::Vec3d> &facetNormal, OpenMesh::Vec3d p1, OpenMesh::Vec3d p2);
    void solveWeight(int vsize, int bsize, Eigen::SparseMatrix<double> &idDA, Eigen::VectorXd &H, std::vector<std::vector<int> > &nearestBoneIdx);
	void normalizeWeight(int vsize);

    Eigen::MatrixXd boneweight;
};