#include "Weight.hh"

double distToSeg_square(OpenMesh::Vec3d p, OpenMesh::Vec3d a, OpenMesh::Vec3d b)
{
    OpenMesh::Vec3d seg = b - a;
    OpenMesh::Vec3d diffbp = b - p;
    
    if((diffbp | seg) <= 0)
        return diffbp | diffbp;
    
    OpenMesh::Vec3d diffpa = p - a;
    double dot_diffpa_seg = diffpa | seg;
    if(dot_diffpa_seg <= 0)
        return diffpa | diffpa;
    
    return max(0.0, (diffpa | diffpa) - (dot_diffpa_seg * dot_diffpa_seg) / (seg | seg));
}

double cot(OpenMesh::Vec3d v1, OpenMesh::Vec3d v2)
{
    double dot_v = v1 | v2;
    return dot_v / (1e-6 + sqrt((v1 | v1) * (v2 | v2) - dot_v * dot_v));
}

void Weight::computeBoneWeight(TriMesh *mesh, Skeleton *skeleton)
{
//	if (!boneweight.empty())
//	{
//		for(size_t i = 0; i < boneweight.size(); i++)
//		{
//			delete [] boneweight[i];
//			boneweight[i] = NULL;
//		}
//		boneweight.clear();
//	}
    boneweight = Eigen::MatrixXd(skeleton->bones.size(), mesh->n_vertices());
//	for(size_t i = 0; i < skeleton->bones.size(); i++)
//	{
//        Eigen::VectorXd w = Eigen::VectorXd(mesh->n_vertices());
//		boneweight.push_back(w);
//	}
	
//    Eigen::MatrixXd idDA = Eigen::MatrixXd(mesh->n_vertices(), mesh->n_vertices());
    Eigen::SparseMatrix<double> idDA(mesh->n_vertices(), mesh->n_vertices());
    Eigen::VectorXd H = Eigen::VectorXd(mesh->n_vertices()); // diagonal matrix
    std::vector<std::vector<int> > nearestBoneIdx;
	setMatrix(mesh, skeleton->bones, idDA, H, nearestBoneIdx);
	solveWeight(mesh->n_vertices(), skeleton->bones.size(), idDA, H, nearestBoneIdx);
    normalizeWeight(mesh->n_vertices());

//	ReleaseMatrix(idDA);
//	delete [] H;
	for(size_t i = 0; i < nearestBoneIdx.size(); i++)
	{
		nearestBoneIdx[i].clear();
	}
	nearestBoneIdx.clear();
}

Weight::~Weight()
{
//	for(size_t i = 0; i < boneweight.size(); i++)
//	{
//		delete [] boneweight[i];
//		boneweight[i] = NULL;
//	}
//	boneweight.clear();
//	DeinitTaucsInterface();
}

Eigen::MatrixXd Weight::getWeight()
{
	return boneweight;
}

void Weight::setMatrix(TriMesh *mesh, const vector<Bone> &bones,
						   Eigen::SparseMatrix<double> &idDA, Eigen::VectorXd &H, vector< vector<int> > &nearestBoneIdx)
{
	// cotagent Laplacian
	// theta_i = Laplacian of v_i
	// theta_i = [sum of edges connected by i and j (sum w_ij * v_j)] - v_i
	// w_ij = omega_ij / [sum of edges connected by i and k (sum omega_ik)]
	// omega_ij = cot(alpha) + cot(beta)
	// therefore, L_ij = -1   if i == j
	//                 = w_ij if (i,j) is edge
	//                 = 0    otherwise

	size_t i = 0;
//    Eigen::SparseMatrix<double> idA(mesh->n_vertices(), mesh->n_vertices());
    Eigen::VectorXd D(mesh->n_vertices());
	
	for (i = 0; i < mesh->n_vertices(); i++)
	{
		H[i] = 0.0;
		D[i] = 0.0;
	}

	// start set D = diagonal matrix, d(i) = 1/area
//	for (i = 0; i < mesh->OpenMesh::BaseKernel::n_faces(); i++)
//	{
//		triArea[i] = sqrt(facetNormal[i] | facetNormal[i]);
//	}

	i = 0;
	for (FaceIter f_it = mesh->faces_begin(); f_it != mesh->faces_end(); f_it++, i++)
	{
        OpenMesh::Vec3d normal = mesh->property(mesh->face_normals_pph(), *f_it);
        double triArea = sqrt(normal | normal);

        HalfedgeHandle he = mesh->halfedge_handle(*f_it);
        
		for (int j = 0; j < 3; j++)
		{
			D(mesh->to_vertex_handle(he).idx()) += triArea; // d(i) = area
			he = mesh->next_halfedge_handle(he);
		}
	}
	// end set D = diagonal matrix, d(i) = 1/area

    
    // Prepare for filling matrix
    std::vector< Eigen::Triplet<double> > tripletList;
	i = 0;
	for(VertexIter v_it = mesh->vertices_begin(); v_it != mesh->vertices_end(); v_it++, i++)
	{
        HalfedgeHandle he = mesh->halfedge_handle(*v_it);
        TriMesh::Point pi = mesh->point(*v_it);
        
		// set diagonal matrix H
		// H = kc/d(j)^2
		// d(j) = the distance from vertex j to the nearest bone
		double* piToBone_square = new double[bones.size()];
		double minDist_square = 1e37;
		for(size_t j = 0; j < bones.size(); j++)
		{
			piToBone_square[j] = distToSeg_square(pi, bones[j].getA(), bones[j].getB());
			minDist_square = min(piToBone_square[j], minDist_square);			
		}
		double threshold = minDist_square * 1.000001;
		if(minDist_square < 1e-8)
		{
			minDist_square = 1e-8;
		}
		vector <int> tempNearestBoneIdx;
		for(size_t j = 0; j < bones.size(); j++)
		{
			//the reason we don't just pick the closest bone is so that if two are
			//equally close, both are factored in.
			if (piToBone_square[j] <= threshold)
			{
				//if (!collision(P, facetNormal, pi, 
				//	projToSeg(pi, *bones[j]->getA(), *bones[j]->getB())))
				{
					H[i] += 1.0 / minDist_square;
					tempNearestBoneIdx.push_back(j);
				}
			}	
		}
		nearestBoneIdx.push_back(tempNearestBoneIdx);
		delete[] piToBone_square;
		tempNearestBoneIdx.clear();
		// end set diagonal matrix H

		// -L+H = DA
		// L = cotangent matrix
		// H = diagonal matrix
		// D = diagonal matrix, d(i) = 1/area
		// A = (sum omega_ik+H(i)/D(i)) if i == j
		//   = -omega_ij                 if (i,j) is edge
		//   = 0                         otherwise
		// start set matrix A and diagonal matrix D
		double sumOmega = 0.0;
		do{
//            if(he->next() != NULL)
			if(true)
				he = mesh->next_halfedge_handle(he);
			else
				break;

			// start calculate omega_ij = cot(alpha) + cot(beta)
			VertexHandle vj = mesh->to_vertex_handle(he); // incident vertex, vi's neighbor
            TriMesh::Point pj = mesh->point(vj);
			double cot_alpha = 0;
			double cot_beta = 0;
//			if(he->next() != NULL)
            if (true)
			{
				VertexHandle vij_1 = mesh->to_vertex_handle(mesh->next_halfedge_handle(he));
                TriMesh::Point pij_1 = mesh->point(vij_1);
                cot_alpha = cot(pi - pij_1, pj - pij_1);
			}
//			if(he->opposite() != NULL)
            if (true)
			{
//				if(he->opposite()->next() != NULL)
                if (true)
				{
                    VertexHandle vij_2 = mesh->to_vertex_handle(mesh->next_halfedge_handle(mesh->opposite_halfedge_handle(he)));
                    TriMesh::Point pij_2 = mesh->point(vij_2);
                    cot_beta = cot(pi - pij_2, pj - pij_2);
				}
			}
			double omega_ij = cot_alpha + cot_beta;
			// end calculate omega_ij = cot(alpha) + cot(beta)

//			SetMatrixEntry(idA, i, vj->id(), -omega_ij);
            tripletList.push_back(Eigen::Triplet<double>(i, vj.idx(), -omega_ij));
			sumOmega += omega_ij;

//			if(he->opposite() != NULL)
            if (true) {
				he = mesh->opposite_halfedge_handle(he);
            } else {
				break;
            }
		} while(mesh->halfedge_handle(*v_it) != he);

		if(D(i) < 1e-7)
			D(i) = 1e-7;
		D(i) = 1.0 / D(i); // d(i) = 1/area

//		idA(i, i) += (sumOmega + H(i) / D(i));
        tripletList.push_back(Eigen::Triplet<double>(i, i, sumOmega + H(i) / D(i)));
		// end set matrix A and diagonal matrix D
	}

//	if(!MultiplyDiagMatrixMatrix(idA, D, idDA))
    idDA.setFromTriplets(tripletList.begin(), tripletList.end());
//    idDA = Eigen::MatrixXd(idA * D.asDiagonal());
    if (false)
	{
		std::cout << "failed to multiply diagonal matrix matrix.\n";
		return;
	}
}

// collision detection
bool Weight::collision(TriMesh *mesh, const vector<OpenMesh::Vec3d> &facetNormal, TriMesh::Point p1, TriMesh::Point p2)
{
	int i = 0;
	for (FaceIter f_it = mesh->faces_begin(); f_it != mesh->faces_end(); f_it++, i++)
	{
		HalfedgeHandle he = mesh->halfedge_handle(*f_it);
		TriMesh::Point tri[3];
		for (int j = 0; j < 3; j++)
		{
            tri[j] = mesh->point(mesh->to_vertex_handle(he));
			he = mesh->next_halfedge_handle(he);
		}

		// line segment p1 to p2 intersect with plane of triangle 
		// compute distances of p1 to plane and p2 to plane
		double d = 0.0;
		for (int j = 0; j < 3; ++j)
		{
			d += facetNormal[i][j] * tri[0][j];
		}
		double d1 = -d;
		double d2 = -d;
		for (int j = 0; j < 3; ++j)
		{
			d1 += facetNormal[i][j] * p1[j];
			d2 += facetNormal[i][j] * p2[j];
		}
		if(d1 * d2 >= 0) // p1/p2 on the plane or on the same side of the plane
		{
			continue;
		}
		// get the intersection point
		TriMesh::Point pintersect = ((p1 * d2) + (p2 * d1)) / (d1 + d2);

		// check whether pintersect is inside the triangle
		TriMesh::Point vu = tri[2] - tri[0];
		TriMesh::Point vv = tri[1] - tri[0];
		TriMesh::Point vo = pintersect - tri[0];
		double dd = (vu[0] * vv[1]) - (vu[1] * vv[0]);
		double su = ((vo[0] * vv[1]) - (vo[1] * vv[0])) / dd;
		if(su < 0 || su > 1) // outside the triangle
		{
			continue;
		}
		double sv = ((vo[0] * vu[1]) - (vo[1] * vu[0])) / (-dd);
		if(sv < 0 || sv > 1) // outside the triangle
		{
			continue;
		}
		double sum = su + sv;
		if(sum <= 1 && sum >= 0) // inside the triangle
		{
			return true;
		}
	}
	return false;
}

void Weight::solveWeight(int vsize, int bsize, Eigen::SparseMatrix<double> &idDA, Eigen::VectorXd &H, vector< vector<int> > &nearestBoneIdx)
{
    Eigen::SimplicialLDLT< Eigen::SparseMatrix<double> > solver;
    solver.compute(idDA);
    //        ldlt.compute(idDA.transpose() * idDA);
    if (solver.info() != Eigen::Success)
    {
        std::cout << "Failed to decompose matrix " << std::endl;
        return;
    }

	// idDA * wi = Hpi
	Eigen::VectorXd Hp = Eigen::VectorXd(vsize);
	for(int i = 0; i < bsize; i++)
	{
		// start set vector Hpi
		for (int j = 0; j < vsize; j++)
		{
			Hp[j] = 0.0;
			for(size_t k = 0; k < nearestBoneIdx[j].size(); k++)
			{
				if(nearestBoneIdx[j][k] == i)
				{
					Hp[j] = H[j] / nearestBoneIdx[j].size();
					break;
				}
			}	
		}
		// end set vector Hpi

//		if(!SolveATA(idDA, Hp, boneweight[i], 1))
//        boneweight[i] = Eigen::VectorXd((idDA->transpose() * (*idDA)).ldlt().solve(idDA->transpose() * Hp));
        boneweight.block(i,0,1,vsize) = Eigen::VectorXd(solver.solve(Hp)).transpose();
	}
}

void Weight::normalizeWeight(int vsize)
{
	for(int i = 0; i < vsize; i++) // vertex
	{
		double sum = 0.0;
		for (int j = 0; j < boneweight.rows(); j++) // bone
		{
			if(boneweight(j,i) > 1)
			{
				boneweight(j,i) = 1;
			}
			else if(boneweight(j,i) < 0)
			{
				boneweight(j,i) = 0;
			}
			sum += boneweight(j,i);
		}
		for (int j = 0; j < boneweight.rows(); j++) // bone
		{
			boneweight(j,i) = boneweight(j,i) / sum;
		}
	}
}