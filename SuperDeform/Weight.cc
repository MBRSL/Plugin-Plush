#include "Weight.hh"
//
//void Weight::computeBoneWeight(TriMesh *P, const std::vector<OpenMesh::Vec3d> &facetNormal, Skeleton *skeleton)
//{
//	if (!boneweight.empty())
//	{
//		for(int i = 0; i < boneweight.size(); i++)
//		{
//			delete [] boneweight[i];
//			boneweight[i] = NULL;
//		}
//		boneweight.clear();
//	}
//	for(size_t i = 0; i < skeleton->bones.size(); i++)
//	{
//		taucsType *w = new taucsType[P->size_of_vertices()];
//		boneweight.push_back(w);
//	}
//	
//	int idDA = CreateMatrix(P->size_of_vertices(), P->size_of_vertices());
//	taucsType *H = new taucsType[P->size_of_vertices()]; // diagonal matrix
//    std::vector<std::vector<int> > nearestBoneIdx;
//	setMatrix(P, facetNormal, skeleton->bones, idDA, H, nearestBoneIdx);
//	solveWeight(P->size_of_vertices(), skeleton->bones.size(), idDA, H, nearestBoneIdx);
//	normalizeWeight(P->size_of_vertices());
//
//	ReleaseMatrix(idDA);
//	delete [] H;
//	for(size_t i = 0; i < nearestBoneIdx.size(); i++)
//	{
//		nearestBoneIdx[i].clear();
//	}
//	nearestBoneIdx.clear();
//}
//
//Weight::~Weight()
//{
//	for(int i = 0; i < boneweight.size(); i++)
//	{
//		delete [] boneweight[i];
//		boneweight[i] = NULL;
//	}
//	boneweight.clear();
//	DeinitTaucsInterface();
//}
//
//taucsType* Weight::getWeight(int idx)
//{
//	return boneweight[idx];
//}
//
//void Weight::setMatrix(Polyhedron *P, const vector<Vector> &facetNormal,const vector<Bone *> &bones, 
//						   int idDA, taucsType *H, vector<vector<int>> &nearestBoneIdx)
//{
//	// cotagent Laplacian
//	// theta_i = Laplacian of v_i
//	// theta_i = [sum of edges connected by i and j (sum w_ij * v_j)] - v_i
//	// w_ij = omega_ij / [sum of edges connected by i and k (sum omega_ik)]
//	// omega_ij = cot(alpha) + cot(beta)
//	// therefore, L_ij = -1   if i == j
//	//                 = w_ij if (i,j) is edge
//	//                 = 0    otherwise
//
//	int i = 0;
//	const int idA = CreateMatrix(P->size_of_vertices(), P->size_of_vertices());
//	taucsType *D = new taucsType[P->size_of_vertices()];
//	
//	for (i = 0; i < P->size_of_vertices(); i++)
//	{
//		H[i] = 0.0;
//		D[i] = 0.0;
//	}
//
//	// start set D = diagonal matrix, d(i) = 1/area
//	double *triArea = new double[facetNormal.size()];
//	for (i = 0; i < facetNormal.size(); i++)
//	{
//		triArea[i] = sqrt(facetNormal[i] * facetNormal[i]);
//	}
//
//	i = 0;
//	for (Facet_iterator fi = P->facets_begin(); fi != P->facets_end(); fi++, i++)
//	{
//		Halfedge_handle he = fi->halfedge();
//
//		for (int j = 0; j < 3; j++)
//		{
//			D[he->vertex()->id()] += triArea[i]; // d(i) = area
//			he = he->next();
//		}
//	}
//	delete [] triArea;
//	// end set D = diagonal matrix, d(i) = 1/area
//
//	i = 0;
//	for(Vertex_iterator vi = P->vertices_begin(); vi != P->vertices_end(); vi++, i++)
//	{
//		Halfedge_handle he = vi->halfedge();
//		Point3D pi = Point3D(vi->point()[0], vi->point()[1], vi->point()[2]);
//
//		// set diagonal matrix H
//		// H = kc/d(j)^2
//		// d(j) = the distance from vertex j to the nearest bone
//		double* piToBone_square = new double[bones.size()];
//		double minDist_square = 1e37;
//		for(int j = 0; j < bones.size(); j++)
//		{
//			piToBone_square[j] = distToSeg_square(pi, *bones[j]->getA(), *bones[j]->getB());
//			minDist_square = min(piToBone_square[j], minDist_square);			
//		}
//		double threshold = minDist_square * 1.000001;
//		if(minDist_square < 1e-8)
//		{
//			minDist_square = 1e-8;
//		}
//		vector <int> tempNearestBoneIdx;
//		for(int j = 0; j < bones.size(); j++)
//		{
//			//the reason we don't just pick the closest bone is so that if two are
//			//equally close, both are factored in.
//			if (piToBone_square[j] <= threshold)
//			{
//				//if (!collision(P, facetNormal, pi, 
//				//	projToSeg(pi, *bones[j]->getA(), *bones[j]->getB())))
//				{
//					H[i] += 1.0 / minDist_square;
//					tempNearestBoneIdx.push_back(j);
//				}
//			}	
//		}
//		nearestBoneIdx.push_back(tempNearestBoneIdx);;
//		delete[] piToBone_square;
//		tempNearestBoneIdx.clear();
//		// end set diagonal matrix H
//
//		// -L+H = DA
//		// L = cotangent matrix
//		// H = diagonal matrix
//		// D = diagonal matrix, d(i) = 1/area
//		// A = (sum omega_ik+H(i)/D(i)) if i == j
//		//   = -omega_ij                 if (i,j) is edge
//		//   = 0                         otherwise
//		// start set matrix A and diagonal matrix D
//		double sumOmega = 0.0;
//		do{
//			if(he->next() != NULL)
//				he = he->next();
//			else
//				break;
//
//			// start calculate omega_ij = cot(alpha) + cot(beta)
//			Vertex_handle vj = he->vertex(); // incident vertex, vi's neighbor
//			Point3D pj = Point3D(vj->point()[0], vj->point()[1], vj->point()[2]);
//			double cot_alpha = 0;
//			double cot_beta = 0;
//			if(he->next() != NULL)
//			{
//				Vertex_handle vij_1 = he->next()->vertex();
//				Point3D pij_1 = Point3D(vij_1->point()[0], vij_1->point()[1], vij_1->point()[2]);
//				cot_alpha = cot(pi - pij_1, pj - pij_1);
//			}
//			if(he->opposite() != NULL)
//			{
//				if(he->opposite()->next() != NULL)
//				{
//					Vertex_handle vij_2 = he->opposite()->next()->vertex();
//					Point3D pij_2 = Point3D(vij_2->point()[0], vij_2->point()[1], vij_2->point()[2]);
//					cot_beta = cot(pi - pij_2, pj - pij_2);
//				}
//			}
//			double omega_ij = cot_alpha + cot_beta;
//			// end calculate omega_ij = cot(alpha) + cot(beta)
//
//			SetMatrixEntry(idA, i, vj->id(), -omega_ij);
//			sumOmega += omega_ij;
//
//			if(he->opposite() != NULL)
//				he = he->opposite();
//			else 
//				break;
//		} while(vi->halfedge() != he);
//
//		if(D[i] < 1e-7)
//			D[i] = 1e-7;
//		D[i] = 1.0 / D[i]; // d(i) = 1/area
//
//		AddToMatrixEntry(idA, i, i, sumOmega + H[i] / D[i]);
//		// end set matrix A and diagonal matrix D
//	}
//
//	if(!MultiplyDiagMatrixMatrix(idA, D, idDA))
//	{
//		std::cout << "failed to multiply diagonal matrix matrix.\n";
//		return;
//	}
//
//	ReleaseMatrix(idA);
//	delete [] D;
//}
//
//// collision detection
//bool Weight::collision(Polyhedron *P, const vector<Vector> &facetNormal, Point3D p1, Point3D p2)
//{
//	int i = 0;
//	for (Facet_iterator fi = P->facets_begin(); fi != P->facets_end(); fi++, i++)
//	{
//		Halfedge_handle he = fi->halfedge();
//		Point3D tri[3];
//		for (int j = 0; j < 3; j++)
//		{
//			tri[j] = Point3D(he->vertex()->point()[0], he->vertex()->point()[1], he->vertex()->point()[2]);
//			he = he->next();
//		}
//
//		// line segment p1 to p2 intersect with plane of triangle 
//		// compute distances of p1 to plane and p2 to plane
//		double d = 0.0;
//		for (int j = 0; j < 3; ++j)
//		{
//			d += facetNormal[i][j] * tri[0][j];
//		}
//		double d1 = -d;
//		double d2 = -d;
//		for (int j = 0; j < 3; ++j)
//		{
//			d1 += facetNormal[i][j] * p1[j];
//			d2 += facetNormal[i][j] * p2[j];
//		}
//		if(d1 * d2 >= 0) // p1/p2 on the plane or on the same side of the plane
//		{
//			continue;
//		}
//		// get the intersection point
//		Point3D pintersect = ((p1 * d2) + (p2 * d1)) / (d1 + d2);
//
//		// check whether pintersect is inside the triangle
//		Vector3D vu = tri[2] - tri[0];
//		Vector3D vv = tri[1] - tri[0];
//		Vector3D vo = pintersect - tri[0];
//		double dd = (vu[0] * vv[1]) - (vu[1] * vv[0]);
//		double su = ((vo[0] * vv[1]) - (vo[1] * vv[0])) / dd;
//		if(su < 0 || su > 1) // outside the triangle
//		{
//			continue;
//		}
//		double sv = ((vo[0] * vu[1]) - (vo[1] * vu[0])) / (-dd);
//		if(sv < 0 || sv > 1) // outside the triangle
//		{
//			continue;
//		}
//		double sum = su + sv;
//		if(sum <= 1 && sum >= 0) // inside the triangle
//		{
//			return true;
//		}
//	}
//	return false;
//}
//
//void Weight::solveWeight(int vsize, int bsize, int idDA, taucsType *H, vector<vector<int>> &nearestBoneIdx)
//{
//	// idDA * wi = Hpi
//	taucsType * Hp = new taucsType[vsize];
//	for(int i = 0; i < bsize; i++)
//	{
//		// start set vector Hpi
//		for (int j = 0; j < vsize; j++)
//		{
//			Hp[j] = 0.0;
//			for(int k = 0; k < nearestBoneIdx[j].size(); k++)
//			{
//				if(nearestBoneIdx[j][k] == i)
//				{
//					Hp[j] = H[j] / nearestBoneIdx[j].size();
//					break;
//				}
//			}	
//		}
//		// end set vector Hpi
//
//		if(!SolveATA(idDA, Hp, boneweight[i], 1))
//		{
//			std::cout << "Failed to solve ATA of weight " << i << std::endl;
//		}
//	}
//
//	delete [] Hp;
//}
//
//void Weight::normalizeWeight(int vsize)
//{
//	for(int i = 0; i < vsize; i++) // vertex
//	{
//		double sum = 0.0;
//		for (int j = 0; j < boneweight.size(); j++) // bone
//		{
//			if(boneweight[j][i] > 1)
//			{
//				boneweight[j][i] = 1;
//			}
//			else if(boneweight[j][i] < 0)
//			{
//				boneweight[j][i] = 0;
//			}
//			sum += boneweight[j][i];
//		}
//		for (int j = 0; j < boneweight.size(); j++) // bone
//		{
//			boneweight[j][i] = boneweight[j][i] / sum;
//		}
//	}
//}