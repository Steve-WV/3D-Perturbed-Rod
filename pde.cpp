// Comment: 

// The following code uses the open source C++-package on cubic spline interpolation included in the class
// spline.h which is made available on: https://kluge.in-chemnitz.de/opensource/spline/

#include <iostream>
#include <fstream>
#include <vector>
#include <iostream>
#include <string>
#include <random>
#include <ctime>        // std::time
#include <cstdlib>     
#include <algorithm> 
#include <cmath> 
#include <cstdio>
#include <vector>
#include<math.h>
#include<stdio.h>
#include<cctype>
using namespace std;

#include <Eigen/Core>
#include <Eigen/SparseCore>

using namespace Eigen;

#include "pde.hpp"
#include "mesh.hpp"
#include "vertex.hpp"
#include "tmv.hpp"
#include "spline.h"
#include "boost/tokenizer.hpp"

#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>




	// Generators for the distributions 
		std::default_random_engine generator_norm1;
		std::default_random_engine generator_norm2;

	// generator for standard normal distribution
 	std::normal_distribution<double> distribution_norm1(0.0,1.0);
					


inline double pos_part(double s)
{
		if (s>=0) {return s;} 

		else

		return 0.0; 
	}


inline double GKernel_Smoothing(std::vector<double> th1, std::vector<double> th2, double s, double b) {

		int length = th2.size();		

		double M=0.0;
		double N= 0.0; 	
		double scal=1/b;	

			for (int i=0; i < length; ++i) {
				
				M+= exp(-0.5*sqr(th1[i]-s)/sqr(scal));
				N+= exp(-0.5*sqr(th1[i]-s)/sqr(scal))*th2[i];
			} 					
			
				N*=(1/M);
				//N*=(1/length);
			
		return N; 

}


void ImplicitPDESystem::Init_Rad(double rad1) { 

		double *r1p;
		r1p=&h2;
		*r1p= rad1;
				
}	

void ImplicitPDESystem::Init_Corr_Len(double l1, double l2) { 

		double *l1p,*l2p;
		l1p=&l_c1;
		*l1p= l1;
		
		l2p=&l_c2;
		*l2p=l2;
				
}	

void ImplicitPDESystem::Init_Sigma(double sigma1, double sigma2) { 

		double *sigma1p, *sigma2p;
		sigma1p=&sigma_new1;
		*sigma1p= sigma1;
		
		sigma2p=&sigma_new2;
		*sigma2p=sigma2;
		
}	


void ImplicitPDESystem::Init_Coords(){

	X_Coord=Vec::Zero(m->VertexSize());	
	Y_Coord=Vec::Zero(m->VertexSize());
	Z_Coord=Vec::Zero(m->VertexSize());


auto vi = m->VertexBegin(), ve = m->VertexEnd();


	for (; vi!=ve; ++vi) {
		
	        double x = (*vi)->x();
	        double y = (*vi)->y();
					double z = (*vi)->z();
					

					int idx = (*vi)->Index_ref();
					X_Coord[idx]=x;
					Y_Coord[idx]=y;
					Z_Coord[idx]=z;

				}

}


void ImplicitPDESystem::SetBC() {
	
	Mesh::VertexIt vi = m->BdryVertexBegin(), ve = m->BdryVertexEnd();
	bool fix1 = false, fix2 = false;
	for (; vi!=ve; ++vi ) {
		double z = (*vi)->z();
		if (z<bottom+bottom_var) {
			if ( (!fix2) & fix1 ) {
				(*vi)->MarkDirBdry1();
				fix2 = true;
			}
			if ( !fix1 ) {			
				(*vi)->MarkDirBdry1();
				(*vi)->MarkDirBdry2();
				fix1 = true;
			}
			//(*vi)->MarkDirBdry1();
			//(*vi)->MarkDirBdry2();
			(*vi)->MarkDirBdry3();
		}
		if (z>top-top_var) {
			(*vi)->MarkDirBdry3();
		}
	}
	
	Mesh::FaceIt fi = m->BdryFaceBegin(), fe = m->BdryFaceEnd();
	for (; fi!=fe; ++fi) 
		for (int j=0; j<3; ++j)
			if (!(*fi)->a()->DirBdry(j) || !(*fi)->b()->DirBdry(j) || !(*fi)->c()->DirBdry(j))
				(*fi)->MarkNeuBdry(j);
	
}

Vec3 ImplicitPDESystem::NeumannBC(Vec3 coord) const {
	return Vec3(0.0,0.0,0.0);
}

void ImplicitPDESystem::InitBC() {
	
	Mesh::VertexIt vi = m->BdryVertexBegin(), ve = m->BdryVertexEnd();
	for (; vi!=ve; ++vi ) {
		double z = (*vi)->z();
		if (z<bottom+2.0*bottom_var) {
			if ((*vi)->DirBdry1()) (*vi)->u1() = 0.0;
			if ((*vi)->DirBdry2()) (*vi)->u2() = 0.0;
			if ((*vi)->DirBdry3()) (*vi)->u3() = 0.0;
		} else {
			if ((*vi)->DirBdry1()) (*vi)->u1() = 0.0;
			if ((*vi)->DirBdry2()) (*vi)->u2() = 0.0;
			if ((*vi)->DirBdry3()) (*vi)->u3() = -0.5;
		}
		
	}
}

void ImplicitPDESystem::Init(std::string config_file) {
	
	//assert(config_file!="");
	// first thing, if we are given a file to init from, we load this data
	LoadOptions(config_file);
	
	Mesh::VertexIt vi = m->VertexBegin(), ve = m->VertexEnd();
	Vec3 bl = (*vi)->Coord(), tr = (*vi)->Coord();
	for (; vi!=ve; ++vi) {
		Vec3 coord = (*vi)->Coord();
		for (int j=0; j<3; ++j) {
			bl[j] = std::min(bl[j], coord[j]);
			tr[j] = std::max(tr[j], coord[j]);
		}
	}
	std::cerr << "Top-right  : " << tr[0] << ", " << tr[1] << ", " << tr[2] << "." << std::endl;
	std::cerr << "Bottom-left: " << bl[0] << ", " << bl[1] << ", " << bl[2] << "." << std::endl;
	std::cerr << "Difference : " << tr[0]-bl[0] << ", " << tr[1]-bl[1] << ", " << tr[2]-bl[2] << "." << std::endl;
	std::cerr << "Area: " << (area=sqr(0.25*((tr[0]-bl[0])+(tr[1]-bl[1])))*pi) << "." << std::endl;
	
	
	// lets set bc indicators
	top = tr[2];
	bottom = bl[2];
	len = top-bottom;
	std::cerr << "len is:" << len << std::endl;
	bottom_var *= len;
	top_var *= len;
	SetBC(); // we need boundary conditions.
	m_size=20; // for computation of the number of z-elements = 5*m_size

	Coords=Eigen::VectorXd::LinSpaced(m_size,bottom,top);

	// now we can index...
	std::cerr << "Indexing " << m->VertexSize() << " vertices, ";
	vi = m->VertexBegin();
	N_dof = 0;
	// count basis functions
	for (; vi!=ve; ++vi)
		for (int j=0; j<3; ++j) if ( !(*vi)->DirBdry(j) ) ++N_dof;
	
	idx_mutex = std::vector<std::mutex>(N_dof);
	std::cerr << N_dof << " interior degrees of freedom, " << 3*m->VertexSize()-N_dof << " Dirichlet boundary... ";
    
	int idx = 0, idx_bdry = 0,idx_ref=0;
	vi = m->VertexBegin();
	for (; vi!=ve; ++vi) {
				(*vi)->SetIndex_ref(idx_ref);
			++idx_ref;		
		for (int j=0; j<3; ++j) {
	if ( !(*vi)->DirBdry(j) ) {
		(*vi)->SetIndex(j, idx);
		++idx;
	} else 	{
		(*vi)->SetIndex(j, idx_bdry);
		++idx_bdry;
	 }
 }
}
	std::cerr << "done; ";
	
	
	std::cerr << "Preparing " << m->TetSize() << " Tets... ";
	N_tet = 0;
	Mesh::TetIt ti = m->TetBegin(), te = m->TetEnd();
	for (; ti!=te; ++ti) {
		(*ti)->Init_vol();
		(*ti)->Init_grad_s();
		mesh_sizes.push_back(	(*ti)-> Init_mesh_size());
		++N_tet;
	}

		
	size_t elpth = ceil( (double)N_tet/numThreads);
		
	std::cerr << "Preparing " << numThreads << " threads... ";
	ti = m->TetBegin();
	for ( int i=0; i<numThreads; ++i ) {
		std::pair< Mesh::TetIt, Mesh::TetIt > bounds;
		bounds.first = ti;
		int j = 0;
		while (j<elpth && ti!=te) {++j; ++ti;}
		bounds.second = ti;
		tet_th.push_back(bounds);
	}

#ifdef EIGEN_HAS_OPENMP
	std::cerr << "Eigen is running in parallel... ";
	Eigen::setNbThreads(numThreads);
	Eigen::initParallel();
	omp_set_num_threads(numThreads); 
#endif 
	std::cerr << "ok; ";
	

	std::cerr << "zustandsvektor... ";
	Z = Vec::Zero(N_dof);
	U_Db = Vec::Zero(3*m->VertexSize()-N_dof);
	FNeu = Vec::Zero(N_dof);
	
	vi = m->VertexBegin();


	for (; vi!=ve; ++vi)
		{
		
		for (int j=0; j<3; ++j) 
	{
			if (!(*vi)->DirBdry(j)) (*vi)->Attach_u(j, &Z[(*vi)->Index(j)]);
	else (*vi)->Attach_u(j, &U_Db[(*vi)->Index(j)] );

	}
}
	
	std::cerr << "ok; ";
	InitBC();
	
}


void ImplicitPDESystem::Perform_Cholesky() {



	Mesh::VertexIt vi = m->VertexBegin(), ve = m->VertexEnd();
	Vec3 bl = (*vi)->Coord(), tr = (*vi)->Coord();
	for (; vi!=ve; ++vi) {
		Vec3 coord = (*vi)->Coord();
	}
		// now we can index...
	std::cerr << "Indexing " << m->VertexSize() << " vertices, ";
	vi = m->VertexBegin();
	N_dof = 0;
	// count basis functions
	for (; vi!=ve; ++vi)
		{
		for (int j=0; j<3; ++j) 
			{

					if ( !(*vi)->DirBdry(j) ) 
					{
							++N_dof;
					}
		  }
	 }	
	idx_mutex = std::vector<std::mutex>(N_dof);
	std::cerr << N_dof << " interior degrees of freedom, " << 3*m->VertexSize()-N_dof << " Dirichlet boundary... ";
    
	int idx = 0, idx_bdry = 0,idx_ref=0;
	vi = m->VertexBegin();
	for (; vi!=ve; ++vi) {
				(*vi)->SetIndex_ref(idx_ref);
			++idx_ref;		
		for (int j=0; j<3; ++j) {
	if ( !(*vi)->DirBdry(j) ) {
		(*vi)->SetIndex(j, idx);
		++idx;
	} else 	{
		(*vi)->SetIndex(j, idx_bdry);
		++idx_bdry;
	 }
 }
}
	std::cerr << "done; ";

	X1 = Vec::Zero(m->VertexSize());
 	X2 = Vec::Zero(m->VertexSize());
	X3 = Vec::Zero(m->VertexSize());
	Init_Coords();	
	
#ifdef EIGEN_HAS_OPENMP
	std::cerr << "Eigen is running in parallel... ";
	Eigen::setNbThreads(numThreads);
	Eigen::initParallel(); 
#endif 

std::cerr<< "Performing Cholesky-Decompostion" << std::endl;
CholeskyLoop();
std::cerr<< "Did it work?" << std::endl;

// Performing Cholesky Decomposition for geometric perturbations; the decomposition needs to be performed only once
	for (int i=1; i<2; ++i){
	
	LLT<Mat> lltOfA(KTL_new1);
	KTL_Cholnew1= lltOfA.matrixL();  
	
	 	}
	 	
	  LLT<Mat> lltOfA(KTL_new2);
	  KTL_Cholnew2= lltOfA.matrixL();

std::cerr<< "It worked!" << std::endl;
}

void ImplicitPDESystem::CholeskyLoop () {

std::vector<T> vals_ktl1, vals_ktl2;
double x1,y1,z1;
double x2,y2,z2;
KTL_new1=Mat(5*m_size,5*m_size);
KTL_new2=Mat(5*m_size,5*m_size);
double b1= 5*m_size-1;	


// Computing entries of the covariance matrix for geomertic perturbations
	for (int i=0; i < 5*m_size; ++i) {

			double b=i;
			z1=(b/b1);
					
				for (int j=0;j < 5*m_size; ++j) {
				
					double c=j;
					z2= c/b1;
				
				double ktl1 = exp(-(abs(z1-z2)/l_c1));
				double ktl2 = exp(-(abs(z1-z2)/l_c2));
								
				KTL_new1(i,j)=ktl1;
				KTL_new2(i,j)=ktl2;
							
					}
											

				}	
		
}


void ImplicitPDESystem::Perform_Perturbation() {

		std::vector<double> X,Y,Z,XYZ,Knots, HKLT,XX,YY,XX1,XX2;
		//std::vector<double> X(m_size),Y(m_size),Z(m_size),Knots(m_size);
		Vec eins =Vec::Ones(m-> VertexSize());	
		std::cerr<< "Performing Geometric Perturbations..";
		Prep_Norm();
		Vec KTL_Chol1;
		Vec KTL_Chol2;
		
		KTL_Chol1= sigma_new1*KTL_Cholnew1*C_norm1; 
		KTL_Chol2 = sigma_new2*KTL_Cholnew1*C_norm2;
		
		double b22=5*m_size-1;

		for (int i=0; i<5*m_size; ++i){
		
				double b11=i;
				X.push_back(KTL_Chol1[i]);
				Y.push_back(KTL_Chol2[i]);
				Knots.push_back(b11/b22);
	
			}
		
	
		X1 = Vec::Zero(5*m_size);
	  	X2 = Vec::Zero(5*m_size);

		double ll2=5*m_size-1;		
	
 for (int i=0; i< 5*m_size; ++i){ 
 	
 	X1[i]+= (X[i]-X[0]);
 	X2[i]+= (Y[i]-Y[0]);
 	
 	XX1.push_back(X[i]-X[0]);
 	XX2.push_back(Y[i]-Y[0]);

 }
 
	tk::spline spline1,spline2,spline3,spline4; 
  	spline1.set_points(Knots,XX1);
 	spline2.set_points(Knots,XX2);
 
		X1_new = Vec::Zero(m-> VertexSize());
	  	X2_new = Vec::Zero(m-> VertexSize());
		X3_new = Vec::Zero(m-> VertexSize());
 
  // Perform cubic spline interpolation; currently it is required that X is already sorted
 for (int i=0; i <  m-> VertexSize(); ++i){
 
	 X1_new[i] = X_Coord[i]+ spline1(Z_Coord[i]);
	 X2_new[i] = Y_Coord[i]+ spline2(Z_Coord[i]);
 
 } 

int idx=0;
Mesh::VertexIt vi = m->VertexBegin(), ve = m->VertexEnd();
double ll=5*m_size-1;

// change the values of the coordinate values; only for plane coordinates x1 and x2
	for (; vi!=ve; ++vi)
		{
		
		double x1=(*vi)-> m_coord.x();
		double x2=(*vi)-> m_coord.y();
		double * point;
		int j = floor(ll*((*vi)-> m_coord.z()));
		point = &((*vi)-> m_coord.x());
		*point= X1_new[idx];

		point = &((*vi)-> m_coord.y());
		*point= X2_new[idx];
	
		++idx;
		
		}
}


void ImplicitPDESystem::Assemble() {
	std::cerr << "K... ";
	PrepK();
	std::cerr << "ok; "; 
}


// Generate random numbers V^i for computing density variation and geometric perturbation
// via L*V^i using Cholesk decomposition
void ImplicitPDESystem::Prep_Norm() {

			double trans1, trans2,trans3,trans4;	

			C_norm1=Vec::Zero(5*m_size);
			C_norm2=Vec::Zero(5*m_size);

	for (int i=0; i < 5*m_size; ++i) {

					// numbers for geomertic perturbation
					trans1 = distribution_norm1(generator_norm1);
					trans2 = distribution_norm1(generator_norm2);
					
					C_norm1[i]+= trans1; 
					C_norm2[i]+= trans2;
													
					}
}
	
		

void ImplicitPDESystem::KLoop ( SpMat* th_k, Vec* th_b, Mesh::TetIt ti, Mesh::TetIt te ) {
	
	std::vector<T> vals_k;
	*th_b = Vec::Zero(N_dof);
	
	Vec3 grad_sj, grad_si;
	
	for (; ti!=te; ++ti) {
		for (int bj=0; bj<4; ++bj) {
			
			(*ti)->Calc_grad_s(bj, grad_sj);
			grad_sj[0]*=(1/h2);
			grad_sj[1]*=(1/h2);
			Vertex* vj = (*ti)->v(static_cast<VertexName>(bj));
			
			for (int j=0; j<3; ++j) {
				int idx_j = vj->Index(j);
				
				if (!vj->DirBdry(j)) {
					for (int bi=0; bi<4; ++bi) {
						
						(*ti)->Calc_grad_s(bi, grad_si);
						Vertex* vi = (*ti)->v(static_cast<VertexName>(bi));
						grad_si[0]*=(1/h2);
						grad_si[1]*=(1/h2);
						
						for (int i=0; i<3; ++i) {
							if (vi->DirBdry(i)) continue;
							int idx_i = vi->Index(i);
							//if (idx_i < idx_j) continue;
			
							double k = 0.0;
					
							for (int l=0; l<3; ++l) {
								for (int m=0; m<3; ++m) {
									k += C(j,l,i,m)*grad_sj[l]*grad_si[m] * (*ti)->Vol();
								}
							}
							
							//k += j==i ? dot(grad_sj,grad_si) * (*ti)->Vol() : 0.0;
							vals_k.push_back( T(idx_i, idx_j, k) );
						}
					}
				} else {
					double c_j = vj->u(j);
					for (int bi=0; bi<4; ++bi) {
						
						(*ti)->Calc_grad_s(bi, grad_si);
						grad_si[0]*=(1/h2);
						grad_si[1]*=(1/h2);
						Vertex* vi = (*ti)->v(static_cast<VertexName>(bi));
						
						for (int i=0; i<3; ++i) {
						
							if (vi->DirBdry(i)) continue;
							int idx_i = vi->Index(i);
					
							double b = 0.0;
							
							for (int l=0; l<3; ++l) {
								for (int m=0; m<3; ++m) {
									b -= c_j*C(j,l,i,m)*grad_sj[l]*grad_si[m] * (*ti)->Vol();
								}
							}
							
							(*th_b)[ idx_i ] += b;
						}
					
					}
				}
			}
		
		}
	}

	(*th_k).setFromTriplets( vals_k.begin(), vals_k.end() );	
	
}

void ImplicitPDESystem::PrepK() {
	std::vector<SpMat> th_k;
	std::vector<Vec> th_b;

	for ( int j = 0; j < numThreads; ++j ) {
		th_k.push_back( SpMat(N_dof, N_dof) );
		th_b.push_back( Vec() );
	}
  
	std::vector<std::thread> threads;
	for ( int j = 0; j < numThreads; ++j ) 
		threads.push_back( std::thread( &ImplicitPDESystem::KLoop, this, 
	&th_k[j], &th_b[j], tet_th[j].first, tet_th[j].second)  );
	// join threads
	for (auto &thread : threads) thread.join();
	K = SpMat(N_dof,N_dof);
	B = Vec::Zero(N_dof);
	for (int j = 0; j < numThreads; ++j) {
		K += th_k[j];
		B += th_b[j];
	}

	K.makeCompressed();
}

void ImplicitPDESystem::CalcFNeu() {
	
	Mesh::FaceIt fi = m->BdryFaceBegin(), fe = m->BdryFaceEnd();
	for (; fi!=fe; ++fi) {
		for (int j=0; j<3; ++j) {
			if (!(*fi)->NeuBdry(j)) continue;
			Vertex *va = (*fi)->a(), *vb = (*fi)->b(), *vc = (*fi)->c();
			Vec3 ca = va->Coord(), cb = vb->Coord(), cc = vc->Coord();
			Vec3 fa = NeumannBC(ca), fb = NeumannBC(cb), fc = NeumannBC(cc);
			Vec3 f = (*fi)->Area()*1.0/3.0*(fa+fb+fc);
			
			if (!(va->DirBdry(j))) FNeu[va->Index(j)] += 1.0/3.0 * f[j];
			if (!(vb->DirBdry(j))) FNeu[vb->Index(j)] += 1.0/3.0 * f[j];
			if (!(vc->DirBdry(j))) FNeu[vc->Index(j)] += 1.0/3.0 * f[j];
			
			
			
		}
		
	}
	
	
}


void ImplicitPDESystem::SolveWithGuess() {
	
	std::cerr << "Calc Neumann-force...";
	CalcFNeu();
	B += FNeu;
	std::cerr << "done." << std::endl;
	double t1,t2; 
	t1= omp_get_wtime();
	std::cerr << "Solver... ";
	Eigen::ConjugateGradient<SpMat, Eigen::Lower | Eigen::Upper> cg;
	cg.setTolerance(1e-5);
	cg.compute(K);
	
	std::cerr << "Solving... ";
	Z = cg.solveWithGuess(B,Z0);
	std::cerr << "done." << std::endl;

	t2 = omp_get_wtime();
	std::cout<< "time needed=" << t2-t1<< std::endl; 
	
}

void ImplicitPDESystem::Solve() {
	
	std::cerr << "Calc Neumann-force...";
	CalcFNeu();
	B += FNeu;
	std::cerr << "done." << std::endl;
	double t1,t2; 
	t1= omp_get_wtime();
	std::cerr << "Solver... ";
	Eigen::ConjugateGradient<SpMat, Eigen::Lower | Eigen::Upper> cg;
	cg.setTolerance(1e-5);
	cg.compute(K);
	
	std::cerr << "Solving... ";
	Z = cg.solve(B);
	std::cerr << "done." << std::endl;

	t2 = omp_get_wtime();
	std::cout<< "time needed="<< t2-t1<< std::endl; 

	Z0=Vec::Zero(N_dof);	
	Z0=Z;
	
}


void ImplicitPDESystem::ELoop ( double* th_e, Mesh::TetIt ti, Mesh::TetIt te ) {
	*th_e = 0;
	Vec3 grad_sj, grad_si;
	
	for (; ti!=te; ++ti) {
		for (int bj=0; bj<4; ++bj) {
			
			(*ti)->Calc_grad_s(bj, grad_sj);
			grad_sj[0]*=(1/h2);
			grad_sj[1]*=(1/h2);
			Vertex* vj = (*ti)->v(static_cast<VertexName>(bj));
			
			for (int j=0; j<3; ++j) {
			
					for (int bi=0; bi<4; ++bi) {
						
						(*ti)->Calc_grad_s(bi, grad_si);
						grad_si[0]*=(1/h2);
						grad_si[1]*=(1/h2);
						Vertex* vi = (*ti)->v(static_cast<VertexName>(bi));
						
						for (int i=0; i<3; ++i) {
							
							for (int l=0; l<3; ++l) {
								for (int m=0; m<3; ++m) {
									*th_e += 0.5*C(j,l,i,m)*grad_sj[l]*grad_si[m]*vi->u(i)*vj->u(j) * (*ti)->Vol();
								}
							}
							
						}
					}

			}
		
		}
	}
}

// compute elastic energy
double ImplicitPDESystem::CalcE() {
	std::vector<double> th_e;
	for ( int j = 0; j < numThreads; ++j ) {
		th_e.push_back( 0.0 );
	}
  
	std::vector<std::thread> threads;
	for ( int j = 0; j < numThreads; ++j ) 
		threads.push_back( std::thread( &ImplicitPDESystem::ELoop, this, 
	&th_e[j], tet_th[j].first, tet_th[j].second)  );
	// join threads
	for (auto &thread : threads) thread.join();
  
	double  e = 0.0;
	for (int j = 0; j < numThreads; ++j) {
		e += th_e[j];
	}
	return e;
}

// compute effective Young's modulus
void ImplicitPDESystem::CalcYoungs(int k) {
	
	double e = CalcE();
	double test=0.0;

	std::cerr << "Energy: " << e << "." << std::endl;
	std::cerr << "Effective Young's modulus: " << 2.0*len*e/(area*sqr(0.5)) << "." << std::endl;

}

// load options regarding shear and bulk modulus etc.
void ImplicitPDESystem::LoadOptions(std::string fname) {
	std::ifstream is;
	is.open(fname.c_str());

	double K, mu;
	
	std::string line;
	while( std::getline(is, line) ){
		
		boost::char_separator<char> sep(" ");
		boost::tokenizer< boost::char_separator<char> > tok(line, sep);
		boost::tokenizer< boost::char_separator<char> >::iterator ti = tok.begin(), te = tok.end();
		
		if (ti==te) continue; //empty line.
		
		if ( *ti == "th" ) numThreads = std::stoi(*(++ti));
		if ( *ti == "K" ) K = std::stod(*(++ti));
		if ( *ti == "mu" ) mu = std::stod(*(++ti));
		if ( *ti == "t_var" ) top_var = std::stod(*(++ti));
		if ( *ti == "b_var" ) bottom_var = std::stod(*(++ti));
	}
	C = ElMod(K,mu);
}



