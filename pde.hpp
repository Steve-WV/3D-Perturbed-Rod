#ifndef __PDE_HPP__
#define __PDE_HPP__

#include <string>
#include "mesh.hpp"

#include <thread>
#include <mutex>

#include <Eigen/Sparse>
#include<Eigen/SparseCholesky>
#include <Eigen/Core>
#include <unsupported/Eigen/Splines>

struct ElMod {
	
	inline ElMod() {}
	inline ElMod(double, double);
	
	inline double& operator ()(int i, int j, int k, int l) { return C[i][j][k][l]; }
	inline double operator ()(int i, int j, int k, int l) const { return C[i][j][k][l]; }
	
	double C[3][3][3][3];

};

class ImplicitPDESystem {
public:
	
    typedef Eigen::VectorXd Vec;
   // typedef Eigen::SparseMatrix<double, Eigen::RowMajor> SpMat;
    typedef Eigen::SparseMatrix<double, Eigen::ColMajor> SpMat;
    typedef Eigen::MatrixXd Mat;
    typedef Eigen::Triplet<double> T;
    typedef Eigen::SimplicialLDLT<SpMat> SpSolver; 
    typedef Eigen::Spline<double,1> Spline1d;
    //typedef Eigen::SplineFitting<Spline1d>::Spline1dFitting;
     inline ImplicitPDESystem( Mesh* mesh ) : m(mesh) { }
    
  void Init(std::string);
  void Init_Coords();
  void Assemble();
  void Perform_Perturbation();
  void Perform_Cholesky();
  void Solve();
  void SolveWithGuess();
  void CalcYoungs(int);
  void Export_Vol(std::string) const;
  void Export_Surface(std::string) const;
  void SaveFace(std::string) const;
  void Init_Corr_Len(double,double);
  void Init_Sigma(double,double);
  void Init_Rad(double);
   
  
  
private:
	void CalcFNeu();
  	void SetBC();
	void InitBC();
	Vec3 NeumannBC(Vec3) const; 
	void LoadOptions(std::string);
  	void PrepK();
  	void KLoop( SpMat*, Vec*, Mesh::TetIt, Mesh::TetIt);
	void Prep_Norm();
	void CholeskyLoop();
	
 
  void ELoop ( double*, Mesh::TetIt, Mesh::TetIt );
  double CalcE();

  // Access to the basis functions, etc.
  Mesh* m;
    
  // basis functions (that are not constrained) and elements
  int N_dof, N_tet;
  
  // mutexes and threading
  mutable std::vector<std::mutex> idx_mutex;
  mutable std::mutex mut;
  std::vector<std::pair<Mesh::TetIt, Mesh::TetIt> > tet_th;
 
  // general params
  int numThreads=4; // number of threads
  //int numThreads=1; // number of threads
  
  ElMod C;
  double bottom, top, bottom_var, top_var;
  double area, len, kernel_b;
  int bdry_size;
  int m_size, num_ev,num_ev1;
  double l_c1,l_c2,sigma_new1,sigma_new2,h2;
  
  //minimization data
  Vec B, FNeu;
  Vec Z, U_Db,X1,X2,X3,X_Coord,Y_Coord,Z_Coord,Z0,C_norm1,C_norm2,C_norm3,C_norm4,Coords,KTL_eig,C_norm5,C_norm6;
  Vec X1_new,X2_new,X3_new;
  Vec X1_bdry,X2_bdry,X3_bdry,X_Coord_bdry,Y_Coord_bdry,Z_Coord_bdry;
  
  SpMat K,KTL,Test_Matrix,SIG1,SIG2,SIG3;
  Mat KTL_Chol,KTL_Choleig;
  Mat KTL1,KTL_new1, KTL_new2,KTL_Cholnew1,KTL_Cholnew2 ;

  std::vector<double> mesh_sizes; 
  //std::vector<double> long_clear;
  //int long_len;

SpSolver solver;

};



inline ElMod::ElMod(double K, double mu) {
	C[0][0][0][0] = K+mu+mu-2.0/3.0*mu;
	C[0][0][0][1] = 0.0;
	C[0][0][0][2] = 0.0;
	C[0][0][1][0] = 0.0;
	C[0][0][1][1] = K-2.0/3.0*mu;
	C[0][0][1][2] = 0.0;
	C[0][0][2][0] = 0.0;
	C[0][0][2][1] = 0.0;
	C[0][0][2][2] = K-2.0/3.0*mu;
	
	C[0][1][0][0] = 0.0;
	C[0][1][0][1] = mu;
	C[0][1][0][2] = 0.0;
	C[0][1][1][0] = mu;
	C[0][1][1][1] = 0.0;
	C[0][1][1][2] = 0.0;
	C[0][1][2][0] = 0.0;
	C[0][1][2][1] = 0.0;
	C[0][1][2][2] = 0.0;
	
	C[0][2][0][0] = 0.0;
	C[0][2][0][1] = 0.0;
	C[0][2][0][2] = mu;
	C[0][2][1][0] = 0.0;
	C[0][2][1][1] = 0.0;
	C[0][2][1][2] = 0.0;
	C[0][2][2][0] = mu;
	C[0][2][2][1] = 0.0;
	C[0][2][2][2] = 0.0;
	
	C[1][0][0][0] = 0.0;
	C[1][0][0][1] = mu;
	C[1][0][0][2] = 0.0;
	C[1][0][1][0] = mu;
	C[1][0][1][1] = 0.0;
	C[1][0][1][2] = 0.0;
	C[1][0][2][0] = 0.0;
	C[1][0][2][1] = 0.0;
	C[1][0][2][2] = 0.0;
	
	C[1][1][0][0] = K-2.0/3.0*mu;
	C[1][1][0][1] = 0.0;
	C[1][1][0][2] = 0.0;
	C[1][1][1][0] = 0.0;
	C[1][1][1][1] = K+mu+mu-2.0/3.0*mu;
	C[1][1][1][2] = 0.0;
	C[1][1][2][0] = 0.0;
	C[1][1][2][1] = 0.0;
	C[1][1][2][2] = K-2.0/3.0*mu;
	
	C[1][2][0][0] = 0.0;
	C[1][2][0][1] = 0.0;
	C[1][2][0][2] = 0.0;
	C[1][2][1][0] = 0.0;
	C[1][2][1][1] = 0.0;
	C[1][2][1][2] = mu;
	C[1][2][2][0] = 0.0;
	C[1][2][2][1] = mu;
	C[1][2][2][2] = 0.0;
	
	C[2][0][0][0] = 0.0;
	C[2][0][0][1] = 0.0;
	C[2][0][0][2] = mu;
	C[2][0][1][0] = 0.0;
	C[2][0][1][1] = 0.0;
	C[2][0][1][2] = 0.0;
	C[2][0][2][0] = mu;
	C[2][0][2][1] = 0.0;
	C[2][0][2][2] = 0.0;
	
	C[2][1][0][0] = 0.0;
	C[2][1][0][1] = 0.0;
	C[2][1][0][2] = 0.0;
	C[2][1][1][0] = 0.0;
	C[2][1][1][1] = 0.0;
	C[2][1][1][2] = mu;
	C[2][1][2][0] = 0.0;
	C[2][1][2][1] = mu;
	C[2][1][2][2] = 0.0;
	
	C[2][2][0][0] = K-2.0/3.0*mu;
	C[2][2][0][1] = 0.0;
	C[2][2][0][2] = 0.0;
	C[2][2][1][0] = 0.0;
	C[2][2][1][1] = K-2.0/3.0*mu;
	C[2][2][1][2] = 0.0;
	C[2][2][2][0] = 0.0;
	C[2][2][2][1] = 0.0;
	C[2][2][2][2] = K+mu+mu-2.0/3.0*mu;
	
	/*	C[0][0][0][0] = 1.0;
		C[0][0][0][1] = 0.0;
		C[0][0][0][2] = 0.0;
		C[0][0][1][0] = 0.0;
		C[0][0][1][1] = 0.0;
		C[0][0][1][2] = 0.0;
		C[0][0][2][0] = 0.0;
		C[0][0][2][1] = 0.0;
		C[0][0][2][2] = 0.0;
	
		C[0][1][0][0] = 0.0;
		C[0][1][0][1] = 1.0;
		C[0][1][0][2] = 0.0;
		C[0][1][1][0] = 0.0;
		C[0][1][1][1] = 0.0;
		C[0][1][1][2] = 0.0;
		C[0][1][2][0] = 0.0;
		C[0][1][2][1] = 0.0;
		C[0][1][2][2] = 0.0;
	
		C[0][2][0][0] = 0.0;
		C[0][2][0][1] = 0.0;
		C[0][2][0][2] = 1.0;
		C[0][2][1][0] = 0.0;
		C[0][2][1][1] = 0.0;
		C[0][2][1][2] = 0.0;
		C[0][2][2][0] = 0.0;
		C[0][2][2][1] = 0.0;
		C[0][2][2][2] = 0.0;
	
		C[1][0][0][0] = 0.0;
		C[1][0][0][1] = 0.0;
		C[1][0][0][2] = 0.0;
		C[1][0][1][0] = 1.0;
		C[1][0][1][1] = 0.0;
		C[1][0][1][2] = 0.0;
		C[1][0][2][0] = 0.0;
		C[1][0][2][1] = 0.0;
		C[1][0][2][2] = 0.0;
	
		C[1][1][0][0] = 0.0;
		C[1][1][0][1] = 0.0;
		C[1][1][0][2] = 0.0;
		C[1][1][1][0] = 0.0;
		C[1][1][1][1] = 1.0;
		C[1][1][1][2] = 0.0;
		C[1][1][2][0] = 0.0;
		C[1][1][2][1] = 0.0;
		C[1][1][2][2] = 0.0;
	
		C[1][2][0][0] = 0.0;
		C[1][2][0][1] = 0.0;
		C[1][2][0][2] = 0.0;
		C[1][2][1][0] = 0.0;
		C[1][2][1][1] = 0.0;
		C[1][2][1][2] = 1.0;
		C[1][2][2][0] = 0.0;
		C[1][2][2][1] = 0.0;
		C[1][2][2][2] = 0.0;
	
		C[2][0][0][0] = 0.0;
		C[2][0][0][1] = 0.0;
		C[2][0][0][2] = 0.0;
		C[2][0][1][0] = 0.0;
		C[2][0][1][1] = 0.0;
		C[2][0][1][2] = 0.0;
		C[2][0][2][0] = 1.0;
		C[2][0][2][1] = 0.0;
		C[2][0][2][2] = 0.0;
	
		C[2][1][0][0] = 0.0;
		C[2][1][0][1] = 0.0;
		C[2][1][0][2] = 0.0;
		C[2][1][1][0] = 0.0;
		C[2][1][1][1] = 0.0;
		C[2][1][1][2] = 0.0;
		C[2][1][2][0] = 0.0;
		C[2][1][2][1] = 1.0;
		C[2][1][2][2] = 0.0;
	
		C[2][2][0][0] = 0.0;
		C[2][2][0][1] = 0.0;
		C[2][2][0][2] = 0.0;
		C[2][2][1][0] = 0.0;
		C[2][2][1][1] = 0.0;
		C[2][2][1][2] = 0.0;
		C[2][2][2][0] = 0.0;
		C[2][2][2][1] = 0.0;
		C[2][2][2][2] = 1.0; */
	
}

#endif
