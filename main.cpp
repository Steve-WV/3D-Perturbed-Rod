#include "mesh.hpp"
#include "vertex.hpp"
#include "tet.hpp"
#include "pde.hpp"
#include <iomanip>
#include <getopt.h>
#include <iostream>
#include <fstream>
#include <cmath>

#include "tmv.hpp"
const double h1= 1.0;
const double eps1=0.05;
const double eps2=0.05;

const double scaler=0.3;
//const double scaler=0.2;
//const double scaler=0.4;

const double radi[5]= {0.25,0.1,0.05,0.025,0.01};

//const double vece01[9]={0.00144, 0.00335, 0.001663, 0.00309, 0.00188, 0.00284, 0.00211, 0.00258, 0.00233};
//const double vece02[9]={0.00209, 0.00454, 0.00233, 0.00416, 0.00257, 0.00380, 0.00281, 0.00342, 0.00305};

int main(int argc, char *argv[]) {
	
	std::string obj_fname, opt_fname="", surf_fname;
	bool exp_surf = false, calc_vf = false;
	double crop_vf = 0.0;
	int c;
	while ((c = getopt (argc, argv, "s:f:c:v:")) != -1)
		switch (c)
	{

		case 'f':
		obj_fname = std::string(optarg);
		break;
		case 'c':
		opt_fname = std::string(optarg);
		break;
		case 's':
		surf_fname = std::string(optarg);
		exp_surf = true;
		break;
		case 'v':
		calc_vf = true;
		crop_vf = std::stod(optarg);
		break;
		case '?':
		if (optopt == 'f' || optopt == 'c')
			std::cerr << "Option -" << optopt <<" requires an argument.\n";
		else if (isprint (optopt))
			std::cerr << "Unknown option `-" << optopt << "'.\n";
		else
			std::cerr << "Unknown option character `" << optopt << "'.\n";
		return 1;
		default:
		return(-1);
	}
    
	if (argc-optind != 0 ||  obj_fname == "") {
		std::cerr << "Problem with options.";
		return(-1);
	}
	
			Mesh* m = new Mesh();
			m->Load(obj_fname);
			ImplicitPDESystem* pde = new ImplicitPDESystem(m);	

	for (int k=0; k<5; ++k){
	
			    pde-> Init_Corr_Len(eps1,eps1); 
			    pde -> Init_Rad(radi[k]);
			    pde-> Init_Sigma(eps2*scaler,eps2*scaler);
			    pde->Init(opt_fname);
			    
			    // Perform the Cholesky decomposition only once
			    pde->Perform_Cholesky();
			    
			    pde->Assemble();
			    pde->Solve();
			    pde->CalcYoungs(k);

		// Calculating Number samples of deformations
		int Numb = 500;	

		for (int i=1; i<Numb; ++i)
		{
	
	           if (calc_vf) {
		
		  std::cerr << "Calculating Volume: ";
		  auto vi = m->VertexBegin(), ve = m->VertexEnd();
		  Vec3 bl = (*vi)->Coord(), tr = (*vi)->Coord();
		  for (; vi!=ve; ++vi) {
			Vec3 coord = (*vi)->Coord();
			for (int j=0; j<3; ++j) {
				bl[j] = std::min(bl[j], coord[j]);
				tr[j] = std::max(tr[j], coord[j]);
			   }
			}
			Vec3 bl_crop = bl + crop_vf * (tr-bl);
			Vec3 tr_crop = tr - crop_vf * (tr-bl);
			
			double vol = 0.0;
			auto ti = m->TetBegin(), te = m->TetEnd();
			for (; ti!=te; ++ti) {
			    Vec3 a = (*ti)->a()->Coord(); 
				Vec3 b = (*ti)->b()->Coord();
				Vec3 c = (*ti)->c()->Coord(); 
				Vec3 d = (*ti)->d()->Coord();
				
				a = crop( a, bl_crop, tr_crop );
				b = crop( b, bl_crop, tr_crop );
				c = crop( c, bl_crop, tr_crop );
				d = crop( d, bl_crop, tr_crop );
    		
		    Vec3 DA = d - a;
    	            Vec3 DB = d - b;
		    Vec3 DC = d - c;
    
		    vol += 1.0/6.0 * dot( DA, cross(DB, DC));
		}
		
		std::cerr << vol << std::endl;
		Vec3 box = tr_crop-bl_crop;
		double box_size = box.x()*box.y()*box.z();
		std::cerr << "Box Size: " << box_size << std::endl;
		std::cerr << "Volume Fraction: " << vol/box_size << std::endl;
		
		std::ofstream os;
		os.open(obj_fname+"_volinfo.txt");
		os << std::setprecision(16);
		os << "Volum Information" << std::endl;
		os << "Crop: " << crop_vf << std::endl;
		os << "Box volume: " << box_size << std::endl;
		os << "Volume: " << vol << std::endl;
		os << "Volume Fraction: " << vol/box_size << std::endl;
		os.close();
		
	}
	
	if (exp_surf) pde->SaveFace(surf_fname);
	else {

//Perform geometric perturbations
		pde -> Perform_Perturbation();
		pde->Init(opt_fname);
		
		//Calculat effective Young's Modulus of the deformed cylinder
		pde->Assemble();
		pde->SolveWithGuess();
		pde->CalcYoungs(k);
		
	//pde->~ImplicitPDESystem();
	//m->~Mesh();

		}
	}
    
  }
	return(0);
}
