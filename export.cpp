#include <iomanip>
#include <fstream>
#include <unordered_map>
#include "pde.hpp"
#include "mesh.hpp"
#include "vertex.hpp"
#include "tmv.hpp"



#include <Eigen/Core>
#include <Eigen/SparseCore>

#ifdef HAVE_VTK
#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkTetra.h>
#include <vtkCellArray.h>
#include <vtkDataSetMapper.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPointData.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkFloatArray.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkPolyData.h>
#include <vtkPoints.h>

void ImplicitPDESystem::Export_Surface(std::string fname) const 
{ 	
	std::lock_guard<std::mutex> lock(mut);
	
	auto points = vtkSmartPointer<vtkPoints>::New();
	points->SetNumberOfPoints(m->BdryVertexSize());
	Mesh::const_VertexIt vi = m->BdryVertexBegin(), ve = m->BdryVertexEnd(); 
	float p[3]; 
	for (int j=0; vi!=ve; ++vi) {
		p[0] = (*vi)->x();
		p[1] = (*vi)->y();
		p[2] = (*vi)->z();
		points->SetPoint(j++, p);
	}
	
	// we need some indexes
	std::unordered_map<Vertex*, int> Vt_list_idx;
	vi = m->BdryVertexBegin(); 
	for (int j = 0; vi!=ve; ++vi) Vt_list_idx[*vi] = j++;
 	
	auto Conn = vtkSmartPointer<vtkIdTypeArray>::New();
	Conn->SetNumberOfValues(4*m->BdryFaceSize());
	
	auto fi = m->BdryFaceBegin(), fe = m->BdryFaceEnd(); 
	for (int j = 0; fi!=fe; ++fi) {
		Conn->SetValue( j++, 3 );
		Conn->SetValue( j++, Vt_list_idx[(*fi)->a()] );
		Conn->SetValue( j++, Vt_list_idx[(*fi)->b()] );
		Conn->SetValue( j++, Vt_list_idx[(*fi)->c()] );
	}

	auto Faces = vtkSmartPointer<vtkCellArray>::New();
	Faces->SetCells(m->BdryFaceSize(),Conn);
	
	auto poly = vtkSmartPointer<vtkPolyData>::New();
	poly->SetPoints(points);
	poly->SetPolys(Faces);
	
	auto u = vtkSmartPointer<vtkFloatArray>::New();
	u->SetNumberOfComponents(3);
	u->SetName("u");
	u->SetNumberOfValues(m->BdryVertexSize()*3);
	vi = m->BdryVertexBegin(); 
	for (int j = 0; vi!=ve; ++vi) u->SetTuple3(j++, (*vi)->u1(), (*vi)->u2(), (*vi)->u3());
	poly->GetPointData()->AddArray(u);
		
	auto DirBC = vtkSmartPointer<vtkFloatArray>::New();
	DirBC->SetNumberOfComponents(3);
	DirBC->SetName("DirBC");
	DirBC->SetNumberOfValues(m->BdryVertexSize()*3);
	vi = m->BdryVertexBegin(); 
	for (int j = 0; vi!=ve; ++vi) DirBC->SetTuple3(j++, 
	(*vi)->DirBdry(0) ? 1.0 : 0.0, 
	(*vi)->DirBdry(1) ? 1.0 : 0.0, 
	(*vi)->DirBdry(2) ? 1.0 : 0.0 );
	poly->GetPointData()->AddArray(DirBC);
		
	auto NeuBC = vtkSmartPointer<vtkFloatArray>::New();
	NeuBC->SetNumberOfComponents(3);
	NeuBC->SetName("NeuBC");
	NeuBC->SetNumberOfValues(m->BdryVertexSize()*3);
	vi = m->BdryVertexBegin(); 
	for (int j = 0; vi!=ve; ++vi) {
		Vec3 coord = (*vi)->Coord();
		NeuBC->SetTuple3(j++,
		NeumannBC(coord)[0],
		NeumannBC(coord)[1],
		NeumannBC(coord)[2] );
	}
	poly->GetPointData()->AddArray(NeuBC);
	
	// Write file
	auto writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
	writer->SetDataModeToAppended();
	writer->EncodeAppendedDataOff();
	std::string fn = fname+"_surf.vtp";
	writer->SetFileName(fn.c_str());
#if VTK_MAJOR_VERSION <= 5
	writer->SetInput(poly);
#else
	writer->SetInputData(poly);
#endif
	writer->Write();
}


void ImplicitPDESystem::Export_Vol(std::string fname) const 
{ 
	std::lock_guard<std::mutex> lock(mut);
	
	auto points = vtkSmartPointer<vtkPoints>::New();
	points->SetNumberOfPoints(m->VertexSize());
	Mesh::const_VertexIt vi = m->VertexBegin(), ve = m->VertexEnd(); 
	float p[3]; 
	for (int j=0; vi!=ve; ++vi) {
		p[0] = (*vi)->x();
		p[1] = (*vi)->y();
		p[2] = (*vi)->z();
		points->SetPoint(j++, p);
	}
	
	// we need some indexes
	std::unordered_map<Vertex*, int> Vt_list_idx;
	vi = m->VertexBegin(); 
	for (int j = 0; vi!=ve; ++vi) Vt_list_idx[*vi] = j++;
	
	auto Conn = vtkSmartPointer<vtkIdTypeArray>::New();
	Conn->SetNumberOfValues(5*m->TetSize());
	
	auto ti = m->TetBegin(), te = m->TetEnd(); 
	for (int j = 0; ti!=te; ++ti) {
		Conn->SetValue( j++, 4 );
		Conn->SetValue( j++, Vt_list_idx[(*ti)->a()] );
		Conn->SetValue( j++, Vt_list_idx[(*ti)->b()] );
		Conn->SetValue( j++, Vt_list_idx[(*ti)->c()] );
		Conn->SetValue( j++, Vt_list_idx[(*ti)->d()] );
	}

	auto Tets = vtkSmartPointer<vtkCellArray>::New();
	Tets->SetCells(m->TetSize(),Conn);
	
	auto unstructuredGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
	unstructuredGrid->SetPoints(points);
	unstructuredGrid->SetCells(VTK_TETRA, Tets);
	
 
	auto u = vtkSmartPointer<vtkFloatArray>::New();
	u->SetNumberOfComponents(3);
	u->SetName("u");
	u->SetNumberOfValues(m->VertexSize()*3);
	vi = m->VertexBegin(); 
	for (int j = 0; vi!=ve; ++vi) u->SetTuple3(j++, (*vi)->u1(), (*vi)->u2(), (*vi)->u3());
	unstructuredGrid->GetPointData()->AddArray(u);
	
	// Write file
	auto writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
	writer->SetDataModeToAppended();
	writer->EncodeAppendedDataOff();
	std::string fn = fname+"_vol.vtu";
	writer->SetFileName(fn.c_str());
#if VTK_MAJOR_VERSION <= 5
	writer->SetInput(unstructuredGrid);
#else
	writer->SetInputData(unstructuredGrid);
#endif
	writer->Write();
}

#endif
#ifndef HAVE_VTK
void ImplicitPDESystem::Export_Surface(std::string fname) const {
	
	std::ofstream os;
	os.open(fname+"_surf.vtk");
	os << std::setprecision(16);
	os << "# vtk DataFile Version 2.0" << std::endl;
	os << "s3d output data" << std::endl;
	os << "ASCII" << std::endl;
	os << "DATASET POLYDATA" << std::endl;
		
	// we need some indexes
	std::map<Vertex*, int> Vt_list_idx;
	Mesh::const_VertexIt vi = m->BdryVertexBegin(), ve = m->BdryVertexEnd();
	int idx = 0;
	for (; vi!=ve; ++vi) {
		Vt_list_idx[*vi] = idx;
		++idx;
	}
	
	os << "POINTS " << m->BdryVertexSize() << " float" << std::endl;
	vi = m->BdryVertexBegin(); 
	for (; vi!=ve; ++vi)
		os << (*vi)->x() << " " << (*vi)->y() << " " << (*vi)->z() << std::endl;

	os << "POLYGONS " << m->BdryFaceSize() << " " << m->BdryFaceSize()*4 << std::endl;
	Mesh::const_FaceIt fi = m->BdryFaceBegin(), fe = m->BdryFaceEnd();
	for (; fi!=fe; ++fi)
		os << "3 " << Vt_list_idx[(*fi)->a()] << " " << Vt_list_idx[(*fi)->b()] << " "
			<< Vt_list_idx[(*fi)->c()] << std::endl;
    
	os << "POINT_DATA " << m->BdryVertexSize() << std::endl;
	os << "SCALARS u float 3" << std::endl;
	os << "LOOKUP_TABLE default" << std::endl;
	vi = m->BdryVertexBegin();
	for (; vi!=ve; ++vi)
		os << (*vi)->u1() << " " << (*vi)->u2() << " " << (*vi)->u3() << std::endl;
	
	os << "SCALARS dir_bc float 3" << std::endl;
	os << "LOOKUP_TABLE default" << std::endl;
	vi = m->BdryVertexBegin();
	for (; vi!=ve; ++vi) {
		for (int j=0; j<3; ++j) {
			os << ( (*vi)->DirBdry(j) ? 1.0 : 0.0 ) << " ";
		}
		os << std::endl;
	}
	os << "SCALARS neu_bc float 3" << std::endl;
	os << "LOOKUP_TABLE default" << std::endl;
	vi = m->BdryVertexBegin();
	for (; vi!=ve; ++vi) {
		Vec3 coord = (*vi)->Coord();
		for (int j=0; j<3; ++j) {
			os << NeumannBC(coord)[j] << " ";
		}
		os << std::endl;
	}
	
	/*os << "CELL_DATA " << m->TetSize() << std::endl;
	os << "SCALARS Tet_Volume float" << std::endl;
	os << "LOOKUP_TABLE default" << std::endl;
	ti = m->TetBegin();
	for (; ti!=te; ++ti)
	os << (*ti)->Vol() << std::endl;
	
	os << "SCALARS grad_u float 9" << std::endl;
	os << "LOOKUP_TABLE default" << std::endl;
	ti = m->TetBegin();
	for (; ti!=te; ++ti) {
	Mat33 grad_u;
	(*ti)->Calc_grad_u(grad_u);
	os << grad_u.row1().x() << " " << grad_u.row1().y() << " " << grad_u.row1().z() << " ";
	os << grad_u.row2().x() << " " << grad_u.row2().y() << " " << grad_u.row2().z() << " ";
	os << grad_u.row3().x() << " " << grad_u.row3().y() << " " << grad_u.row3().z() << std::endl;
	} */
	os.close();
	
    
}
void ImplicitPDESystem::Export_Vol(std::string fname) const {
	
	
	std::ofstream os;
	os.open(fname+"_vol.vtk");
	os << std::setprecision(16);
	os << "# vtk DataFile Version 2.0" << std::endl;
	os << "s3d output data" << std::endl;
	os << "ASCII" << std::endl;
	os << "DATASET UNSTRUCTURED_GRID" << std::endl;
	
	os << "POINTS " << m->VertexSize() << " float" << std::endl;
	Mesh::const_VertexIt vi = m->VertexBegin(), ve = m->VertexEnd();
	for (; vi!=ve; ++vi)
		os << (*vi)->x() << " " << (*vi)->y() << " " << (*vi)->z() << std::endl;

	// we need some indexes
	std::map<Vertex*, int> Vt_list_idx;
	vi = m->VertexBegin(); int idx = 0;
	for (; vi!=ve; ++vi) {
		Vt_list_idx[*vi] = idx;
		++idx;
	}

	os << "CELLS " << m->TetSize() << " " << m->TetSize()*5 << std::endl;
	Mesh::const_TetIt ti = m->TetBegin(), te = m->TetEnd();
	for (; ti!=te; ++ti)
		os << "4 " << Vt_list_idx[(*ti)->a()] << " " << Vt_list_idx[(*ti)->b()] << " "
			<< Vt_list_idx[(*ti)->c()] << " " << Vt_list_idx[(*ti)->d()] << std::endl;
    
	os << "CELL_TYPES " << m->TetSize() << std::endl;
	ti = m->TetBegin();
	for (; ti!=te; ++ti)
		os << "10" << std::endl;
    
	os << "POINT_DATA " << m->VertexSize() << std::endl;
	os << "SCALARS u float 3" << std::endl;
	os << "LOOKUP_TABLE default" << std::endl;
	vi = m->VertexBegin();
	for (; vi!=ve; ++vi)
		os << (*vi)->u1() << " " << (*vi)->u2() << " " << (*vi)->u3() << std::endl;
	
	//os << "SCALARS bc float 3" << std::endl;
	//os << "LOOKUP_TABLE default" << std::endl;
	//vi = m->VertexBegin();
	//for (; vi!=ve; ++vi) {
	//	Vec3 coord = (*vi)->Coord();
	//	for (int j=0; j<3; ++j) {
	//		os << ( (*vi)->DirBdry(j) ? -1000.0 : ( (*vi)->Boundary() ? NeumannBC(coord)[j] : -2000.0) ) << " ";
	//	}
	//	os << std::endl;
	//}
	
	os << "CELL_DATA " << m->TetSize() << std::endl;
	os << "SCALARS Tet_Volume float" << std::endl;
	os << "LOOKUP_TABLE default" << std::endl;
	ti = m->TetBegin();
	for (; ti!=te; ++ti)
		os << (*ti)->Vol() << std::endl;
	
	os << "SCALARS grad_u float 9" << std::endl;
	os << "LOOKUP_TABLE default" << std::endl;
	ti = m->TetBegin();
	for (; ti!=te; ++ti) {
		Mat33 grad_u;
		(*ti)->Calc_grad_u(grad_u);
		os << grad_u.row1().x() << " " << grad_u.row1().y() << " " << grad_u.row1().z() << " ";
		os << grad_u.row2().x() << " " << grad_u.row2().y() << " " << grad_u.row2().z() << " ";
		os << grad_u.row3().x() << " " << grad_u.row3().y() << " " << grad_u.row3().z() << std::endl;
	}
	os << "SCALARS strain float 6" << std::endl;
	os << "LOOKUP_TABLE default" << std::endl;
	ti = m->TetBegin();
	for (; ti!=te; ++ti) {
		Mat33 grad_u;
		(*ti)->Calc_grad_u(grad_u);
		os << grad_u.row1().x() << " " << grad_u.row2().y() << " " << grad_u.row3().z() << " ";
		os << 0.5*(grad_u.row3().y()+grad_u.row2().z()) << " ";
		os << 0.5*(grad_u.row3().x()+grad_u.row1().z()) << " ";
		os << 0.5*(grad_u.row2().x()+grad_u.row1().y()) << std::endl;
	} 
		
	os << "SCALARS E float 1" << std::endl;
	os << "LOOKUP_TABLE default" << std::endl;
	ti = m->TetBegin();
	for (; ti!=te; ++ti) {
		Mat33 grad_u;
		(*ti)->Calc_grad_u(grad_u);
		double E = 0.0;
		for (int l=0; l<3; ++l)
			for (int k=0; k<3; ++k)
				for (int j=0; j<3; ++j)
					for (int i=0; i<3; ++i)
						E += grad_u.row(i)[j]*C(i,j,k,l)*grad_u.row(l)[k];
		os << E << std::endl;
	} 
		
		
		os << "SCALARS P float 1" << std::endl;
	os << "LOOKUP_TABLE default" << std::endl;
	ti = m->TetBegin();
	for (; ti!=te; ++ti) {
		//std::cout<< "Hallo"<< std::endl;
		Mat33 grad_u;
		(*ti)->Calc_grad_u(grad_u);
		 Eigen::MatrixXd e_strain(3, 3);
		 Eigen::MatrixXd stress(3, 3);
		//Mat33 e_strain,stress;

		e_strain(0,0)= grad_u.row1().x();
		e_strain(1,1)= grad_u.row2().y();
		e_strain(2,2)= grad_u.row3().z();

		e_strain(2,1)= 0.5*(grad_u.row3().y()+grad_u.row2().z());
		e_strain(1,2)=e_strain(2,1);

		e_strain(2,0)= 0.5*(grad_u.row3().x()+grad_u.row1().z());
		e_strain(0,2)=e_strain(2,0);

		e_strain(1,0)= 0.5*(grad_u.row2().x()+grad_u.row1().y());
		e_strain(0,1)=e_strain(1,0);

	for (int i=0; i<3; ++i){
		for (int j=0; j<3; ++j){
			for (int k=0; k<3; ++k){
				for (int l=0; l<3; ++l){
				
					stress(i,j)+= C(i,j,k,l)*e_strain(k,l);

			}
 		}

	}
}

}
		
		
		
		
	os.close();
	
	
}


#endif


void ImplicitPDESystem::SaveFace(std::string fname) const {
	
	
	std::ofstream os;
	os.open(fname);
	os << std::setprecision(16);
	os << "#obj file of surface" << std::endl;
	
	Mesh::const_VertexIt vi = m->BdryVertexBegin(), ve = m->BdryVertexEnd();
	for (; vi!=ve; ++vi)
		os << "v " << (*vi)->x() << " " << (*vi)->y() << " " << (*vi)->z() << std::endl;
	
	// we need some indexes
	std::map<Vertex*, int> Vt_list_idx;
	vi = m->BdryVertexBegin(); int idx = 0;
	for (; vi!=ve; ++vi) {
		Vt_list_idx[*vi] = idx;
		++idx;
	}
	
	Mesh::const_FaceIt fi = m->BdryFaceBegin(), fe = m->BdryFaceEnd();
	for (; fi!=fe; ++fi)
		os << "f " << Vt_list_idx[(*fi)->a()]+1 << " " << Vt_list_idx[(*fi)->b()]+1 << " "
			<< Vt_list_idx[(*fi)->c()]+1 << std::endl;
    
	os.close();
	
}
