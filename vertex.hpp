#ifndef	__VERTEX_HPP__
#define	__VERTEX_HPP__

// include files
#include "names.h"
#include "utils.hpp"
//#include "tet.hpp"
#include "tmv.hpp"

class Mesh;
/// The vertices. They are pretty dumb.
class Vertex {
public:
  
	inline Vertex() : m_boundary(false) { 
		m_dirbd[0] = false;
		m_dirbd[1] = false;
		m_dirbd[2] = false; 
	}
  
	inline Vertex( Vec3 coord ) : m_coord( coord ) { 
		m_dirbd[0] = false;
		m_dirbd[1] = false;
		m_dirbd[2] = false; 
	}
  
	inline Vertex( double x, double y, double z ) : m_boundary(false) { 
		m_dirbd[0] = false;
		m_dirbd[1] = false;
		m_dirbd[2] = false; 
		m_coord.x() = x; 
		m_coord.y() = y; 
		m_coord.z() = z; 
	}

	inline ~Vertex() {}

	// Accessors
	inline Vec3& Coord() { return m_coord; }
	inline Vec3 Coord() const { return m_coord; }
	inline double& x() { return m_coord.x(); }
	inline double x() const { return m_coord.x(); }
	inline double& x1() { return m_coord.x(); }
	inline double x1() const { return m_coord.x(); }
	inline double& y() { return m_coord.y(); }
	inline double y() const { return m_coord.y(); }
	inline double& z() { return m_coord.z(); }
	inline double z() const { return m_coord.z(); }
  
	inline int Index_ref() const { return m_index_ref; }
  inline void SetIndex_ref(int i) { m_index_ref = i; }
	inline int Index1() const { return m_index[0]; }
	inline void SetIndex1(int i) { m_index[0] = i; }
	inline int Index2() const { return m_index[1]; }
	inline void SetIndex2(int i) { m_index[1] = i; }
	inline int Index3() const { return m_index[2]; }
	inline void SetIndex3(int i) { m_index[2] = i; }
	inline int Index(int j) const { return m_index[j]; }
	inline void SetIndex(int j, int i) { m_index[j] = i; }
	
  
	//inline Tet*& Tet() { return m_ti; }
	//inline Tet* Tet() const { return m_ti; }
	//inline int Valence() const { return m_valence; }

	// Is this a boundary vertex?
	inline bool Boundary() const { return m_boundary; }
	inline bool DirBdry1() const { return m_dirbd[0]; }
	inline bool DirBdry2() const { return m_dirbd[1]; }
	inline bool DirBdry3() const { return m_dirbd[2]; }
	inline bool DirBdry(int j) const { return m_dirbd[j]; }
  
	inline void MarkDirBdry1() { m_dirbd[0] = true; }
	inline void MarkDirBdry2() { m_dirbd[1] = true; }
	inline void MarkDirBdry3() { m_dirbd[2] = true; }
	inline void MarkDirBdry(int j) { m_dirbd[j] = true; }
  
	inline double& u1() {return *(m_u[0]);}
	inline double u1() const {return *(m_u[0]);}
	inline double& u2() {return *(m_u[1]);}
	inline double u2() const {return *(m_u[1]);}
	inline double& u3() {return *(m_u[2]);}
	inline double u3() const {return *(m_u[2]);}
	inline double& u(int j) {return *(m_u[j]);}
	inline double u(int j) const {return *(m_u[j]);}
	inline double& u_scal() {return *(m_u_scal);}
	inline double u_scal() const {return *(m_u_scal);}
  
	inline Vec3 u() const { return Vec3(*(m_u[0]), *(m_u[1]), *(m_u[2])); }
	
	/*inline double& f() {return m_f;}
	inline double f() const {return m_f;}
	inline double& mass(){return m_mass;}
	inline double mass() const {return m_mass;}*/

	//inline void AttachX1( double u_add) { x_val=u_add; }
	//inline void AttachX2( double u_add ) { m_coord.y()=u_add; }
	//inline void AttachX3( double u_add ) { m_coord.z()=u_add;  }
  
	inline void Attach_u1( double* u_addr ) { m_u[0] = u_addr; }
	inline void Attach_u2( double* u_addr ) { m_u[1] = u_addr; }
	inline void Attach_u3( double* u_addr ) { m_u[2] = u_addr; }
	inline void Attach_u( int j, double* u_addr ) { m_u[j] = u_addr; }
	inline void Attach_u_scal( double* u_addr ) { m_u_scal = u_addr; }
	inline void Attach_u( double* u_addr1, double* u_addr2, double* u_addr3 ) { 
		m_u[0] = u_addr1;
		m_u[1] = u_addr2;
		m_u[2] = u_addr3;
	}
  
	// Test if this vertex is on the boundary. 
	// The vertex is on a boundary if there is no clockwise triangle,
	// because we set up the triangle adjacency to the most clockwise triangle
	// if we were on a boundary, or to any triangle if we could march around.
	// not sure what to do here yet...
	//inline void MarkBoundary( void ) { _boundary = !_ti->CWTriangle( this ); }

	// also not sure yet.
	//void SetValence();
	Vec3 m_coord;	  

	friend class Mesh;
  
private:
	
	int m_index[3];
	int m_index_ref;
	
	
	// One of the incident tets.
	// Clockwise most incident triangle if the one ring is open. Any face if
	// it's closed.
	// not sure yet
	//Tet* m_ti;       

	/// valence of this vertex;
	//int m_valence;

	bool m_boundary;
	bool m_dirbd[3];
	inline void MarkBoundary() { m_boundary = true; }
	/// is this vertex on the boundary?
  
	double *m_xi;
	double *m_u[3];
	double *m_u_scal;
};

#endif	/* __VERTEX_HPP__ */
