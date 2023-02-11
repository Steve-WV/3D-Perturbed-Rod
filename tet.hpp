#ifndef __TET_HPP__
#define __TET_HPP__

#include <iostream> 
#include <algorithm> 
#include <vector>

#include "names.h"
#include "utils.hpp"
#include "tmv.hpp"
using namespace std; 


class Vertex;

/// Topological Triangles.
class Tet{
public:

  /// Default constructor
  inline Tet() { 
	/*m_cb[A] = m_cb[B]= m_cb[C] = m_cb[D] = false;*/
    m_ts[A] = m_ts[B]= m_ts[C] = m_ts[D] = NULL; }
  /// Constructor from three vertices.
  inline Tet( Vertex*, Vertex*, Vertex*, Vertex* );
  /// Copy constructor.
  inline Tet( const Tet& t );
  /// Destructor removes the triangle connections.
  inline ~Tet();

  /// Query routines
  inline bool         Contains( const Vertex* ) const;
  inline bool         IsAdj( const Tet* ) const;
  inline VertexName   VName( const Vertex* ) const;
  inline VertexName   OtherVName( const Vertex*, const Vertex*, const Vertex* ) const;
  inline VertexName   VNameAcross( const Tet* ) const;
  
  /// Adjacency information
  inline Tet*&   Across( VertexName );
  inline Tet*    Across( VertexName ) const;
  inline Tet*&   Across( const Vertex*, const Vertex*, const Vertex* );
  inline Tet*    Across( const Vertex*, const Vertex*, const Vertex* ) const;
  inline Tet*&   Across( const Vertex* );
  inline Tet*    Across( const Vertex* ) const;
  
  inline Vertex*&     OtherV( const Vertex*, const Vertex*, const Vertex* );
  inline Vertex*      OtherV( const Vertex*, const Vertex*, const Vertex* ) const;
  inline Vertex*&     Across( const Tet* );
  inline Vertex*      Across( const Tet* ) const;
  
  /// Connect Bits related stuff
  inline void         Unlink( const Tet* );
  inline void         SetConnected( const Vertex * );
  inline bool         AllVerticesConnected() const;

  /// Accessors
  inline Vertex*             a( void ) const { return m_vs[A]; }
  inline Vertex*             b( void ) const { return m_vs[B]; }
  inline Vertex*             c( void ) const { return m_vs[C]; }
  inline Vertex*             d( void ) const { return m_vs[D]; }
  inline Vertex*             v( const VertexName vn ) const { assert( vn != VNONE ); return m_vs[vn]; }

  inline int Valence( const VertexName ) const;
  inline int Valence( const Vertex* ) const;
  //inline bool BoundaryAcross( const VertexName ) const;
  //inline bool BoundaryAcross( const Vertex* ) const;
  //inline bool Boundary() const;
    
    inline double Vol() const { return m_volume; }
    inline double Init_vol();
		inline double Init_mesh_size();
	
	/*inline void Calc_Ftinv();
	inline void Calc_grad_u_loc(Vec3&) const;
	inline void Calc_grad_v_loc(Vec3&) const;
	*/
    
	inline void Calc_grad_u1(Vec3&) const;
	inline void Calc_grad_u2(Vec3&) const;
	inline void Calc_grad_u3(Vec3&) const;
	inline void Calc_grad_u(int, Vec3&) const;
	inline void	Calc_grad_u(Mat33& grad_u) const;
	
	inline void Calc_u1(double*) const;
	inline void Calc_u2(double*) const;
	inline void Calc_u3(double*) const;
	inline void Calc_u_scal(double*) const;
	inline void Calc_u(Vec3*) const;
	inline void Calc_u(int, double*) const;
	
	
	inline void Calc_s(int, double*) const;
	inline void Init_grad_s();
	inline void Calc_grad_s( int, Vec3&) const;

private:
  /// Vertices of the tet	
  Vertex*		m_vs[4];
  /// Triangle across from vs[i]	
  Tet*			m_ts[4];
    
    
    // some properties of tets
    double m_volume;
	Mat33 Ftinv;
    Vec3 grad_s[4];
	
  // is there a face connected across that vertex?
  /*bool			m_cb[4];*/
};

//#include "vertex.hpp"

inline Tet::Tet( Vertex* a, Vertex* b, Vertex* c, Vertex* d )
{
  assert( a && b && c && d );
  assert( ( a != b ) && ( b != c ) && ( c != a ) && ( d != a ) && ( d != b ) && ( d != c ) );
  // set the vertices
  m_vs[A] = a;  m_vs[B] = b;  m_vs[C] = c; m_vs[D] = d;
  // set all the connect bits to false
  /*m_cb[A] = m_cb[B] = m_cb[C] = m_cb[D] = false;*/
  // set the adjacencies to null for now
  m_ts[A] = m_ts[B] = m_ts[C] = m_ts[D] = NULL;
}

inline Tet::Tet( const Tet& t ) {
  m_vs[A] = t.m_vs[A]; m_vs[B] = t.m_vs[B]; m_vs[C] = t.m_vs[C]; m_vs[D] = t.m_vs[D];
  m_ts[A] = t.m_ts[A]; m_ts[B] = t.m_ts[B]; m_ts[C] = t.m_ts[C]; m_ts[D] = t.m_ts[D];
  /*m_cb[A] = t.m_cb[A]; m_cb[B] = t.m_cb[B]; m_cb[C] = t.m_cb[C]; m_cb[D] = t.m_cb[D];*/
  
}

inline Tet::~Tet( void ) {
  // unlink from neighbors
  if( m_ts[A] ){ m_ts[A]->Unlink( this ); m_ts[A] = NULL; }
  if( m_ts[B] ){ m_ts[B]->Unlink( this ); m_ts[B] = NULL; }
  if( m_ts[C] ){ m_ts[C]->Unlink( this ); m_ts[C] = NULL; }
  if( m_ts[D] ){ m_ts[D]->Unlink( this ); m_ts[D] = NULL; }
}

/// is this triangle surrounded by triangles?
/*inline bool
Tet::Boundary() const {
  return  m_vs[A]->Boundary() || m_vs[B]->Boundary() || m_vs[C]->Boundary() || m_vs[D]->Boundary() ;
}*/

/// Is the given vertex one of the three that make up the face
inline bool
Tet::Contains( const Vertex *v ) const
{
  return ( m_vs[A] == v ) || ( m_vs[B] == v ) || ( m_vs[C] == v ) || ( m_vs[D] == v );
}

inline bool
Tet::IsAdj( const Tet *t ) const
{
  return ( m_ts[A] == t ) || ( m_ts[B] == t ) || ( m_ts[C] == t ) || ( m_ts[D] == t );
}

inline VertexName
Tet::VName( const Vertex *v) const
{	
	assert( Contains(v) );
    if ( m_vs[A] == v ) return A;
    if ( m_vs[B] == v ) return B;
    if ( m_vs[C] == v ) return C;
    return D; //else
}

/// Given two of the vertices in the face, return the name of the third vertex
inline VertexName     
Tet::OtherVName( const Vertex* v1, const Vertex* v2, const Vertex* v3 ) const {
  assert( v1 && v2 && v3 ); assert( v1 != v2 && v2 != v3 && v3 != v1 );
  assert( Contains( v1 ) && Contains( v2 ) && Contains( v3 ) );
	
  VertexName a = VName( v1 );
  VertexName b = VName( v2 );
  VertexName c = VName( v3 );
  
  return OtherVN(a,b,c); 
}

inline VertexName
Tet::VNameAcross( const Tet* t) const {
  assert( IsAdj(t) );
  
  if ( m_ts[A] == t ) return A;
  if ( m_ts[B] == t ) return B;
  if ( m_ts[C] == t ) return C;
  return D; //else
}

inline Tet*&
Tet::Across( VertexName i ) {
  assert( i<3 );
  return m_ts[i];
}

inline Tet*
Tet::Across( VertexName i ) const {
  assert( i<3 );
  return m_ts[i];
}

/// Return the tet that is across the face defined by the three
/// Vertices
inline Tet*&
Tet::Across( const Vertex *v1, const Vertex *v2, const Vertex *v3  ) {
	return m_ts[OtherVName( v1, v2, v3 )];
}

/// Return the face that is across the edge defined by the three
/// Vertices. This one is const so that others that don't need to 
/// change Tet can use it.
inline Tet*
Tet::Across( const Vertex* v1, const Vertex* v2, const Vertex* v3 ) const {
  return m_ts[OtherVName( v1, v2, v3 )];
}

inline Tet*&
Tet::Across( const Vertex* v ) {
  return m_ts[VName( v )];
}

inline Tet*
Tet::Across( const Vertex* v ) const {
  return m_ts[VName( v )];
}

inline Vertex*&
Tet::OtherV( const Vertex* v1, const Vertex* v2, const Vertex* v3 ) {
	return m_vs[OtherVName( v1, v2, v3 )];
}

inline Vertex*
Tet::OtherV( const Vertex* v1, const Vertex* v2, const Vertex* v3 ) const {
	return m_vs[OtherVName( v1, v2, v3 )];
}

inline Vertex*&
Tet::Across( const Tet* t ) {
	return m_vs[VNameAcross( t )];
}
inline Vertex*
Tet::Across( const Tet* t ) const {
	return m_vs[VNameAcross( t )];
}

/// Remove the connections for this triangle.
inline void
Tet::Unlink( const Tet* t ){
  if( m_ts[A] == t ) m_ts[A] = NULL;
  else if( m_ts[B] == t ) m_ts[B] = NULL;
  else if( m_ts[C] == t ) m_ts[C] = NULL;
  else if( m_ts[D] == t ) m_ts[D] = NULL;
  else die();
}


#include "vertex.hpp"
inline double
Tet::Init_vol() {
   
    Vertex *a = m_vs[0], *b = m_vs[1], *c = m_vs[2], *d = m_vs[3];
    
    Vec3 DA = d->Coord() - a->Coord();
    Vec3 DB = d->Coord() - b->Coord();
    Vec3 DC = d->Coord() - c->Coord();
    
    m_volume = 1.0/6.0 * dot( DA, cross(DB, DC));

    return m_volume;
}

inline double
Tet::Init_mesh_size() {
   
    Vertex *a = m_vs[0], *b = m_vs[1], *c = m_vs[2], *d = m_vs[3];
    
    Vec3 DA = d->Coord() - a->Coord();
    Vec3 DB = d->Coord() - b->Coord();
    Vec3 DC = d->Coord() - c->Coord();
    
		std::vector<double> mesh_vec; 
    mesh_vec.push_back(DA.norm()),mesh_vec.push_back(DB.norm()),mesh_vec.push_back(DC.norm());
		double mesh_size = *std::max_element(mesh_vec.begin(),mesh_vec.end());
    return mesh_size;
}

inline void
Tet::Calc_grad_u1(Vec3& grad_u) const {

    grad_u = a()->u1() * grad_s[0];
    grad_u += b()->u1() * grad_s[1];
    grad_u += c()->u1() * grad_s[2];
    grad_u += d()->u1() * grad_s[3];
}
inline void
Tet::Calc_grad_u2(Vec3& grad_u) const {

    grad_u = a()->u2() * grad_s[0];
    grad_u += b()->u2() * grad_s[1];
    grad_u += c()->u2() * grad_s[2];
    grad_u += d()->u2() * grad_s[3];
}
inline void
Tet::Calc_grad_u3(Vec3& grad_u) const {

    grad_u = a()->u3() * grad_s[0];
    grad_u += b()->u3() * grad_s[1];
    grad_u += c()->u3() * grad_s[2];
    grad_u += d()->u3() * grad_s[3];
}
inline void
Tet::Calc_grad_u(int j, Vec3& grad_u) const {

    grad_u = a()->u(j) * grad_s[0];
    grad_u += b()->u(j) * grad_s[1];
    grad_u += c()->u(j) * grad_s[2];
    grad_u += d()->u(j) * grad_s[3];
}
inline void
Tet::Calc_grad_u(Mat33& grad_u) const {

    Vec3 grad_u1, grad_u2, grad_u3;
	Calc_grad_u1(grad_u1);
	Calc_grad_u2(grad_u2);
	Calc_grad_u3(grad_u3);
	grad_u = Mat33(grad_u1, grad_u2, grad_u3).t();
	
}

inline void
Tet::Calc_u1(double* u) const {
	for (int k = 0; k<NumIntPts; ++k) {
		
		u[k] = GaussPoints[k][0] * a()->u1() + 
			GaussPoints[k][1] * b()->u1() +
				GaussPoints[k][2] * c()->u1() +
					GaussPoints[k][3] * d()->u1();
	} 	
}
inline void
Tet::Calc_u2(double* u) const {
	for (int k = 0; k<NumIntPts; ++k) {
		
		u[k] = GaussPoints[k][0] * a()->u2() + 
			GaussPoints[k][1] * b()->u2() +
				GaussPoints[k][2] * c()->u2() +
					GaussPoints[k][3] * d()->u2();
	} 	
}
inline void
Tet::Calc_u3(double* u) const {
	for (int k = 0; k<NumIntPts; ++k) {
		
		u[k] = GaussPoints[k][0] * a()->u3() + 
			GaussPoints[k][1] * b()->u3() +
				GaussPoints[k][2] * c()->u3() +
					GaussPoints[k][3] * d()->u3();
	} 	
}

inline void
Tet::Calc_u_scal(double* u) const {
	for (int k = 0; k<NumIntPts; ++k) {
		
		u[k] = GaussPoints[k][0] * a()->u_scal() + 
			GaussPoints[k][1] * b()->u_scal() +
				GaussPoints[k][2] * c()->u_scal() +
					GaussPoints[k][3] * d()->u_scal();
	} 	
}

inline void
Tet::Calc_u(int j, double* u) const {
	for (int k = 0; k<NumIntPts; ++k) {
		
		u[k] = GaussPoints[k][0] * a()->u(j) + 
			GaussPoints[k][1] * b()->u(j) +
				GaussPoints[k][2] * c()->u(j) +
					GaussPoints[k][3] * d()->u(j);
	} 	
}
inline void
Tet::Calc_u(Vec3* u) const {
	for (int k = 0; k<NumIntPts; ++k) {
		
		u[k] = GaussPoints[k][0] * a()->u() + 
			GaussPoints[k][1] * b()->u() +
				GaussPoints[k][2] * c()->u() +
					GaussPoints[k][3] * d()->u();
	} 	
}

inline void 
Tet::Calc_s(int j, double* s) const {
    
    for (int k= 0; k<NumIntPts; ++k){
        s[k] = GaussPoints[k][j];
    }
	
}

inline void
Tet::Init_grad_s() {
	
    Mat33 F (a()->Coord() - d()->Coord(), b()->Coord() - d()->Coord(), c()->Coord() - d()->Coord() );
    F = F.inv();
    grad_s[0] = F.row1();
    grad_s[1] = F.row2();
    grad_s[2] = F.row3();
    grad_s[3] = -1.0* ( grad_s[0] + grad_s[1] + grad_s[2] );
	
}

inline void
Tet::Calc_grad_s(int j, Vec3& v) const {
	v = grad_s[j];
}

#endif	/* __TET_HPP__ */
