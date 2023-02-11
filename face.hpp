#ifndef __FACE_HPP__
#define __FACE_HPP__

#include <functional>
#include "utils.hpp"

class Tet;
class Vertex;

/// Faces in the triangulation.
class Face{
public:
	/// Construction from two vertices, and maybe an adjacent triangle.
	inline Face( Vertex*, Vertex*, Vertex*, Tet* = NULL );
	/// Copy constructor. Would prefer to do without..
	inline Face( const Face& );

	/// Topological information
	inline bool         Boundary() const;

	// Information about the vertices
	inline Vertex*      OtherV( const Vertex*, const Vertex* ) const;
	inline bool         Contains( const Vertex* ) const;
	inline bool         IsFace( const Vertex*, const Vertex*, const Vertex* ) const;

	// Accessors
	inline Vertex*&            a() { return m_vi[0]; }
	inline Vertex*             a() const { return m_vi[0]; }
	inline Vertex*&            b() { return m_vi[1]; }
	inline Vertex*             b() const { return m_vi[1]; }
	inline Vertex*&            c() { return m_vi[2]; }
	inline Vertex*             c() const { return m_vi[2]; }
	inline Tet*                t() const { return m_t; }
	inline Tet*&               t() { return m_t; }
	
	inline double Area() const;
	
	inline void MarkNeuBdry1() { m_neubd[0] = true; }
	inline void MarkNeuBdry2() { m_neubd[1] = true; }
	inline void MarkNeuBdry3() { m_neubd[2] = true; }
	inline void MarkNeuBdry(int j) { m_neubd[j] = true; }
	
	inline bool NeuBdry1() const { return m_neubd[0]; }
	inline bool NeuBdry2() const { return m_neubd[1]; }
	inline bool NeuBdry3() const { return m_neubd[2]; }
	inline bool NeuBdry(int j) const { return m_neubd[j]; }
	
	
	// Comparison function for use in Face sets, lexicographical order
	struct face_comp : public std::binary_function<Face, Face, bool> {
		face_comp() {}
		bool operator()( const Face f1, const Face f2 ) const {
			const Vertex *min1 = min( f1.a(), f1.b(), f1.c() );
			const Vertex *min2 = min( f2.a(), f2.b(), f2.c() );
			const Vertex *max1 = max( f1.a(), f1.b(), f1.c() );
			const Vertex *max2 = max( f2.a(), f2.b(), f2.c() );
			const Vertex *mid1 = f1.OtherV( min1, max1 );
            const Vertex *mid2 = f2.OtherV( min2, max2 );
            
            assert( mid1 != min1 && mid1 != max1 && min1 != max1 );
            assert( mid2 != min2 && mid2 != max2 && min2 != max2 );
            
            if ( min1 <  min2 ) return true;
			if ( min1 == min2 && mid1 <  mid2 ) return true;
			if ( min1 == min2 && mid1 == mid2 && max1 < max2 ) return true;
			return false;
		}
	};
	
private:
	
	
	//@{
	/// Vertices defining end point	
	Vertex*       m_vi[3];
	/// one of the tets incident with this Face	
	Tet*          m_t;
	
	bool m_neubd[3];
	
};

#include "tet.hpp"

inline
	Face::Face( Vertex* a, Vertex* b, Vertex* c, Tet* t )
		: m_t( t )
{
	assert( a && b && c ); assert( a != b && b != c && c != a );
	assert( !t || ( t->Contains(a) && t->Contains(b) && t->Contains(c) ) );
  
	m_vi[0] = a;
	m_vi[1] = b;
	m_vi[2] = c;
}

inline
	Face::Face( const Face& e )
	 	: m_t( e.m_t )
{
	assert( e.m_vi[0] && e.m_vi[1] && e.m_vi[2] );
	assert( e.m_vi[0] != e.m_vi[1] && e.m_vi[1] != e.m_vi[2] && e.m_vi[2] != e.m_vi[0] );
	assert( !e.m_t || ( e.m_t->Contains(e.m_vi[0]) && e.m_t->Contains(e.m_vi[1]) && e.m_t->Contains(e.m_vi[2]) ) );
  
	m_vi[0] = e.m_vi[0];
	m_vi[1] = e.m_vi[1];
	m_vi[2] = e.m_vi[2];
}


/// Is the Face on a boundary?
inline bool
Face::Boundary() const{
	// if the adjacent face does not have a neighbor pointer than this is a 
	// boundary.
	return !m_t->Across( m_vi[0], m_vi[1], m_vi[2] );
}

/// Given one vertex, return the other one.
inline Vertex*
Face::OtherV( const Vertex* v1, const Vertex* v2 ) const{
	assert( v1 && v2 ); assert( Contains( v1 ) && Contains( v2 ) );
    assert( v1 != v2 );
	if ( m_vi[0] == v1 && m_vi[1] == v2 ) return m_vi[2];
	if ( m_vi[1] == v1 && m_vi[0] == v2 ) return m_vi[2];
	if ( m_vi[0] == v1 && m_vi[2] == v2 ) return m_vi[1];
	if ( m_vi[2] == v1 && m_vi[0] == v2 ) return m_vi[1];
    assert( (m_vi[1] == v1 && m_vi[2] == v2) || (m_vi[2] == v1 && m_vi[1] == v2) );
	return m_vi[0];
}

/// Is the given vertex one of this faces's points.
inline bool
Face::Contains( const Vertex* v ) const{
	assert( v );
	return v == m_vi[0] || v == m_vi[1] || v == m_vi[2];
}

/// Does this Face contain the three vertices
inline bool
	Face::IsFace( const Vertex* v1, const Vertex* v2, const Vertex* v3 ) const
{
	assert( v1 && v2 && v3 );
	if ( m_vi[0] == v1 && m_vi[1] == v2 && m_vi[2] == v3 ) return true;
	if ( m_vi[1] == v1 && m_vi[2] == v2 && m_vi[0] == v3 ) return true;
	if ( m_vi[2] == v1 && m_vi[0] == v2 && m_vi[1] == v3 ) return true;
	if ( m_vi[1] == v1 && m_vi[0] == v2 && m_vi[2] == v3 ) return true;
	if ( m_vi[2] == v1 && m_vi[1] == v2 && m_vi[0] == v3 ) return true;
	if ( m_vi[0] == v1 && m_vi[2] == v2 && m_vi[1] == v3 ) return true;
	return false;
}

#include "vertex.hpp"
inline double 
	Face::Area() const {
		
		Vec3 a = (m_vi[0])->Coord();
		Vec3 b = (m_vi[1])->Coord();
		Vec3 c = (m_vi[2])->Coord();
		
		return 0.5*(cross(b-a,c-a)).norm();
		
	}

#endif	/* __Face_HPP__ */
