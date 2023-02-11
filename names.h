#ifndef	__NAMES_H__
#define	__NAMES_H__

enum VertexName { A = 0, B = 1, C = 2, D = 3, VNONE = 4 };

// turn one type of number (rightmost letter) into another kind of number
// (left letter)
inline VertexName   OtherVN( const VertexName v1, const VertexName v2, const VertexName v3 ){
	assert( v1 != v2 && v1 != v3 && v2 != v3 ); return VertexName(6 - v1 - v2 - v3);
}

inline VertexName   NextVN( const VertexName v ) { return VertexName((v+1)%4); } 
inline VertexName   PrevVN( const VertexName v ) { return VertexName((v+3)%4); } 

#endif	/* __NAMES_H__ */
