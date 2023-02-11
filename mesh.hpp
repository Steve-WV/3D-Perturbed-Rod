#ifndef __MESH_HPP__
#define __MESH_HPP__

#include <vector>
#include <set>

#include "face.hpp"
#include "utils.hpp"

class Tet;
class Vertex;

// The mesh.
class Mesh  {
	public:
  
  //@{
  /// Typedefs for vertex and triangle containers in the mesh.
  typedef             std::vector<Vertex*> VertexCt;
  typedef             std::vector<Tet*> TetCt;
  typedef			  std::vector<Face*> FaceCt;
  	
  typedef             TetCt::iterator TetIt;
  typedef             TetCt::const_iterator const_TetIt;
  typedef             VertexCt::iterator   VertexIt;
  typedef             VertexCt::const_iterator   const_VertexIt;
  typedef             FaceCt::iterator FaceIt;
  typedef			  FaceCt::const_iterator const_FaceIt;
 
  
  typedef             std::set<Face, Face::face_comp>  FaceSet;
  typedef             std::set<Vertex*>  VertexSet;
  	
  /// Constructor initializes empty mesh.
  inline Mesh( void ) : m_vc(), m_tc() {}
  /// Destructor cleans up vertex and triangle containers.
  //inline ~Mesh( void );
	 ~Mesh( void );

  void Load(std::string);
	
  //@{
  /// Iterators for the vertex and triangle containers.
  inline TetIt          TetBegin() { return m_tc.begin(); }
  inline TetIt          TetEnd() { return m_tc.end(); }
  inline const_TetIt    TetBegin() const { return m_tc.begin(); }
  inline const_TetIt    TetEnd() const { return m_tc.end(); }
  inline int                    TetSize() const { return static_cast<int>(m_tc.size()); }
  inline VertexIt            VertexBegin() { return m_vc.begin(); }
  inline VertexIt            VertexEnd() { return m_vc.end(); }
  inline VertexIt            BdryVertexBegin() { return m_bv.begin(); }
  inline VertexIt            BdryVertexEnd() { return m_bv.end(); }
  inline const_VertexIt      VertexBegin() const { return m_vc.begin(); }
  inline const_VertexIt      VertexEnd() const { return m_vc.end(); }
  inline const_VertexIt            BdryVertexBegin() const { return m_bv.begin(); }
  inline const_VertexIt            BdryVertexEnd() const { return m_bv.end(); }
  inline FaceIt			BdryFaceBegin() {return m_bf.begin();} 
  inline FaceIt			BdryFaceEnd() { return  m_bf.end(); }
  inline const_FaceIt			BdryFaceBegin() const { return  m_bf.begin(); }
  inline const_FaceIt			BdryFaceEnd() const { return  m_bf.end(); }
  
  inline int              VertexSize() const { return static_cast<int>(m_vc.size()); }
  
  inline int              BdryFaceSize() const { return static_cast<int>(m_bf.size()); }
  inline int              BdryVertexSize() const { return static_cast<int>(m_bv.size()); }
  //@}

  //@{
  /// Add vertices and triangles to the respective containers.
  inline Vertex* addVertex( Vertex* v ) { m_vc.push_back( v ); return v; }
  inline Tet* addTet( Tet* t ) { m_tc.push_back( t ); return t; }
  //@}

  //@{
  /// Direct accessors to the vertex and triangle containers.
  inline Vertex* &V( int i ) { return m_vc[i]; }
  inline Vertex* V( int i ) const { return m_vc[i]; }
  inline Tet* &T( int i ) { return m_tc[i]; }
  inline Tet* T( int i ) const { return m_tc[i]; }
  //@}
  
  // Calls a function to set up topological data & also mark boundaries & set valences
  inline void Init() {
      // connect up triangles
      BuildTetTopology();
      // intialize structures needed by vertex iterators
      //BuildVertexTopology();
  }

 
private:
  /// All vertices.
  VertexCt    m_vc;
  VertexCt    m_bv;
  /// All triangles.
  TetCt   m_tc;
  
  FaceCt m_bf;
	
  void AddAndOrMatch( FaceSet&, Vertex*, Vertex*, Vertex*, Tet* );
  void SetBoundary( FaceSet& );
  //void BuildRing( Tet* t, Vertex *v );
  //void BuildVertexTopology();
  void BuildTetTopology();
  
};


#endif /* __MESH_HPP__ */
