#include <vector>
#include <istream>
#include <iostream>
#include <string>
#include <random>
#include <ctime>        // std::time
#include <cstdlib>     
#include <algorithm> 
#include <cmath> 

#include "boost/tokenizer.hpp"

#include "utils.hpp"
#include "mesh.hpp"
#include "tet.hpp"
#include "vertex.hpp"
#include "face.hpp"

#include <Eigen/Eigenvalues>
#include <Eigen/Cholesky>



	
  

inline double dot_vec(std::vector<double> &x, std::vector<double> &y, int n)
{
    double res=0.0;
    int i;
    for (i = 0; i < n; i++)
    {
        res += x[i] * y[i];
    }
    return res;
}

double mat_vec(std::vector<std::vector<double>> &x, std::vector<double> &y, int n)
{
    double res=0.0;
    int i;
    for (i = 0; i < n; i++)
    {
        //res += x[i] * y[i];
				
    }
    return res;
}

/// This function sets the adjacency pointers for all triangles.
/// It does so by iterating over triangles and calling a helper function
/// that deals with the edges given by the respective vertex pairs.
void Mesh::BuildTetTopology()
{
    // build edge list
    FaceSet fs;
    TetIt ti = TetBegin(), te = TetEnd();
    for( ; ti != te ; ++ti ) {
        Tet* t = *ti;
        AddAndOrMatch( fs, t->a(), t->d(), t->c(), t );
        AddAndOrMatch( fs, t->b(), t->c(), t->d(), t );
        AddAndOrMatch( fs, t->a(), t->c(), t->b(), t );
        AddAndOrMatch( fs, t->a(), t->b(), t->d(), t );
    }
    
    // now we can check for boundary vertices
    SetBoundary( fs );
    
}

void Mesh::SetBoundary( Mesh::FaceSet& fs ) {
	
	VertexSet vs;
    
    FaceSet::iterator fi = fs.begin(), fe = fs.end();
    for (; fi!=fe; ++fi ) {
        Face f = *fi;
        if ( f.Boundary() ) {
            (*fi).a()->MarkBoundary();
            vs.insert((*fi).a());
			(*fi).b()->MarkBoundary();
			vs.insert((*fi).b());
			(*fi).c()->MarkBoundary();
			vs.insert((*fi).c());
			
			m_bf.push_back( new Face(f) );
        }
    }
	
	m_bv = VertexCt(vs.begin(), vs.end() );
    
}

/// Helper function to connect triangles.
void Mesh::AddAndOrMatch( FaceSet& fs, Vertex *a, Vertex *b, Vertex *c, Tet* tn )
{
    // For a given edge which has a triangle adjacent try to find the triangle
    // on the other side. The edge set es knows whether the edge we're trying to
    // insert is already contained in the set.
    Face fn = Face( a, b, c, tn );
    std::pair<FaceSet::iterator, bool> i = fs.insert( fn );
    if( !i.second ){ // the edge was already in the set.
        
        // This means that there is a triangle on the other side,
        // and we can set up the adjacency information.
        Face fo = *i.first; Tet *to = fo.t();
        
        // first we test if the other tet already has some connection
        // set up. This would be fatal!
        assert( !(to->Across( fo.a(), fo.b(), fo.c() ) ) );
        
        // Some more testing to make sure that the two triangles
        // have the same endpoints
        
        // XXX maybe check this somehow... idk.
        //assert( ( en.a() == eo.b() ) && ( en.b() == eo.a() ) );
        // edges must originate from opposite triangles
        assert( to != tn );
        
        // set neighbor information for triangles
        to->Across( fo.a(), fo.b(), fo.c() ) = tn;
        tn->Across( fo.a(), fo.b(), fo.c() ) = to;
        
    }
}

#include <fstream>
void Mesh::Load(std::string fname) {




			//std::cerr<< "Translation values:" << trans1 << trans2 << std::endl;
					
	

	

	//	std::cerr<< "The Matrix is:" << Transx-Trans1_x << std::endl;
	//	std::cerr<< "The Matrix is:" << Transy << std::endl;

	std::ifstream is;
	is.open(fname.c_str());

	std::string line;
	while( std::getline(is, line) ){
		
		boost::char_separator<char> sep(" ");
		boost::tokenizer< boost::char_separator<char> > tok(line, sep);
		boost::tokenizer< boost::char_separator<char> >::iterator ti = tok.begin(), te = tok.end();
		
		if (ti==te) continue; //empty line.
		
       	if( *ti == "v" ){
            
            double x = std::stod(*(++ti));
            double y = std::stod(*(++ti));
            double z = std::stod(*(++ti));


 
			
            Vertex* v = new Vertex( x, y, z );
            assert( v );
            addVertex( v );
           
				}
       
		if( *ti == "t" ) {
            
			Tet *t = new Tet( V(std::stoi(*(++ti))-1),
                             V(std::stoi(*(++ti))-1),
                             V(std::stoi(*(++ti))-1),
                             V(std::stoi(*(++ti))-1) );
            assert( t );
            // get anything the Triangle may want (currently empty)
            addTet( t );
        }
   
    }
    // all the initialization code associated with a fully loaded Mesh, i.e.
    // build triangle and vertex topology, set boundary and set vertex valences
    Init();
}


// Destructor
Mesh::~Mesh( void )
{
    // be a good citizen and clean up after yourself...
    for_each( VertexBegin(), VertexEnd(), Delete<Vertex>() );
    for_each( TetBegin(), TetEnd(), Delete<Tet>() );
}




/*
 /// Set up the vertex topology. This is done by going through  all the triangles
 /// and using the helper function BuildRing to march around each vertex starting
 /// with the current triangle.
 void Mesh::BuildVertexTopology( void )
 {
 // Set the Vertex Face incidence pointers and the connect bits.
	
 // the correctness of this depends on the fact that BuildRing() sets
 // the v->Tri() unless it is already set. Consequently if there are
 // multiple connected components meeting at a vertex only the first
 // one will have the connect bits set for this vertex.
 TriangleIterator ti = TriangleBegin(), te = TriangleEnd();
 for( ; ti != te; ++ti ){
 Triangle *t = *ti;
 BuildRing( t, t->a() );
 BuildRing( t, t->b() );
 BuildRing( t, t->c() );
 }
	
 // is every vertex referenced by somebody? (no lone vertices
 // allowed in a proper 2-manifold)
 VertexIterator vi = VertexBegin(), ve = VertexEnd();
 for( ; vi != ve; ++vi ){
 if( !(*vi)->Tri() ){
 std::cerr <<
 "Fatal error: vertex found whose neighborhood is zero-dimensional\n";
 die();
 }
 }
	
 // Do we have a manifold mesh without vertices whose neighborhood is
 // multiply connected?
 for( ti = TriangleBegin(); ti != te; ti++ ){
 if( !(*ti)->AllVerticesConnected() ){
 std::cerr <<
 "Fatal error: vertex found whose neighborhood is not a (half-)disk\n";
 die();
 }
 }
 }
 
 
 /// March around the vertex to set up the connection to a triangle.
 void Mesh::BuildRing( Triangle* t, Vertex *v )
 {
 if( !v->Tri() ){
 // link it up to current triangle
 v->Tri() = t;
 // now march it CW until boundary or once around
 Triangle *tn = t;
 while( ( tn = tn->CWTriangle( v ) ) && tn != t ) v->Tri() = tn;
 // vertex iterator now valid; use it to set connected bit
 Vertex::TriangleIterator ti = v->TriangleBegin(), te = v->TriangleEnd();
 for( ; ti != te; ++ti ) (*ti)->SetConnected( v );
 }
 }
 
 // Sets up boundary information for each vertex in the mesh
 void Mesh::SetBoundary() {
 
 Mesh::VertexIterator vi = VertexBegin(), ve = VertexEnd();
 for (; vi!=ve; ++vi) (*vi)->MarkBoundary();
 
 }
 
 
 // Set valences for all the vertices in mesh 
 void Mesh::SetValences() {
 
 Mesh::VertexIterator vi = VertexBegin(), ve = VertexEnd();
 for (; vi!=ve; ++vi) (*vi)->SetValence();
 }
 
 

 */




