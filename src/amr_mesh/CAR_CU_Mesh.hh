//----------------------------------*-C++-*----------------------------------//
// CAR_CU_Mesh.hh
// B.T. Adams (bta@lanl.gov)
// 18 May 99
//---------------------------------------------------------------------------//
// @> CAR_CU_Mesh mesh class header file
//---------------------------------------------------------------------------//

#ifndef __mc_CAR_CU_Mesh_hh__
#define __mc_CAR_CU_Mesh_hh__

//===========================================================================//
// class CAR_CU_Mesh - 
//
// Purpose : Continuous Adaptive Refinement Cartesion Unstructured Mesh Class
//
// revision history:
// -----------------
//  0) original (developed from OS_Mesh.hh).
// 
//===========================================================================//

#include "mc/Coord_sys.hh"
#include "Layout.hh"
#include "rng/Sprng.hh"
#include "ds++/SP.hh"
#include "ds++/Assert.hh"
#include <vector>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <map>

namespace rtt_mc 
{
// stl namespaces
using std::vector;
using std::fill;
using std::min_element;
using std::max_element;
using std::ostream;
using std::endl;
using std::string;
using std::multimap;

// draco namespaces
using rtt_rng::Sprng;
using dsxx::SP;
    
class CAR_CU_Mesh
{
  public:
    // The cell centered fields are used only outside of the MESH-type class

    // class definitions of the cell-centered fields: neither of these classes
    // require copy constructors or assignment operators as the SP<> and 
    // vector<> classes can do assignment
    template<class T>
    class CCSF
    {
      private:
	// SP back to CAR_CU_Mesh 
	SP<CAR_CU_Mesh> mesh;
	// data in field, (num_cells)
	vector<T> data;

      public:
	// inline explicit constructor
	inline explicit CCSF(SP<CAR_CU_Mesh>);

	// additional constructors
	inline CCSF(SP<CAR_CU_Mesh>, const vector<T> &);

	// return reference to mesh
	const CAR_CU_Mesh & get_Mesh() const { return * mesh; }

	// subscripting
	const T & operator()(int cell) const { return data[cell-1]; }
	T & operator()(int cell) { return data[cell-1]; }
    };  

    template<class T>
    class CCVF
    {
      private:

	// SP back to CAR_CU_Mesh
	SP<CAR_CU_Mesh> mesh;
	// 2-D field vector, (field width size, num_cells)
	vector< vector<T> > data;

      public:
	// inline explicit constructor (default the leading index vector size
        // to the problem geometry dimension)
	inline explicit CCVF(SP<CAR_CU_Mesh>);

	// inline explicit constructor (arbitry leading index vector size)
	inline explicit CCVF(SP<CAR_CU_Mesh>, int size);

	// additional constructors
	inline CCVF(SP<CAR_CU_Mesh>, const vector<vector<T> > &);

	// return reference to mesh
	const CAR_CU_Mesh & get_Mesh() const { return * mesh; }

	// subscripting
	inline const T & operator()(int, int) const;
	inline T & operator()(int, int);

	// getting a CC vector
	inline vector<T> operator()(int) const;

        // return the size of the CCVF leading index
        int get_size() {return data.size();}
    };  

    // class definitions of the face-centered fields.
    template<class T>
    class FCSF
    {
      private:
	// SP back to CAR_CU_Mesh 
	SP<CAR_CU_Mesh> mesh;
	// data in field, (num_faces)
	vector<T> data;

      public:
	// inline explicit constructor
	inline explicit FCSF(SP<CAR_CU_Mesh>);

	// additional constructors
	inline FCSF(SP<CAR_CU_Mesh>, const vector<T> &);

	// return reference to mesh
	const CAR_CU_Mesh & get_Mesh() const { return * mesh; }

	// subscripting
	inline const T & operator()(int, int) const;
	inline T & operator()(int, int);

 	inline const T & operator()(int) const;
	inline T & operator()(int);

   };  

    template<class T>
    class FCDSF
    {
      private:

	// SP back to CAR_CU_Mesh
	SP<CAR_CU_Mesh> mesh;
	// 2-D field vector, (dimension, num_cells)
	vector< vector<T> > data;

      public:
	// inline explicit constructor
	inline explicit FCDSF(SP<CAR_CU_Mesh>);

	// additional constructors
	inline FCDSF(SP<CAR_CU_Mesh>, const vector<vector<T> > &);

	// return reference to mesh
	const CAR_CU_Mesh & get_Mesh() const { return * mesh; }

	// subscripting
	inline const T & operator()(int, int) const;
	inline T & operator()(int, int);

	// getting a CC vector
	inline vector<T> operator()(int) const;
    };  

    // class definitions of the node-centered fields.
    template<class T>
    class NCSF
    {
      private:
	// SP back to CAR_CU_Mesh 
	SP<CAR_CU_Mesh> mesh;
	// data in field, (num_nodes)
	vector<T> data;

      public:
	// inline explicit constructor (default vector size to the number of
        // nodes)
	inline explicit NCSF(SP<CAR_CU_Mesh>);

	// inline explicit constructor (semi-arbitrary vector size to allow 
        // exclusion of the face-centered nodes)
	inline explicit NCSF(SP<CAR_CU_Mesh>, int size_);

	// additional constructors
	inline NCSF(SP<CAR_CU_Mesh>, const vector<T> &);

	// return reference to mesh
	const CAR_CU_Mesh & get_Mesh() const { return * mesh; }

	// subscripting
        const T & operator()(int node) const { return data[node - 1];}
        T & operator()(int node) { return data[node - 1]; }

        // return the size of the NCSF (allows exclusion of the face-centered
        // nodes)
        int get_size() {return data.size();}
    };  

    template<class T>
    class NCVF
    {
      private:

	// SP back to CAR_CU_Mesh
	SP<CAR_CU_Mesh> mesh;
	// 2-D field vector, (num_nodes, field width size)
	vector< vector<T> > data;

      public:
	// inline explicit constructor (default vector size to the number of
        // nodes by the problem geometry dimension)
	inline explicit NCVF(SP<CAR_CU_Mesh>);

	// inline explicit constructor (semi-arbitrary leading vector size, 
        // default trailing vector size to the number of geometry dimensions
	inline explicit NCVF(SP<CAR_CU_Mesh>, int size_1);

	// inline explicit constructor (semi-arbitrary leading vector size 
        // and arbitrary trailing vector size)
	inline explicit NCVF(SP<CAR_CU_Mesh>, int size_1, int size_2);

	// additional constructors
	inline NCVF(SP<CAR_CU_Mesh>, const vector<vector<T> > &);

	// return reference to mesh
	const CAR_CU_Mesh & get_Mesh() const { return * mesh; }

	// subscripting
	inline const T & operator()(int, int) const;
	inline T & operator()(int, int);

	// getting a NCVF vector
	inline vector<T> operator()(int) const;

        // return the size of the NCSF leading index
        int get_size_1() {return data.size();}

        // return the size of the NCSF trailing index
        int get_size_2() {return data[0].size();}
    };

    template<class K, class T>
    class LNCVF
    {
      private:

	// SP back to CAR_CU_Mesh
	SP<CAR_CU_Mesh> mesh;
	// 2-D field vector, (dimension, num_nodes)
	multimap<K, vector<T> > data;

      public:
	// inline explicit constructor
	inline explicit LNCVF(SP<CAR_CU_Mesh>);

	// additional constructors
	inline LNCVF(SP<CAR_CU_Mesh>, const  multimap<K, vector<T> > &);

	// return reference to mesh
	const CAR_CU_Mesh & get_Mesh() const { return * mesh; }

	// subscripting
	inline const T & operator()(int, int) const;
	inline T & operator()(int, int);

	// getting a NC vector
	inline vector<T> operator()(int) const;
    };  

    // useful typedefs when working with a mesh
    typedef vector<double> CCSF_d;
    typedef vector<int> CCSF_i;
    typedef vector<bool> CCSF_b;

    typedef vector< vector<double> > CCVF_d;
    typedef vector< vector<int> > CCVF_i;

    typedef vector<double> FCSF_d;
    typedef vector<int> FCSF_i;

    typedef vector<double> FCDSF_d;
    typedef vector<int> FCDSF_i;

    typedef vector<double> NCSF_d;
    typedef vector<int> NCSF_i;

    typedef vector< vector<double> > NCVF_d;
    typedef vector< vector<int> > NCVF_i;

    typedef multimap<int, vector<double> > LNCVF_d;

    // temporary typedefs for compiling code until KCC 3.3+ is released
    // (retained herein for compatablity with historical codes).
    typedef CCSF<double> CCSF_double;
    typedef CCSF<int> CCSF_int;
    typedef CCSF<bool> CCSF_bool;
    typedef CCSF<string> CCSF_string;
    typedef CCVF<double> CCVF_double;
    typedef CCVF<int> CCVF_int;

  private:
    // base class reference to a derived coord class
    SP<Coord_sys> coord;
    // layout of mesh
    Layout layout;
    // vertices in mesh
    NCVF_d vertex;
    // cell-pairings of cell to its vertices
    CCVF_i cell_pair;
    // area of surfaces on each dimension
    CCVF_d sur;
    // indicator whether this is a submesh
    bool submesh;
    // indicator that the cell has been refined
    CCSF_b has_kids;
    // indicator for the level of cell refinement, with the original coarse
    // mesh input by the user assigned as zero. 
    CCSF_i generation;

    // private functions

    // compare real values for equality
    bool compReal(const double & low_val, const double & high_val) const;

    // calculate a surface array from the vertices of the mesh
    void calc_surface();

    // private copy and assignment operators; can't copy or assign a mesh
    CAR_CU_Mesh(const CAR_CU_Mesh &);
    CAR_CU_Mesh & operator=(const CAR_CU_Mesh &);

    // Begin_Doc os_mesh-int.tex
    // Begin_Verbatim 

  public:
    // generalized constructor for all mesh types
    CAR_CU_Mesh(SP<Coord_sys>, Layout & , NCVF_d &, CCVF_i &, CCSF_i &,
       bool = false); 

    // member functions used by the CAR_CU_Mesh-dependent classes

    // mesh dimensionality functions

    // give the dimension and begin and end return the beginning and ending
    // coordinate along that dimension
    inline double begin(int) const;
    inline double end(int) const;

    // Problem geometry dimension
    int get_ndim() const { return coord->get_dim(); }

    // return number of cells
    int num_cells() const { return layout.num_cells(); }

    // return total number of nodes
    int num_nodes() const { return vertex[0].size(); }

    // return number of cell-corner nodes
    int num_corner_nodes() const 
    { return cell_pair[num_cells() - 1][static_cast<int>(pow(2.0, 
					         coord->get_dim())) - 1]; }

    // return number of face-centered nodes
    int num_face_nodes() const { return num_nodes() - num_corner_nodes(); }

    // cell dimensionality functions

    // find minimum and maximum dimension of cell
    inline double min(int, int) const;
    inline double max(int, int) const;

    // find centerpoint of cell and width of cell
    inline double pos(int, int) const;
    double dim(int d, int cell) const { return max(d, cell) - min(d, cell); } 

    // diagnostic functions
    void print(ostream &) const;
    void print(ostream &, int) const;

    // End_Verbatim 
    // End_Doc 

    // Begin_Doc car_os_mesh-rint.tex
    // Begin_Verbatim 

    // services required by ALL mesh types used in JAYENNE

    // references to imbedded objects and data required for Parallel_Building
    const Layout & get_Layout() const { return layout; }
    const Coord_sys & get_Coord() const { return * coord; }
    SP<Coord_sys> get_SPCoord() const { return coord; }
    const NCVF_d & get_vertex() const { return vertex; }
    const CCVF_i & get_cell_pair() const { return cell_pair; }
    const CCSF_i & get_generation() const { return generation; }

    // required services for transport; 
    int num_adj(int cell_index, int face_index) const
    { return layout.num_adj(cell_index, face_index);}
    int next_cell(int cell, int face, int adjcell = 1) const 
    { return layout(cell, face, adjcell); }
    int get_cell(const vector<double> &) const;
    double get_db(const vector<double> &, const vector<double> &, int, 
		  int &) const;
    int cell_node(int cell, int node) const
    { return cell_pair[cell - 1][node - 1];}
    int cell_face_centered_node(int cell, int face) const
    {
        int offset = static_cast<int>(pow(2.0,get_ndim()));
	return cell_pair[cell - 1][offset + face - 1];
    }
    inline vector<int> cell_nodes(int cell) const;
    inline vector<int> cell_corner_nodes(int cell) const;
    inline vector<int> cell_face_centered_nodes(int cell) const;
    inline vector<int> cell_face_nodes(int cell, int face) const;
    inline vector<double> get_normal(int, int) const;
    inline vector<double> get_normal_in(int, int) const;
    inline double volume(int) const;
    inline double face_area(int, int) const;
    vector<int> get_surcells(string) const;
    int get_bndface(string, int) const;
    inline vector<double> get_vertex(int node) const;
    inline NCVF_d get_vertices(int, int) const;
    inline NCVF_d get_vertices(int) const;
    inline int get_generation(int) const;
    inline vector<double> sample_pos(int, Sprng &) const;
    inline vector<double> sample_pos(int, Sprng &, vector<double>, 
				     double) const; 
    inline vector<double> sample_pos_on_face(int, int, Sprng &)	const; 
    inline bool check_on_face(vector<double> &, int &, int &) const; 

    // overloaded operators
    bool operator==(const CAR_CU_Mesh &) const;
    bool operator!=(const CAR_CU_Mesh & rhs) const { return !(*this == rhs); }

    // End_Verbatim 
    // End_Doc 
};

//---------------------------------------------------------------------------//
// overloaded operators
//---------------------------------------------------------------------------//

inline ostream & operator<<(ostream & output, const CAR_CU_Mesh & object)
{
    object.print(output);
    return output;
}

//---------------------------------------------------------------------------//
// CAR_CU_Mesh::CCSF inline functions
//---------------------------------------------------------------------------//
// CCSF explicit constructor

template<class T>
inline CAR_CU_Mesh::CCSF<T>::CCSF(SP<CAR_CU_Mesh> mesh_) 
    : mesh(mesh_), data(mesh->num_cells()) 
{
    Require (mesh);
}

//---------------------------------------------------------------------------//
// constructor for automatic initialization

template<class T>
inline CAR_CU_Mesh::CCSF<T>::CCSF(SP<CAR_CU_Mesh> mesh_, 
   const vector<T> & array) : mesh(mesh_), data(array)
{
    // make sure things are kosher
    Ensure (data.size() == mesh->num_cells());
}

//---------------------------------------------------------------------------//
// CAR_CU_Mesh::CCVF inline functions
//---------------------------------------------------------------------------//
// CCVF explicit constructor (default vector leading index size to the 
// problem geometry dimension)

template<class T>
inline CAR_CU_Mesh::CCVF<T>::CCVF(SP<CAR_CU_Mesh> mesh_) 
    : mesh(mesh_), data(mesh->get_ndim())
{
    Require (mesh);

    // initialize data array
    for (int i = 0; i <  mesh->get_ndim(); i++)
	data[i].resize(mesh->num_cells());
}

// CCVF explicit constructor (arbitrary vector leading index size)

template<class T>
inline CAR_CU_Mesh::CCVF<T>::CCVF(SP<CAR_CU_Mesh> mesh_, int size_) 
    : mesh(mesh_), data(size_)
{
    Require (mesh);

    // initialize data array
    for (int i = 0; i <  data.size(); i++)
	data[i].resize(mesh->num_cells());
}

//---------------------------------------------------------------------------//
//---------------------------------------------------------------------------//
// constructor for automatic initialization

template<class T>
inline CAR_CU_Mesh::CCVF<T>::CCVF(SP<CAR_CU_Mesh> mesh_,
				  const vector<vector<T> > & array)
    : mesh(mesh_), data(array)
{
    Require (mesh);

    for (int s = 0; s < data.size(); s++)
	Ensure (data[s].size() == mesh->num_cells());
}

//---------------------------------------------------------------------------//
// constant overloaded ()

template<class T>
inline const T & CAR_CU_Mesh::CCVF<T>::operator()(int dim, int cell) const 
{
    return data[dim-1][cell-1]; 
}

//---------------------------------------------------------------------------//
// assignment overloaded ()

template<class T>
inline T & CAR_CU_Mesh::CCVF<T>::operator()(int dim, int cell)
{
    return data[dim-1][cell-1];
}

//---------------------------------------------------------------------------//
// vector return overload()

template<class T>
inline vector<T> CAR_CU_Mesh::CCVF<T>::operator()(int cell) const
{
    // declare return vector
    vector<T> x;
    
    // loop through dimensions and make return vector for this cell
    for (int i = 0; i < data.size(); i++)
	x.push_back(data[i][cell-1]);

    // return
    Ensure (x.size() == data.size());
    return x;
}

//---------------------------------------------------------------------------//
// CAR_CU_Mesh::FCSF inline functions
//---------------------------------------------------------------------------//
// FCSF explicit constructor

template<class T>
inline CAR_CU_Mesh::FCSF<T>::FCSF(SP<CAR_CU_Mesh> mesh_)
    : mesh(mesh_), data(mesh->num_face_nodes())
{
    Require (mesh);
}

//---------------------------------------------------------------------------//
// constructor for automatic initialization

template<class T>
inline CAR_CU_Mesh::FCSF<T>::FCSF(SP<CAR_CU_Mesh> mesh_, 
			      const vector<T> & array)
    : mesh(mesh_), data(array)
{
    Require (mesh);
    // check things out
    Ensure (data.size() == mesh->num_face_nodes());
}

//---------------------------------------------------------------------------//
// constant overloaded ()

template<class T>
inline const T & CAR_CU_Mesh::FCSF<T>::operator()(int face) const 
{
    return data[face - 1]; 
}

//---------------------------------------------------------------------------//
// constant overloaded ()

template<class T>
inline const T & CAR_CU_Mesh::FCSF<T>::operator()(int cell, int face) const 
{
    int index = mesh->cell_node(cell, face + static_cast<int>(pow(2.0, 
        mesh->get_Coord().get_dim())));
    return data[index - mesh->num_corner_nodes() - 1]; 
}

//---------------------------------------------------------------------------//
// assignment overloaded ()

template<class T>
inline T & CAR_CU_Mesh::FCSF<T>::operator()(int face)
{
    return data[face - 1];
}

//---------------------------------------------------------------------------//
// assignment overloaded ()

template<class T>
inline T & CAR_CU_Mesh::FCSF<T>::operator()(int cell, int face)
{
    int index = mesh->cell_node(cell, face + static_cast<int>(pow(2.0, 
        mesh->get_Coord().get_dim())));
    return data[index - mesh->num_corner_nodes() - 1];
}

//---------------------------------------------------------------------------//
// CAR_CU_Mesh::FCDSF inline functions
//---------------------------------------------------------------------------//
// FCDSF explicit constructor

template<class T>
inline CAR_CU_Mesh::FCDSF<T>::FCDSF(SP<CAR_CU_Mesh> mesh_)
    : mesh(mesh_), data(mesh->num_cells())
{
    Require (mesh);

    // initialize data array
    for (int i = 0; i <mesh->num_cells(); i++)
	data[i].resize(2.0 * mesh->get_Coord().get_dim());
}

//---------------------------------------------------------------------------//
// constructor for automatic initialization

template<class T>
inline CAR_CU_Mesh::FCDSF<T>::FCDSF(SP<CAR_CU_Mesh> mesh_, 
			      const vector<vector<T> > & array)
    : mesh(mesh_), data(array)
{
    // check things out
    Ensure (data.size() == mesh->num_cells());
    for (int i = 0; i < mesh->num_cells(); i++)
	Ensure (data[i].size() == 2.0 * mesh->get_Coord().get_dim());
}

//---------------------------------------------------------------------------//
// constant overloaded ()

template<class T>
inline const T & CAR_CU_Mesh::FCDSF<T>::operator()(int cell, int face) const 
{
    return data[cell - 1][face - 1]; 
}

//---------------------------------------------------------------------------//
// assignment overloaded ()

template<class T>
inline T & CAR_CU_Mesh::FCDSF<T>::operator()(int cell, int face)
{
    return data[cell - 1][face - 1];
}

//---------------------------------------------------------------------------//
// vector return overload()

template<class T>
inline vector<T> CAR_CU_Mesh::FCDSF<T>::operator()(int cell) const
{
    // declare return vector
    vector<T> x;
    
    // loop through faces and make return vector for this cell
    for (int i = 0; i < data[cell-1].size(); i++)
	x.push_back(data[cell-1][i]);

    // return
    Ensure (x.size() == data[cell-1].size());
    return x;
}

//---------------------------------------------------------------------------//
// CAR_CU_Mesh::NCSF inline functions
//---------------------------------------------------------------------------//
// NCSF explicit constructor (default vector size to the number of nodes)

template<class T>
inline CAR_CU_Mesh::NCSF<T>::NCSF(SP<CAR_CU_Mesh> mesh_) 
    : mesh(mesh_), data(mesh->num_nodes()) 
{
    Require (mesh);
}

//---------------------------------------------------------------------------//
// NCSF explicit constructor (semi-arbitrary vector size - allows excluding 
// the face-centered nodes if they are not needed for this field)

template<class T>
inline CAR_CU_Mesh::NCSF<T>::NCSF(SP<CAR_CU_Mesh> mesh_, int size_) 
    : mesh(mesh_), data(size_) 
{
    Require (mesh);

    Ensure (data.size() == mesh->num_nodes() ||
	    data.size() == mesh->num_corner_nodes())
}

//---------------------------------------------------------------------------//
// constructor for automatic initialization

template<class T>
inline CAR_CU_Mesh::NCSF<T>::NCSF(SP<CAR_CU_Mesh> mesh_, 
   const vector<T> & array) : mesh(mesh_), data(array)
{
    Require (mesh);

    Ensure (data.size() == mesh->num_nodes() ||
	    data.size() == mesh->num_corner_nodes())
}

//---------------------------------------------------------------------------//
// CAR_CU_Mesh::NCVF inline functions
//---------------------------------------------------------------------------//
// NCVF explicit constructor (default vector sizes to the number of nodes by
// the number of geometry dimensions)

template<class T>
inline CAR_CU_Mesh::NCVF<T>::NCVF(SP<CAR_CU_Mesh> mesh_) 
    : mesh(mesh_), data(mesh->num_nodes())
{
    Require (mesh);

    // initialize data array
    for (int i = 0; i < size_; i++)
	data[i].resize(mesh->get_ndim());
}

//---------------------------------------------------------------------------//
// NCVF explicit constructor (semi-arbitrary leading vector size, default 
// trailing vector size to the number of geometry dimensions)
template<class T>
inline CAR_CU_Mesh::NCVF<T>::NCVF(SP<CAR_CU_Mesh> mesh_, int size_1) 
    : mesh(mesh_), data(size_1)
{
    Require (mesh);

    // initialize data array
    for (int i = 0; i < size_1; i++)
	data[i].resize(mesh->get_ndim());

    Ensure (size_1 == mesh->num_nodes() ||
	    size_1 == mesh->num_corner_nodes())
}

//---------------------------------------------------------------------------//
// NCVF explicit constructor (semi-arbitrary leading vector size, arbitrary
// trailing vector size)
template<class T>
inline CAR_CU_Mesh::NCVF<T>::NCVF(SP<CAR_CU_Mesh> mesh_, int size_1, 
				  int size_2) 
    : mesh(mesh_), data(size_1)
{
    Require (mesh);

    Ensure (size_1 == mesh->num_nodes() ||
	    size_1 == mesh->num_corner_nodes())

    // initialize data array
    for (int i = 0; i < size_1; i++)
	data[i].resize(size_2);

}

//---------------------------------------------------------------------------//
// constructor for automatic initialization

template<class T>
inline CAR_CU_Mesh::NCVF<T>::NCVF(SP<CAR_CU_Mesh> mesh_, 
			      const vector<vector<T> > & array)
    : mesh(mesh_), data(array)
{
    Require (mesh);

    Ensure (data.size() == mesh->num_nodes() ||
	    data.size() == mesh->num_corner_nodes())
}

//---------------------------------------------------------------------------//
// constant overloaded ()

template<class T>
inline const T & CAR_CU_Mesh::NCVF<T>::operator()(int node, int dim) const 
{
    return data[node-1][dim-1]; 
}

//---------------------------------------------------------------------------//
// assignment overloaded ()

template<class T>
inline T & CAR_CU_Mesh::NCVF<T>::operator()(int node, int dim)
{
    return data[node-1][dim-1];
}

//---------------------------------------------------------------------------//
// vector return overload()

template<class T>
inline vector<T> CAR_CU_Mesh::NCVF<T>::operator()(int node) const
{
    // declare return vector
    vector<T> x;
    
    // loop through dimensions and make return vector for this node
    for (int i = 0; i < data.size(); i++)
	x.push_back(data[node-1][i]);

    // return
    Ensure (x.size() == data.size());
    return x;
}

//---------------------------------------------------------------------------//
// CAR_CU_Mesh::LNCVF inline functions
//---------------------------------------------------------------------------//
// LNCVF explicit constructor

template<class K, class T>
inline CAR_CU_Mesh::LNCVF<K,T>::LNCVF(SP<CAR_CU_Mesh> mesh_)
    : mesh(mesh_), data
{
    Require (mesh);
}

//---------------------------------------------------------------------------//
// constructor for automatic initialization

template<class K, class T>
inline CAR_CU_Mesh::LNCVF<K,T>::LNCVF(SP<CAR_CU_Mesh> mesh_, 
			      const multimap<K, vector<T> > & array)
    : mesh(mesh_), data(array)
{
    // check things out
    Ensure (data.size() == mesh->num_nodes());
    int dim = mesh->get_Coord().get_dim();
    for (multimap<K, vector<T> >::const_iterator niter = data.begin();
	 niter != data.end(); niter++)
	Ensure (niter->second.size() == dim);
}

//---------------------------------------------------------------------------//
// constant overloaded ()

template<class K, class T>
inline const T & CAR_CU_Mesh::LNCVF<K,T>::operator()(int dim, int node) const 
{
    vector<T> local_vector = * data.find(node - 1);
    return local_vector[dim - 1]; 
}

//---------------------------------------------------------------------------//
// assignment overloaded ()

template<class K, class T>
inline T & CAR_CU_Mesh::LNCVF<K,T>::operator()(int dim, int node)
{
    vector<T> local_vector = * data.find(node - 1);
    return local_vector[dim - 1]; 
}

//---------------------------------------------------------------------------//
// vector return overload()

template<class K, class T>
inline vector<T> CAR_CU_Mesh::LNCVF<K,T>::operator()(int node) const
{
    return data.find(node - 1)->second;
}

//---------------------------------------------------------------------------//
// CAR_CU_Mesh inline functions
//---------------------------------------------------------------------------//

inline double CAR_CU_Mesh::begin(int d) const 
{
    // find the minimum surface for d over the whole mesh
    return * min_element(vertex[d - 1].begin(), vertex[d - 1].end()); 
}

//---------------------------------------------------------------------------//

inline double CAR_CU_Mesh::end(int d) const 
{
    // find the maximum surface for d over the whole mesh
    return * max_element(vertex[d - 1].begin(), vertex[d - 1].end()); 
}

//---------------------------------------------------------------------------//

inline double CAR_CU_Mesh::pos(int d, int cell) const
{
    // find center position of cell

    // set return value
    double return_pos = 0.0;

    // loop over all vertices and take average value to get the center
    // point 
    for (int i = 0; i < pow(2.0,coord->get_dim()); i++)
	return_pos += vertex[d-1][cell_pair[cell-1][i]-1];

    // return value
    return return_pos /  pow(2.0,coord->get_dim());     
}

//---------------------------------------------------------------------------//

inline double CAR_CU_Mesh::min(int d, int cell) const 
{	
    // find minimum dimension along d of cell

    // loop over all vertices and find the minimum
    double minimum = vertex[d-1][cell_pair[cell-1][0]-1];
    for (int i = 1; i < pow(2.0,coord->get_dim()); i++)
    {
	double point = vertex[d-1][cell_pair[cell-1][i]-1];

	// update the minimum value point
	if (point < minimum)
	    minimum = point;
    }
	
    // return minimum dimension
    return minimum;
}

//---------------------------------------------------------------------------//

inline double CAR_CU_Mesh::max(int d, int cell) const
{
    // find maximum dimension of cell

    // loop over all vertices and find the maximum
    double maximum = vertex[d-1][cell_pair[cell-1][0]-1];
    for (int i = 1; i < pow(2.0,coord->get_dim()); i++)
    {
	double point = vertex[d-1][cell_pair[cell-1][i]-1];

	// update the maximum value point
	if (point > maximum)
	    maximum = point;
    }

    // return maximum dimension
    return maximum;
}

//---------------------------------------------------------------------------//

inline double CAR_CU_Mesh::volume(int cell) const 
{
    // calculate volume of cell

    // loop through dimensions and get volume
    double volume = 1.0;
    for (int d = 1; d <= coord->get_dim(); d++)
	volume *= dim(d, cell);

    // return volume
    return volume;
}

//---------------------------------------------------------------------------//

inline double CAR_CU_Mesh::face_area(int cell, int face) const 
{
    // calculate area of face on cell

    // loop through dimensions and multiply off-dimension widths
    double face_area = 1.0;
    int dim_face_on = abs(face - coord->get_dim()) + 1 - 
                      face/(coord->get_dim() + 1);

    for (int d = 1; d <= coord->get_dim(); d++)
    {
	if (d != dim_face_on)
	    face_area *= dim(d, cell);
    }

    // return face_area
    return face_area;
}

//---------------------------------------------------------------------------//

inline vector<int> CAR_CU_Mesh::cell_nodes(int cell) const
{
    // Return the set of nodes that make up a cell, including both the corner
    // nodes and the face-centered nodes.
    int nnodes = pow(2.0,coord->get_dim()) + 2.0 * coord->get_dim();
    vector<int> node_set(nnodes);
    for (int node = 0; node < nnodes; node++)
        node_set[node] = cell_pair[cell - 1][node];

    return node_set;
}

//---------------------------------------------------------------------------//

inline vector<int> CAR_CU_Mesh::cell_corner_nodes(int cell) const
{
    // Return the set of corner nodes that make up a cell.
    int nnodes = pow(2.0,coord->get_dim());
    vector<int> node_set(nnodes);
    for (int node = 0; node < nnodes; node++)
        node_set[node] = cell_pair[cell - 1][node];

    return node_set;
}

//---------------------------------------------------------------------------//

inline vector<int> CAR_CU_Mesh::cell_face_centered_nodes(int cell) const
{
    // Return the set of face-centered nodes for a cell.
    int nnodes = pow(2.0,coord->get_dim()) + 2.0 * coord->get_dim();
    vector<int> node_set(2.0 * coord->get_dim());
    for (int node = pow(2.0,coord->get_dim()); node < nnodes; node++)
        node_set[node - pow(2.0,coord->get_dim())] = cell_pair[cell - 1][node];

    return node_set;
}

//---------------------------------------------------------------------------//

inline vector<int> CAR_CU_Mesh::cell_face_nodes(int cell, int face) const
{
    // Return the set of nodes that make up a cell face
    vector<int> node_set(pow(2.0,coord->get_dim() - 1));
    // Correlations are based on cell faces indexed to start at zero
    --face;

    int node = face/3 + face/4 + 2 * (face/5);
    node_set[0] = cell_pair[cell - 1][node];

    node = ((face+1)%2) * (1 + face/2) + (face%2) * (4 + face/2);
    node_set[1] = cell_pair[cell - 1][node];

    if (coord->get_dim() == 3)
    {
        node = ((face+1)%2) * (3 * (1 + face/2) - 2 * (face/4)) +
	       (face%2) * (5 + 2 * (face/3));
        node_set[2] = cell_pair[cell - 1][node];

        node = ((face+1)%2) * (face + 2) + (face%2) * face;
        node_set[3] = cell_pair[cell - 1][node];
    }
    // Reset the face value to the original input
    ++face;

    return node_set;
}

//---------------------------------------------------------------------------//

inline vector<double> CAR_CU_Mesh::get_normal(int cell, int face) const
{
    // CAR_CU_Meshes do not require any functionality from Coord_sys to 
    // calculate the outward normal, do simple return

    // normal always has 3 components, use Get_sdim()
    vector<double> normal(coord->get_sdim(), 0.0);
	
    // calculate normal based on face, (-z, -y, -x, +x, +y, +z), only
    // one coordinate is non-zero    
    int index = abs(face - coord->get_dim()) - face/(coord->get_dim() + 1);
    if (face <= coord->get_dim())
        normal[index] = -1.0;
    else
        normal[index] = 1.0;

    // return the normal
    return normal;
}

//---------------------------------------------------------------------------//

inline vector<double> CAR_CU_Mesh::get_normal_in(int cell, int face) const
{
    // CAR_CU_Meshes do not require any functionality from Coord_sys to 
    // calculate the inward normal, do simple return

    // normal always has 3 components, use Get_sdim()
    vector<double> normal(coord->get_sdim(), 0.0);
	
    // calculate normal based on face, (-z, -y, -x, +x, +y, +z), only
    // one coordinate is non-zero    
    int index = abs(face - coord->get_dim()) - face/(coord->get_dim() + 1);
    if (face > coord->get_dim())
        normal[index] = -1.0;
    else
        normal[index] = 1.0;

    // return the normal
    return normal;
}

//---------------------------------------------------------------------------//
// Return the node vertex

inline vector<double> CAR_CU_Mesh::get_vertex(int node) const
{
    // determine the vertices along a cell-face

    // return vertices
    vector<double> ret_vert(coord->get_dim());

    for (int d = 0; d < coord->get_dim(); d++)
        ret_vert[d] = vertex[d][node-1];

    // return vector of vertices
    return ret_vert;
}

//---------------------------------------------------------------------------//
// calculate the vertices bounding a cell face

inline CAR_CU_Mesh::NCVF_d CAR_CU_Mesh::get_vertices(int cell, int face) const
{
    // determine the vertices along a cell-face

    // return vertices
    NCVF_d ret_vert(coord->get_dim());

    // determine axis dimension of surface (x=1, y=2, z=3)
    int axis = abs(face - coord->get_dim()) + 1 - face/(coord->get_dim() + 1);
    double plane;
    if (face <= coord->get_dim())
	plane = min(axis, cell);
    else
	plane = max(axis, cell);

    // loop over vertices in cell and get the vertices that are in the plane
    for (int i = 0; i < pow(2.0,coord->get_dim()); i++)
	if (plane == vertex[axis-1][cell_pair[cell-1][i]-1])
	    for (int d = 0; d < coord->get_dim(); d++)
		ret_vert[d].push_back(vertex[d][cell_pair[cell-1][i]-1]);

    // return vector of vertices
    return ret_vert;
}

//---------------------------------------------------------------------------//
// calculate the vertices bounding a cell

inline CAR_CU_Mesh::NCVF_d CAR_CU_Mesh::get_vertices(int cell) const
{
    // determine the vertices bounding a cell
    
    // return vertices
    NCVF_d ret_vert(coord->get_dim());

    // loop over cell vertices and build the cell vertices
    for (int i = 0; i < pow(2.0,coord->get_dim()); i++)
	for (int d = 0; d < coord->get_dim(); d++)
	    ret_vert[d].push_back(vertex[d][cell_pair[cell-1][i]-1]);

    // return vector of vertices
    return ret_vert;
}

//---------------------------------------------------------------------------//
// determine the generation level of a cell

inline int CAR_CU_Mesh::get_generation(int cell) const
{
    return generation[cell-1];
}

//---------------------------------------------------------------------------//
// sample the position uniformly in a cell

inline vector<double> CAR_CU_Mesh::sample_pos(int cell, Sprng & random) const
{
    // assign minimums and maximums for cell dimensions
    vector<double> vmin(coord->get_dim());
    vector<double> vmax(coord->get_dim());

    for (int d = 1; d <= coord->get_dim(); d++)
    {
	vmin[d-1] = min(d, cell);
	vmax[d-1] = max(d, cell);
    }

    // use coord_sys to sample the location
    vector<double> r = coord->sample_pos(vmin, vmax, random);

    // return position vector
    return r;
}

//---------------------------------------------------------------------------//
// sample the position in a cell given a tilt (or other slope-like function)

inline vector<double> CAR_CU_Mesh::sample_pos(int cell, Sprng & random, 
					  vector<double> slope, 
					  double center_pt) const
{
    // assign minimums and maximums for cells dimensions
    vector<double> vmin(coord->get_dim());
    vector<double> vmax(coord->get_dim());

    for (int d = 1; d <= coord->get_dim(); d++)
    {
	vmin[d-1] = min(d, cell);
	vmax[d-1] = max(d, cell);
    }

    // use coord_sys to sample the location
    vector<double> r = coord->sample_pos(vmin, vmax, random, slope,
					 center_pt);

    // return position vector
    return r;
}

//---------------------------------------------------------------------------//
// sample a position on a face

inline vector<double> CAR_CU_Mesh::sample_pos_on_face(int cell, int face, 
						  Sprng & random) const
{
    // assign minimums and maximums for cell dimensions
    vector<double> vmin(coord->get_dim());
    vector<double> vmax(coord->get_dim());

    for (int d = 1; d <= coord->get_dim(); d++)
    {
	vmin[d-1] = min(d, cell);
	vmax[d-1] = max(d, cell);
    }

    // use coord_sys to sample the location
    vector<double> r = coord->sample_pos_on_face(vmin, vmax, face, random);

    // return position vector
    return r;
}

//---------------------------------------------------------------------------//
// check if a point lies on a face

inline bool CAR_CU_Mesh::check_on_face(vector<double> & pos, int & cell, 
				       int & face) const
{
    int off_face = abs(face - coord->get_dim()) + 1 - 
                       face/(1 + coord->get_dim());
    int on_face = 1;

    for (int d = 1; d <= coord->get_dim(); d++)
    {
	if (d != off_face)
	{
	    double side_max = max(d, cell);
	    double side_min = min(d, cell);
	    // multiply integer by zero if this coordinate is not on the face
	    on_face *= (pos[d-1] >= side_min && pos[d-1] <= side_max);
	}
    }
    return (on_face != 0);

}

} // end namespace rtt_mc

#endif                          // __mc_CAR_CU_Mesh_hh__

//---------------------------------------------------------------------------//
//                              end of mc/CAR_CU_Mesh.hh
//---------------------------------------------------------------------------//
