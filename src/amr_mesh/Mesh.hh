//----------------------------------*-C++-*----------------------------------//
// Mesh.hh
// B.T. Adams (bta@lanl.gov)
// 18 May 99
/*! 
 * \file   amr_mesh/Mesh.hh
 * \author B.T. Adams
 * \date   Tue May 18 10:33:26 1999
 * \brief  Header file for CAR_CU_Mesh class library.
 */
//---------------------------------------------------------------------------//
// @> CAR_CU_Mesh mesh class header file
//---------------------------------------------------------------------------//

#ifndef __amr_CAR_CU_Mesh_hh__
#define __amr_CAR_CU_Mesh_hh__

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

#include "Layout.hh"
#include "ds++/SP.hh"
#include "ds++/Assert.hh"
#include <vector>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include <map>

namespace rtt_amr 
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
using rtt_dsxx::SP;

/*!
 * \brief  The continuous adaptive refinement (CAR) Cartesion unstructured 
 *         (CU) Mesh class provides adaptive mesh refinement (amr) capability
 *         for use in transport codes.
 *
 *\sa The CAR_CU_Mesh class is typically instantiated using the associated 
 *    CAR_CU_Interface and CAR_CU_Builder class objects. The class contains
 *    numerous nested mesh field classes that can be instantiated and accessed
 *    as needs demand. The \ref amr_overview presents a summary of the 
 *    capabilities and intended usage of this mesh class.
 */     
class CAR_CU_Mesh
{
  public:
    // The cell centered fields are used only outside of the MESH-type class

    // class definitions of the cell-centered fields: neither of these classes
    // require copy constructors or assignment operators as the SP<> and 
    // vector<> classes can do assignment
/*!
 * \brief CAR_CU_Mesh templated nested mesh field class for cell-centered 
 *        scalar field (CCSF) data.
 */
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
/*!
 * \brief Constructs an uninitialized CAR_CU_Mesh cell-centered scalar field 
 *        (CCSF) class object sized to the number of cells in the mesh.
 * \param mesh Smart pointer to the current CAR_CU_Mesh class object.
 */
	inline explicit CCSF(SP<CAR_CU_Mesh> mesh);

	// additional constructors
/*!
 * \brief Constructs an initialized CAR_CU_Mesh cell-centered scalar field 
 *        (CCSF) class object sized to the number of cells in the mesh.
 * \param mesh Smart pointer to the current CAR_CU_Mesh class object.
 * \param array CCSF initialization data of type T (must be sized to the 
 *              number of cells in the mesh).
 */
	inline CCSF(SP<CAR_CU_Mesh> mesh, const vector<T> & array);

	// return reference to mesh
/*!
 * \brief Returns the CAR_CU_Mesh class object associated with the current 
 *        CCSF nested mesh field class object. 
 * \return Reference to the mesh.
 */
	const CAR_CU_Mesh & get_Mesh() const { return * mesh; }

	// subscripting
/*!
 * \brief Overloaded operator to return the CCSF data value for the specified 
 *        cell.
 * \param cell Cell number.
 * \return CCSF data value of type T.
 */
	const T & operator()(int cell) const { return data[cell]; }
/*!
 * \brief Overloaded operator to assign a data value to the specified CCSF
 *        cell.
 * \param cell Cell number.
 * \return CCSF data value of type T.
 */
	T & operator()(int cell) { return data[cell]; }
    };  

/*!
 * \brief CAR_CU_Mesh templated nested mesh field class for cell-centered 
 *        vector field (CCVF) data.
 */
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
/*!
 * \brief Constructs an uninitialized CAR_CU_Mesh cell-centered vector field 
 *        (CCVF) class object with the vector leading index size defaulted to
 *        the number of spatial dimensions and the second index sized to the 
 *        number of cells in the mesh.
 * \param mesh Smart pointer to the current CAR_CU_Mesh class object.
 */
	inline explicit CCVF(SP<CAR_CU_Mesh> mesh);

	// inline explicit constructor (arbitry leading index vector size)
/*!
 * \brief Constructs an uninitialized CAR_CU_Mesh cell-centered vector field 
 *        (CCVF) class object with the vector leading index size arbitrary 
 *        and the second index sized to the number of cells in the mesh.
 * \param mesh Smart pointer to the current CAR_CU_Mesh class object.
 * \param size Vector leading index size.
 */
	inline explicit CCVF(SP<CAR_CU_Mesh> mesh, int size);

	// additional constructors
/*!
 * \brief Constructs an initialized CAR_CU_Mesh cell-centered vector field 
 *        (CCVF) class object with the vector leading index arbitrary and the 
 *        second index sized to the number of cells in the mesh.
 * \param mesh Smart pointer to the current CAR_CU_Mesh class object.
 * \param array CCVF initialization data of type T (second index must be sized
 *              to the number of cells in the mesh).
 */
	inline CCVF(SP<CAR_CU_Mesh> mesh, const vector<vector<T> > & array);

	// return reference to mesh
/*!
 * \brief Returns the CAR_CU_Mesh class object associated with the current 
 *        CCVF nested mesh field class object. 
 * \return Reference to the mesh.
 */
	const CAR_CU_Mesh & get_Mesh() const { return * mesh; }

	// subscripting
/*!
 * \brief Overloaded operator to return the CCVF data value for the specified 
 *        leading index element and cell. 
 * \param dim Leading index element number.
 * \param cell Cell number.
 * \return CCVF data value of type T.
 */
	inline const T & operator()(int dim, int cell) const;
/*!
 * \brief Overloaded operator to assign a data value to the specified CCVF
 *        leading index element and cell. 
 * \param dim Leading index element number.
 * \param cell Cell number.
 * \return CCVF data value of type T.
 */
	inline T & operator()(int dim, int cell);

	// getting a CCVF vector
/*!
 * \brief Overloaded operator to return all of the CCVF leading index data 
 *        values for the specified cell. 
 * \param cell Cell number.
 * \return CCVF data values of type T.
 */
	inline vector<T> operator()(int cell) const;

        // return the size of the CCVF leading index
/*!
 * \brief Returns the size of the CCVF leading index. 
 * \return CCVF leading index size.
 */
        int get_size() {return data.size();}
    };  

    // class definitions of the face-centered fields.
/*!
 * \brief CAR_CU_Mesh templated nested mesh field class for face-centered 
 *        scalar field (FCSF) data.
 */
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
/*!
 * \brief Constructs an uninitialized CAR_CU_Mesh face-centered scalar field 
 *        (FCSF) class object sized to the number of unique cell faces in the 
 *        mesh.
 * \param mesh Smart pointer to the current CAR_CU_Mesh class object
 */
	inline explicit FCSF(SP<CAR_CU_Mesh> mesh);

	// additional constructors
/*!
 * \brief Constructs an initialized CAR_CU_Mesh face-centered scalar field 
 *        (FCSF) class object sized to the number of unique cell faces in 
 *        the mesh.
 * \param mesh Smart pointer to the current CAR_CU_Mesh class object
 * \param array FCSF initialization data of type T (must be sized to the 
 *              number of unique cell faces in the mesh).
 */
	inline FCSF(SP<CAR_CU_Mesh> mesh, const vector<T> & array);

	// return reference to mesh
/*!
 * \brief Returns the CAR_CU_Mesh class object associated with the current 
 *        FCSF nested mesh field class object. 
 * \return Reference to the mesh.
 */
	const CAR_CU_Mesh & get_Mesh() const { return * mesh; }

	// subscripting
/*!
 * \brief Overloaded operator to return the FCSF data value for the specified 
 *        cell face.
 * \param cell Cell number.
 * \param face Face number.
 * \return FCSF data value of type T.
 */
	inline const T & operator()(int cell, int face) const;
/*!
 * \brief Overloaded operator to assign a data value to the specified FCSF
 *        cell face.
 * \param cell Cell number.
 * \param face Face number.
 * \return FCSF data value of type T.
 */
	inline T & operator()(int cell, int face);

/*!
 * \brief Overloaded operator to return the FCSF data value for the specified 
 *        unique cell face.
 * \param face Unique cell face number.
 * \return FCSF data value of type T.
 */
 	inline const T & operator()(int face) const;
/*!
 * \brief Overloaded operator to assign a data value to the specified FCSF
 *        unique cell face.
 * \param face Unique cell face number.
 * \return FCSF data value of type T.
 */
	inline T & operator()(int face);

   };

/*!
 * \brief CAR_CU_Mesh templated nested mesh field class for face-centered 
 *        discontinuous scalar field (FCDSF) data.
 */
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
/*!
 * \brief Constructs an uninitialized CAR_CU_Mesh face-centered discontinuous
 *        scalar field (FCDSF) class object with the leading index sized to 
 *        the number of cells in the mesh and the trailing index sized to the 
 *        number of faces per cell.
 * \param mesh Smart pointer to the current CAR_CU_Mesh class object
 */
	inline explicit FCDSF(SP<CAR_CU_Mesh> mesh);

	// additional constructors
/*!
 * \brief Constructs an initialized CAR_CU_Mesh face-centered discontinuous
 *        scalar field (FCDSF) class object with the leading index sized to
 *        the number of cells in the mesh and the trailing index sized to
 *        the number of faces per cell.
 * \param mesh Smart pointer to the current CAR_CU_Mesh class object
 * \param array FCDSF initialization data of type T (must be sized with the 
 *        leading index equal to the number of cells in the mesh and the 
 *        trailing index sized to the number of faces per cell).
 */
	inline FCDSF(SP<CAR_CU_Mesh> mesh, const vector<vector<T> > & array);

	// return reference to mesh
/*!
 * \brief Returns the CAR_CU_Mesh class object associated with the current 
 *        FCDSF nested mesh field class object. 
 * \return Reference to the mesh.
 */
	const CAR_CU_Mesh & get_Mesh() const { return * mesh; }

	// subscripting
/*!
 * \brief Overloaded operator to return the FCDSF data value for the specified 
 *        cell face.
 * \param cell Cell number.
 * \param face Face number.
 * \return FCDSF data value of type T.
 */
	inline const T & operator()(int cell, int face) const;
/*!
 * \brief Overloaded operator to assign a data value to the specified FCDSF
 *        cell face.
 * \param cell Cell number.
 * \param face Face number.
 * \return FCDSF data value of type T.
 */
	inline T & operator()(int cell, int face);

	// getting a FCDSF vector
/*!
 * \brief Overloaded operator to return all of the FCDSF face values for the 
 *        specified cell. 
 * \param cell Cell number.
 * \return FCDSF data values of type T.
 */
	inline vector<T> operator()(int cell) const;
    };  

/*!
 * \brief CAR_CU_Mesh templated nested mesh field class for face-centered 
 *        vector field (FCVF) data.
 */
    template<class T>
    class FCVF
    {
      private:
	// SP back to CAR_CU_Mesh 
	SP<CAR_CU_Mesh> mesh;
	// data in field, (num_faces, dim)
	vector<vector<T> > data;

      public:
	// inline explicit constructor
/*!
 * \brief Constructs an uninitialized CAR_CU_Mesh face-centered vector field 
 *        (FCVF) class object with the vector leading index sized to the 
 *        number of unique cell faces in the mesh and the second index size 
 *        defaulted to the number of spatial dimensions.
 * \param mesh Smart pointer to the current CAR_CU_Mesh class object.
 */
	inline explicit FCVF(SP<CAR_CU_Mesh> mesh);

	// inline explicit constructor
/*!
 * \brief Constructs an uninitialized CAR_CU_Mesh face-centered vector field 
 *        (FCVF) class object with the vector leading index sized to the 
 *        number of unique cell faces in the mesh and the second index size
 *        arbitrary.
 * \param mesh Smart pointer to the current CAR_CU_Mesh class object.
 * \param vec_size Size of the vector second index.
 */
	inline explicit FCVF(SP<CAR_CU_Mesh> mesh, int vec_size);

	// additional constructors
	// inline explicit constructor
/*!
 * \brief Constructs an initialized CAR_CU_Mesh face-centered vector field 
 *        (FCVF) class object with the vector leading index sized to the 
 *        number of unique cell faces in the mesh and the second index size
 *        arbitrary.
 * \param mesh Smart pointer to the current CAR_CU_Mesh class object.
 * \param array FCVF initialization data of type T (the leading index must be
 *              sized to the number of unique cell faces in the mesh).
 */
	inline FCVF(SP<CAR_CU_Mesh> mesh, const vector<vector<T> > & array);

	// return reference to mesh
/*!
 * \brief Returns the CAR_CU_Mesh class object associated with the current 
 *        FCVF nested mesh field class object. 
 * \return Reference to the mesh.
 */
	const CAR_CU_Mesh & get_Mesh() const { return * mesh; }

	// subscripting
/*!
 * \brief Overloaded operator to return the FCVF data value for the specified 
 *        cell face trailing index element number. 
 * \param cell Cell number.
 * \param face Face number.
 * \param dim Trailing index element number.
 * \return FCVF data value of type T.
 */
	inline const T & operator()(int cell, int face, int dim) const;
/*!
 * \brief Overloaded operator to assign a data value to the specified FCVF 
 *        cell face trailing index element number. 
 * \param cell Cell number.
 * \param face Face number.
 * \param dim Trailing index element number.
 * \return FCVF data value of type T.
 */
	inline T & operator()(int cell, int face, int dim);

/*!
 * \brief Overloaded operator to return all of the FCVF data values for the 
 *        specified cell face. 
 * \param cell Cell number.
 * \param face Face number.
 * \return FCVF data values of type T.
 */
 	inline const vector<T> & operator()(int cell, int face) const;
/*!
 * \brief Overloaded operator to assign all of the data values for the 
 *        specified FCVF cell face. 
 * \param cell Cell number.
 * \param face Face number.
 * \return FCVF data values of type T.
 */
	inline vector<T> & operator()(int cell, int face);

	// return reference to the data
/*!
 * \brief Overloaded operator to return an entire FCVF. 
 * \return FCVF with data type T.
 */
	inline const vector<vector<T> > & operator()() const;
/*!
 * \brief Overloaded operator to assign an entire FCVF. 
 * \return FCVF with data type T.
 */
	inline vector<vector<T> > & operator()();

        // return the size of the FCVF trailing index
/*!
 * \brief Returns the size of the FCVF trailing index. 
 * \return FCVF trailing index size.
 */
        int get_size() {return data[0].size();}
    };  

    // class definitions of the node-centered fields.
/*!
 * \brief CAR_CU_Mesh templated nested mesh field class for node-centered 
 *        scalar field (NCSF) data.
 */
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
/*!
 * \brief Constructs an uninitialized CAR_CU_Mesh node-centered scalar field 
 *        (NCSF) class object sized to the number of nodes (i.e., corner and
 *        face-centered) in the mesh.
 * \param mesh Smart pointer to the current CAR_CU_Mesh class object
 */
	inline explicit NCSF(SP<CAR_CU_Mesh> mesh);

	// inline explicit constructor (semi-arbitrary vector size to allow 
        // exclusion of the face-centered nodes)
/*!
 * \brief Constructs an uninitialized CAR_CU_Mesh node-centered scalar field 
 *        (NCSF) class object of semi-arbitrary size (must be either the total
 *        number of nodes in the mesh or the number of corner nodes in the
 *        mesh).
 * \param mesh Smart pointer to the current CAR_CU_Mesh class object
 * \param size NCSF size (must be either the total number of nodes in the mesh 
 *             or the number of corner nodes in the mesh).
 */
	inline explicit NCSF(SP<CAR_CU_Mesh> mesh, int size);

	// additional constructors
/*!
 * \brief Constructs an initialized CAR_CU_Mesh node-centered scalar field 
 *        (NCSF) class object of semi-arbitrary size (must be either the 
 *        total number of nodes in the mesh or the number of corner nodes in 
 *        the mesh).
 * \param mesh Smart pointer to the current CAR_CU_Mesh class object
 * \param array NCSF initilization data of type T (must be sized to either 
 *              the total number of nodes in the mesh or the number of corner 
 *              nodes in the mesh).
 */
	inline NCSF(SP<CAR_CU_Mesh> mesh, const vector<T> & array);

	// return reference to mesh
/*!
 * \brief Returns the CAR_CU_Mesh class object associated with the current 
 *        NCSF nested mesh field class object. 
 * \return Reference to the mesh.
 */
	const CAR_CU_Mesh & get_Mesh() const { return * mesh; }

	// subscripting
/*!
 * \brief Overloaded operator to return the NCSF data value for the specified 
 *        node. 
 * \param node Node number.
 * \return NCSF data value of type T.
 */
        const T & operator()(int node) const { return data[node];}
/*!
 * \brief Overloaded operator to assign a data value to the specified NCSF
 *        node. 
 * \param node Node number.
 * \return NCSF data value of type T.
 */
        T & operator()(int node) { return data[node]; }

        // return the size of the NCSF (allows exclusion of the face-centered
        // nodes)
/*!
 * \brief Returns the size of the NCSF. 
 * \return NCSF size.
 */
        int get_size() {return data.size();}
    };  

/*!
 * \brief CAR_CU_Mesh templated nested mesh field class for node-centered 
 *        vector field (NCVF) data.
 */
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
/*!
 * \brief Constructs an uninitialized CAR_CU_Mesh node-centered vector field 
 *        (NCVF) class object with the leading index sized to the number of 
 *        nodes (i.e., corner and face-centered) in the mesh and the trailing
 *        index size defaulted to the number of spatial dimensions.
 * \param mesh Smart pointer to the current CAR_CU_Mesh class object
 */
	inline explicit NCVF(SP<CAR_CU_Mesh> mesh);

	// inline explicit constructor (semi-arbitrary leading vector size, 
        // default trailing vector size to the number of geometry dimensions
/*!
 * \brief Constructs an uninitialized CAR_CU_Mesh node-centered vector field 
 *        (NCVF) class object with the leading index semi-arbitrary (i.e., 
 *        sized to the either total number of nodes in the mesh or the number
 *        of corner nodes in the mesh) and the trailing index size defaulted
 *        to the number of spatial dimensions.
 * \param mesh Smart pointer to the current CAR_CU_Mesh class object
 * \param size_1 NCVF leading index size (must be either the total number of 
 *               nodes in the mesh or the number of corner nodes in the mesh).
 */
	inline explicit NCVF(SP<CAR_CU_Mesh> mesh, int size_1);

	// inline explicit constructor (semi-arbitrary leading vector size 
        // and arbitrary trailing vector size)
/*!
 * \brief Constructs an uninitialized CAR_CU_Mesh node-centered vector field 
 *        (NCVF) class object with the leading index semi-arbitrary (i.e., 
 *        sized to the either the total number of nodes in the mesh or the 
 *        number of corner nodes in the mesh) and the trailing index size 
 *        arbitrary.
 * \param mesh Smart pointer to the current CAR_CU_Mesh class object
 * \param size_1 NCVF leading index size (must be either the total number of 
 *               nodes in the mesh or the number of corner nodes in the mesh).
 * \param size_2 NCVF trailing index size.
 */
	inline explicit NCVF(SP<CAR_CU_Mesh> mesh, int size_1, int size_2);

	// additional constructors
/*!
 * \brief Constructs an initialized CAR_CU_Mesh node-centered vector field 
 *        (NCVF) class object with the leading index semi-arbitrary (i.e., 
 *        sized to the either total number of nodes in the mesh or the number
 *        of corner nodes in the mesh) and the trailing index size arbitrary.
 * \param mesh Smart pointer to the current CAR_CU_Mesh class object
 * \param array NCVF initialization data of type T (leading index size must 
 *              be either the total number of nodes in the mesh or the number 
 *              of corner nodes in the mesh).
 */
	inline NCVF(SP<CAR_CU_Mesh> mesh, const vector<vector<T> > & array);

	// return reference to mesh
/*!
 * \brief Returns the CAR_CU_Mesh class object associated with the current 
 *        NCVF nested mesh field class object. 
 * \return Reference to the mesh.
 */
	const CAR_CU_Mesh & get_Mesh() const { return * mesh; }

	// subscripting
/*!
 * \brief Overloaded operator to return the NCVF data value for the specified 
 *        node trailing index element number. 
 * \param node Node number.
 * \param dim Trailing index element number.
 * \return NCVF data value of type T.
 */
	inline const T & operator()(int node, int dim) const;
/*!
 * \brief Overloaded operator to assign a data value to the specified NCVF
 *        node trailing index element number. 
 * \param node Node number.
 * \param dim Trailing index element number.
 * \return NCVF data value of type T.
 */
	inline T & operator()(int node, int dim);

	// getting a NCVF vector
/*!
 * \brief Overloaded operator to return all of the NCVF data values for the 
 *        specified node. 
 * \param node Node number.
 * \return NCVF data values of type T.
 */
	inline vector<T> operator()(int node) const;

        // return the size of the NCVF leading index
/*!
 * \brief Returns the size of the NCVF leading index. 
 * \return NCVF leading index size.
 */
        int get_size_1() {return data.size();}

        // return the size of the NCVF trailing index
/*!
 * \brief Returns the size of the NCVF trailing index. 
 * \return NCVF trailing index size.
 */
        int get_size_2() {return data[0].size();}
    };

/*!
 * \brief CAR_CU_Mesh templated nested mesh field class for localized 
 *        node-centered vector field (LNCVF) data (not currently implemented).
 */
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
/*!
 * \brief Data structure for cell-centered scalar field double data (not a
 *        CAR_CU_Mesh nested mesh field class).
 */
    typedef vector<double> ccsf_d;
/*!
 * \brief Data structure for cell-centered scalar field integer data (not a 
 *        CAR_CU_Mesh nested mesh field class).
 */
    typedef vector<int> ccsf_i;
/*!
 * \brief Data structure for cell-centered scalar field boolean data (not a 
 *        CAR_CU_Mesh nested mesh field class).
 */
    typedef vector<bool> ccsf_b;

/*!
 * \brief Data structure for cell-centered vector field double data (not a 
 *         CAR_CU_Mesh nested mesh field class).
 */
    typedef vector< vector<double> > ccvf_d;
/*!
 * \brief Data structure for cell-centered vector field integer data (not a 
 *        CAR_CU_Mesh nested mesh field class).
 */
    typedef vector< vector<int> > ccvf_i;

/*!
 * \brief Data structure for face-centered scalar field double data (not a 
 *        CAR_CU_Mesh nested mesh field class).
 */
    typedef vector<double> fcsf_d;
/*!
 * \brief Data structure for face-centered vector field integer data (not a 
 *        CAR_CU_Mesh nested mesh field class).
 */
    typedef vector<int> fcsf_i;

/*!
 * \brief Data structure for face-centered discontinuous scalar field double 
 *        data (not a CAR_CU_Mesh nested mesh field class).
 */
    typedef vector<double> fcdsf_d;
/*!
 * \brief Data structure for face-centered discontinuous scalar field integer 
 *        data (not a CAR_CU_Mesh nested mesh field class).
 */
    typedef vector<int> fcdsf_i;

/*!
 * \brief Data structure for node-centered scalar field double data (not a 
 *        CAR_CU_Mesh nested mesh field class).
 */
    typedef vector<double> ncsf_d;
/*!
 * \brief Data structure for node-centered scalar field integer data (not a 
 *        CAR_CU_Mesh nested mesh field class).
 */
    typedef vector<int> ncsf_i;

/*!
 * \brief Data structure for node-centered vector field double data (not a 
 *        CAR_CU_Mesh nested mesh field class).
 */
    typedef vector< vector<double> > ncvf_d;
/*!
 * \brief Data structure for node-centered vector field integer data (not a 
 *        CAR_CU_Mesh nested mesh field class).
 */
    typedef vector< vector<int> > ncvf_i;

/*!
 * \brief Data structure for localized node-centered vector field double 
 *        data (not a CAR_CU_Mesh nested mesh field class, and not yet 
 *        implemented).
 */
    typedef multimap<int, vector<double> > lncvf_d;

    // temporary typedefs for compiling code until KCC 3.3+ is released
    // (retained herein for compatablity with historical codes).
/*!
 * \brief CAR_CU_Mesh nested cell-centered scalar field (CCSF) class with 
 *        double data.
 */
    typedef CCSF<double> CCSF_double;
/*!
 * \brief CAR_CU_Mesh nested cell-centered scalar mesh field (CCSF) class with 
 *        integer data.
 */
    typedef CCSF<int> CCSF_int;
/*!
 * \brief CAR_CU_Mesh nested cell-centered scalar mesh field (CCSF) class with 
 *        boolean data.
 */
    typedef CCSF<bool> CCSF_bool;
/*!
 * \brief CAR_CU_Mesh nested cell-centered scalar mesh field (CCSF) class with 
 *        string data.
 */
    typedef CCSF<string> CCSF_string;
/*!
 * \brief CAR_CU_Mesh nested cell-centered vector mesh field (CCVF) class with 
 *        double data.
 */
    typedef CCVF<double> CCVF_double;
/*!
 * \brief CAR_CU_Mesh nested cell-centered vector mesh field (CCVF) class with 
 *        integer data.
 */
    typedef CCVF<int> CCVF_int;

  private:
    // layout of mesh
/*!
 * \brief The mesh Layout class object.
 */
    Layout layout;
    // vertices in mesh
/*!
 * \brief The node coordinate values.
 */
    ncvf_d node_coords;
    // cell-pairings of cell to its vertices
/*!
 * \brief The connections between cell faces.
 */
    ccvf_i cell_nodes;
    // indicator whether this is a submesh
/*!
 * \brief An indicator that this is a submesh for parallel processing.
 */
    bool submesh;
    // indicator that the cell has been refined
/*!
 * \brief An indicator that the cells have been refined (not used).
 */
    ccsf_b has_kids;
    // indicator for the level of cell refinement, with the original coarse
    // mesh input by the user assigned as zero. 
/*!
 * \brief The cell generation (i.e., refinement) levels.
 */
    ccsf_i generation;

    // private functions

    // compare real values for equality
/*!
 * \brief Compares real values for equality considering machine precision.
 * \param low_val Values that are assumed to be at the lower value.
 * \param high_val Values that are assumed to be at the higher value.
 * \return high_val > low_val && 
 *         abs(low_val- high_val) > epsilon (desired precision).
 */
    bool compReal(const double & low_val, const double & high_val) const;

    // private copy and assignment operators; can't copy or assign a mesh
/*!
 * \brief Private copy operator; can't copy a mesh.
 * \param mesh CAR_CU_Mesh class object.
 */
    CAR_CU_Mesh(const CAR_CU_Mesh & mesh);
/*!
 * \brief Private assignment operator; can't assign a mesh.
 * \param mesh CAR_CU_Mesh class object.
 */
    CAR_CU_Mesh & operator=(const CAR_CU_Mesh & mesh);

    // Begin_Doc os_mesh-int.tex
    // Begin_Verbatim 

  public:
    // generalized constructor for all mesh types
    CAR_CU_Mesh(Layout & layout_, ncvf_d & node_coords_, ccvf_i & cell_nodes_,
		ccsf_i & generation_, bool submesh_ = false); 

    // mesh dimensionality functions
    // The problem geometry dimension (infered from the number of coordinate
    // directions).
/*!
 * \brief Returns the number of spatial dimensions for the mesh.
 * \return Number of spatial dimensions for the mesh. 
 */
    int get_ndim() const { return node_coords[0].size(); }

    // return number of cells
/*!
 * \brief Returns the number of cells in the mesh.
 * \return Number of mesh cells. 
 */
    int get_num_cells() const { return layout.get_num_cells(); }

    // return total number of nodes
/*!
 * \brief Returns the total number of nodes (i.e., corner plus face-centered) 
 *        in the mesh.
 * \return Total number of nodes in the mesh. 
 */
    int get_num_nodes() const { return node_coords.size(); }

    // return number of cell-corner nodes
/*!
 * \brief Returns the number of cell corner nodes in the mesh.
 * \return Number of cell corner nodes in the mesh. 
 */
    int get_num_corner_nodes() const 
    { return cell_nodes[get_num_cells() - 1][static_cast<int>(pow(2.0, 
					         get_ndim())) - 1]; }

    // return number of face-centered nodes
/*!
 * \brief Returns the number of cell face-centered nodes in the mesh.
 * \return Number of cell face-centered nodes in the mesh. 
 */
    int get_num_face_nodes() const 
    { return get_num_nodes() - get_num_corner_nodes(); }

    // cell dimensionality functions
    // give the dimension and begin and end return the beginning and ending
    // coordinate along that dimension
/*!
 * \brief Returns the minimum coordinate value along the specified direction
 *        for the entire mesh.
 * \param dir Coordinate direction (x=1, y=2, z =3).
 * \return Mimimum coordinate value. 
 */
    inline double get_mesh_min_coord(int dir) const;
/*!
 * \brief Returns the maximum coordinate value along the specified direction
 *        for the entire mesh.
 * \param dir Coordinate direction (x=1, y=2, z =3).
 * \return Maximum coordinate value. 
 */
    inline double get_mesh_max_coord(int dir) const;

    // find minimum and maximum dimension of cell
/*!
 * \brief Returns the minimum coordinate value along the specified direction
 *        for the specified cell.
 * \param dir Coordinate direction (x=1, y=2, z =3).
 * \param cell Cell number.
 * \return Mimimum coordinate value. 
 */
    inline double get_cell_min_coord(int dir, int cell) const;
/*!
 * \brief Returns the maximum coordinate value along the specified direction
 *        for the specified cell.
 * \param dir Coordinate direction (x=1, y=2, z =3).
 * \param cell Cell number.
 * \return Maximum coordinate value. 
 */
    inline double get_cell_max_coord(int dir, int cell) const;

    // find centerpoint of cell and width of cell
/*!
 * \brief Returns the center-point coordinate value along the specified
 *        direction for the specified cell.
 * \param dir Coordinate direction (x=1, y=2, z =3).
 * \param cell Cell number.
 * \return Center-point coordinate value. 
 */
    inline double get_cell_center_coord(int dir, int cell) const;
/*!
 * \brief Returns the cell width along the specified direction for the 
 *        specified cell.
 * \param dir Coordinate direction (x=1, y=2, z =3).
 * \param cell Cell number.
 * \return Cell width. 
 */
    double get_cell_width(int dir, int cell) const 
    { return get_cell_max_coord(dir, cell) - get_cell_min_coord(dir, cell);} 

    // diagnostic functions
    void print(ostream & output) const;
    void print(ostream & output, int cell) const;

    // End_Verbatim 
    // End_Doc 

    // Begin_Doc car_os_mesh-rint.tex
    // Begin_Verbatim 

    // member functions used by the CAR_CU_Mesh-dependent classes
    // services required by ALL mesh types used in JAYENNE
    // references to imbedded objects and data required for Parallel_Building
/*!
 * \brief Returns a reference to the mesh Layout class object.
 * \return Reference to the mesh Layout class object. 
 */
    const Layout & get_Layout() const { return layout; }
/*!
 * \brief Returns a reference to the mesh node_coords ncvf_d type object.
 * \return Reference to the mesh node_coords ncvf_d type object. 
 */
    const ncvf_d & get_node_coords() const { return node_coords; }
/*!
 * \brief Returns a reference to the mesh cell_nodes ccvf_i type object.
 * \return Reference to the mesh cell_nodes ccvf_i type object. 
 */
    const ccvf_i & get_cell_nodes() const { return cell_nodes; }
/*!
 * \brief Returns a reference to the mesh generation ccsf_i type object.
 * \return Reference to the mesh generation ccsf_i type object. 
 */
    const ccsf_i & get_generation() const { return generation; }

    // required services for transport; 
/*!
 * \brief Returns the number of cells adjacent to the specified cell face.
 * \param cell_index Cell number.
 * \param face_index Face number.
 * \return Number of adjacent cells. 
 */
    int get_num_adj_cells(int cell_index, int face_index) const
    { return layout.get_num_adj_cells(cell_index, face_index);}
/*!
 * \brief Returns the cell number of the specified adjacent cell for the 
 *        specified cell face.
 * \param cell Cell number.
 * \param face Face number.
 * \param adjcell Adjacent cell (defaults to the first adjacent cell).
 * \return Adjacent cell number. 
 */
    int get_next_cell(int cell, int face, int adjcell = 1) const 
    { return layout(cell, face, adjcell); }
    int get_cell_from_coords(const vector<double> & r) const;
    double get_dist_2_bndry(const vector<double> & r, 
			    const vector<double> & omega, int cell, 
			    int & face) const;
/*!
 * \brief Returns the node number of the specified cell node.
 * \param cell Cell number.
 * \param node Cell node index.
 * \return Node number. 
 */
    int get_cell_node(int cell, int node) const
    { return cell_nodes[cell][node];}
/*!
 * \brief Returns the node number for the specified cell face.
 * \param cell Cell number.
 * \param face Cell face index.
 * \return Node number. 
 */
    int get_cell_face_centered_node(int cell, int face) const
    {
        int offset = static_cast<int>(pow(2.0,get_ndim()));
	return cell_nodes[cell][offset + face];
    }
/*!
 * \brief Returns all of the nodes for the specified cell (i.e, both the 
 *        corner nodes and the face-centered nodes).
 * \param cell Cell number.
 * \return Cell nodes. 
 */
    inline vector<int> get_cell_nodes(int cell) const;
/*!
 * \brief Returns all of the corner nodes for the specified cell.
 * \param cell Cell number.
 * \return Cell corner nodes. 
 */
    inline vector<int> get_cell_corner_nodes(int cell) const;
/*!
 * \brief Returns all of the face-centered nodes for the specified cell.
 * \param cell Cell number.
 * \return Cell face-centered nodes. 
 */
    inline vector<int> get_cell_face_centered_nodes(int cell) const;
/*!
 * \brief Returns all of the corner nodes for the specified cell face.
 * \param cell Cell number.
 * \param face Face number.
 * \return Node number.
 */
    inline vector<int> get_cell_face_nodes(int cell, int face) const;
/*!
 * \brief Returns the normalized components of the outward-directed normal 
 *        for the specified cell face.
 * \param cell Cell number.
 * \param face Face number.
 * \return Outward-directed normal.
 */
    inline vector<double> get_cell_face_normal(int cell, int face) const;
/*!
 * \brief Returns the normalized components of the inward-directed normal 
 *        for the specified cell face.
 * \param cell Cell number.
 * \param face Face number.
 * \return Inward-directed normal.
 */
    inline vector<double> get_cell_face_normal_in(int cell, int face) const;
/*!
 * \brief Returns the volume of the specified cell.
 * \param cell Cell number.
 * \return Cell volume. 
 */
    inline double get_cell_volume(int cell) const;
/*!
 * \brief Returns the area of the specified cell face.
 * \param cell Cell number.
 * \param face Face number.
 * \return Cell face area. 
 */
    inline double get_cell_face_area(int cell, int face) const;
    vector<int> get_surcells(string boundary) const;
    void check_defined_surcells(const string, const vector<int> &) const;
    int get_bndface(string, int) const;
/*!
 * \brief Returns the coordinate values for the specified node.
 * \param node Node number.
 * \return Node coordinate values. 
 */
    inline vector<double> get_node_coords(int node) const;
/*!
 * \brief Returns all of the corner node coordinate values for the specified 
 *        cell face.
 * \param cell Cell number.
 * \param face Face number.
 * \return Coordinate values for the cell face nodes. 
 */
    inline ncvf_d get_cell_nodes_coords(int cell, int face) const;
/*!
 * \brief Returns all of the corner node coordinate values for the specified 
 *        cell.
 * \param cell Cell number.
 * \return Coordinate values for the cell corner nodes. 
 */
    inline ncvf_d get_cell_nodes_coords(int cell) const;
/*!
 * \brief Returns the generation (i.e., refinement) level for the specified 
 *        cell.
 * \param cell Cell number.
 * \return Cell generation level. 
 */
    inline int get_generation(int cell) const;
/*!
 * \brief Determines if a spatial position lies on the specified cell face.
 * \param pos Spatial position.
 * \param cell Cell number.
 * \param face Face number.
 * \return Status of the spatial position relative to the cell face. 
 */
    inline bool check_on_face(vector<double> & pos, int & cell, 
			      int & face) const; 

    // overloaded operators
/*!
 * \brief Compares a mesh for equivalence relative to the current mesh.
 * \param rhs Second mesh
 * \return Status of second mesh = current mesh. 
 */
    bool operator==(const CAR_CU_Mesh & rhs) const;
/*!
 * \brief Compares a mesh for non-equivalence relative to the current mesh.
 * \param rhs Second mesh
 * \return Status of second mesh != current mesh. 
 */
    bool operator!=(const CAR_CU_Mesh & rhs) const { return !(*this == rhs); }

    // End_Verbatim 
    // End_Doc 
};

//---------------------------------------------------------------------------//
// overloaded operators
//---------------------------------------------------------------------------//

/*!
 * \brief Overloaded stream-insertion operator for mesh output.
 * \param output Stream-output class object.
 * \param object CAR_CU_Mesh class object.
 * \return Reference to output.
 */
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
    : mesh(mesh_), data(mesh->get_num_cells()) 
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
    Ensure (data.size() == mesh->get_num_cells());
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
	data[i].resize(mesh->get_num_cells());
}

template<class T>
inline CAR_CU_Mesh::CCVF<T>::CCVF(SP<CAR_CU_Mesh> mesh_, int size_) 
    : mesh(mesh_), data(size_)
{
    Require (mesh);

    // initialize data array
    for (int i = 0; i <  data.size(); i++)
	data[i].resize(mesh->get_num_cells());
}

//---------------------------------------------------------------------------//
//---------------------------------------------------------------------------//
// constructor for automatic initialization
// CCVF explicit constructor (arbitrary vector leading index size)

template<class T>
inline CAR_CU_Mesh::CCVF<T>::CCVF(SP<CAR_CU_Mesh> mesh_,
				  const vector<vector<T> > & array)
    : mesh(mesh_), data(array)
{
    Require (mesh);

    for (int s = 0; s < data.size(); s++)
	Ensure (data[s].size() == mesh->get_num_cells());
}

//---------------------------------------------------------------------------//
// constant overloaded ()

template<class T>
inline const T & CAR_CU_Mesh::CCVF<T>::operator()(int dim, int cell) const 
{
    return data[dim][cell]; 
}

//---------------------------------------------------------------------------//
// assignment overloaded ()

template<class T>
inline T & CAR_CU_Mesh::CCVF<T>::operator()(int dim, int cell)
{
    return data[dim][cell];
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
	x.push_back(data[i][cell]);

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
    : mesh(mesh_), data(mesh->get_num_face_nodes())
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
    Ensure (data.size() == mesh->get_num_face_nodes());
}

//---------------------------------------------------------------------------//
// constant overloaded ()

template<class T>
inline const T & CAR_CU_Mesh::FCSF<T>::operator()(int face) const 
{
    return data[face]; 
}

//---------------------------------------------------------------------------//
// constant overloaded ()

template<class T>
inline const T & CAR_CU_Mesh::FCSF<T>::operator()(int cell, int face) const 
{
    int index = mesh->get_cell_node(cell, face + static_cast<int>(pow(2.0, 
        get_ndim())));
    return data[index - mesh->get_num_corner_nodes()]; 
}

//---------------------------------------------------------------------------//
// assignment overloaded ()

template<class T>
inline T & CAR_CU_Mesh::FCSF<T>::operator()(int face)
{
    return data[face];
}

//---------------------------------------------------------------------------//
// assignment overloaded ()

template<class T>
inline T & CAR_CU_Mesh::FCSF<T>::operator()(int cell, int face)
{
    int index = mesh->get_cell_node(cell, face + static_cast<int>(pow(2.0, 
        mesh->get_ndim())));
    return data[index - mesh->get_num_corner_nodes()];
}

//---------------------------------------------------------------------------//
// CAR_CU_Mesh::FCDSF inline functions
//---------------------------------------------------------------------------//
// FCDSF explicit constructor

template<class T>
inline CAR_CU_Mesh::FCDSF<T>::FCDSF(SP<CAR_CU_Mesh> mesh_)
    : mesh(mesh_), data(mesh->get_num_cells())
{
    Require (mesh);

    // initialize data array
    for (int i = 0; i <mesh->get_num_cells(); i++)
	data[i].resize(2.0 * mesh->get_ndim());
}

//---------------------------------------------------------------------------//
// constructor for automatic initialization

template<class T>
inline CAR_CU_Mesh::FCDSF<T>::FCDSF(SP<CAR_CU_Mesh> mesh_, 
			      const vector<vector<T> > & array)
    : mesh(mesh_), data(array)
{
    // check things out
    Ensure (data.size() == mesh->get_num_cells());
    for (int i = 0; i < mesh->get_num_cells(); i++)
	Ensure (data[i].size() == 2.0 * mesh->get_ndim());
}

//---------------------------------------------------------------------------//
// constant overloaded ()

template<class T>
inline const T & CAR_CU_Mesh::FCDSF<T>::operator()(int cell, int face) const 
{
    return data[cell][face]; 
}

//---------------------------------------------------------------------------//
// assignment overloaded ()

template<class T>
inline T & CAR_CU_Mesh::FCDSF<T>::operator()(int cell, int face)
{
    return data[cell][face];
}

//---------------------------------------------------------------------------//
// vector return overload()

template<class T>
inline vector<T> CAR_CU_Mesh::FCDSF<T>::operator()(int cell) const
{
    // declare return vector
    vector<T> x;
    
    // loop through faces and make return vector for this cell
    for (int i = 0; i < data[cell].size(); i++)
	x.push_back(data[cell][i]);

    // return
    Ensure (x.size() == data[cell].size());
    return x;
}

//---------------------------------------------------------------------------//
// CAR_CU_Mesh::FCVF inline functions
//---------------------------------------------------------------------------//
// FCVF explicit constructor with the size of the second index of the vector 
// defaulted to be equal to the problem geometry dimension.

template<class T>
inline CAR_CU_Mesh::FCVF<T>::FCVF(SP<CAR_CU_Mesh> mesh_)
    : mesh(mesh_), data(mesh->get_num_face_nodes())
{
    Require (mesh);

    // resize the second index of the vector to the default size equal to 
    // the problem geometry dimension
    for (int i = 0; i < data.size(); i++)
        data[i].resize(mesh->get_ndim());
}

// FCVF explicit constructor with the size of the second index of the vector 
// input.

template<class T>
inline CAR_CU_Mesh::FCVF<T>::FCVF(SP<CAR_CU_Mesh> mesh_, int vec_size)
    : mesh(mesh_), data(mesh->get_num_face_nodes())
{
    Require (mesh);

    // resize the second index of the vector to the input size
    for (int i = 0; i < data.size(); i++)
        data[i].resize(vec_size);
}

//---------------------------------------------------------------------------//
// constructor for automatic initialization

template<class T>
inline CAR_CU_Mesh::FCVF<T>::FCVF(SP<CAR_CU_Mesh> mesh_, 
			      const vector<vector<T> > & array)
    : mesh(mesh_), data(array)
{
    Require (mesh);
    // check things out
    Ensure (data.size() == mesh->get_num_face_nodes());
}

//---------------------------------------------------------------------------//
// constant overloaded ()

template<class T>
inline const vector<vector<T> > & CAR_CU_Mesh::FCVF<T>::operator()() const 
{
    return data;
}

//---------------------------------------------------------------------------//
// constant overloaded ()

template<class T>
inline const vector<T> & CAR_CU_Mesh::FCVF<T>::operator()(int cell, 
							  int face) const 
{
    // declare return vector
    vector<T> x;
    
    int index = mesh->get_cell_node(cell, face + static_cast<int>(pow(2.0, 
        get_ndim())));

    // loop through faces and make return vector for this cell
    for (int i = 0; i < data[index - mesh->get_num_corner_nodes()].size();
	 i++)
	x.push_back(data[index - mesh->get_num_corner_nodes()][i]);

    // return
    Ensure (x.size() == data[index - mesh->get_num_corner_nodes()].size());
    return x;
}

//---------------------------------------------------------------------------//
// constant overloaded ()

template<class T>
inline const T & CAR_CU_Mesh::FCVF<T>::operator()(int cell, int face, 
						  int dim) const 
{
    int index = mesh->get_cell_node(cell, face + static_cast<int>(pow(2.0, 
        get_ndim())));
    return data[index - mesh->get_num_corner_nodes()][dim]; 
}

//---------------------------------------------------------------------------//
// assignment overloaded ()

template<class T>
inline vector<vector<T> > & CAR_CU_Mesh::FCVF<T>::operator()()
{
    return data;
}

//---------------------------------------------------------------------------//
// assignment overloaded ()

template<class T>
inline vector<T> & CAR_CU_Mesh::FCVF<T>::operator()(int cell, int face)
{
    // declare return vector
    vector<T> x;
    
    int index = mesh->get_cell_node(cell, face + static_cast<int>(pow(2.0, 
        get_ndim())));

    // loop through faces and make return vector for this cell
    for (int i = 0; i < data[index - mesh->get_num_corner_nodes()].size();
	 i++)
	x.push_back(data[index - mesh->get_num_corner_nodes()][i]);

    // return
    Ensure (x.size() == data[index - mesh->get_num_corner_nodes()].size());
    return x;
}

//---------------------------------------------------------------------------//
// assignment overloaded ()

template<class T>
inline T & CAR_CU_Mesh::FCVF<T>::operator()(int cell, int face, int dim)
{
    int index = mesh->get_cell_node(cell, face + static_cast<int>(pow(2.0, 
        mesh->get_ndim())));
    return data[index - mesh->get_num_corner_nodes()][dim];
}

//---------------------------------------------------------------------------//
// CAR_CU_Mesh::NCSF inline functions
//---------------------------------------------------------------------------//
// NCSF explicit constructor (default vector size to the number of nodes)

template<class T>
inline CAR_CU_Mesh::NCSF<T>::NCSF(SP<CAR_CU_Mesh> mesh_) 
    : mesh(mesh_), data(mesh->get_num_nodes()) 
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

    Ensure (data.size() == mesh->get_num_nodes() ||
	    data.size() == mesh->get_num_corner_nodes())
}

//---------------------------------------------------------------------------//
// constructor for automatic initialization

template<class T>
inline CAR_CU_Mesh::NCSF<T>::NCSF(SP<CAR_CU_Mesh> mesh_, 
   const vector<T> & array) : mesh(mesh_), data(array)
{
    Require (mesh);

    Ensure (data.size() == mesh->get_num_nodes() ||
	    data.size() == mesh->get_num_corner_nodes())
}

//---------------------------------------------------------------------------//
// CAR_CU_Mesh::NCVF inline functions
//---------------------------------------------------------------------------//
// NCVF explicit constructor (default vector sizes to the number of nodes by
// the number of geometry dimensions)

template<class T>
inline CAR_CU_Mesh::NCVF<T>::NCVF(SP<CAR_CU_Mesh> mesh_) 
    : mesh(mesh_), data(mesh->get_num_nodes())
{
    Require (mesh);

    // initialize data array
    for (int i = 0; i < data.size(); i++)
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

    Ensure (size_1 == mesh->get_num_nodes() ||
	    size_1 == mesh->get_num_corner_nodes())
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

    Ensure (size_1 == mesh->get_num_nodes() ||
	    size_1 == mesh->get_num_corner_nodes())

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

    Ensure (data.size() == mesh->get_num_nodes() ||
	    data.size() == mesh->get_num_corner_nodes())
}

//---------------------------------------------------------------------------//
// constant overloaded ()

template<class T>
inline const T & CAR_CU_Mesh::NCVF<T>::operator()(int node, int dim) const 
{
    return data[node][dim]; 
}

//---------------------------------------------------------------------------//
// assignment overloaded ()

template<class T>
inline T & CAR_CU_Mesh::NCVF<T>::operator()(int node, int dim)
{
    return data[node][dim];
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
	x.push_back(data[node][i]);

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
    Ensure (data.size() == mesh->get_num_nodes());
    int dim = get_ndim();
    for (multimap<K, vector<T> >::const_iterator niter = data.begin();
	 niter != data.end(); niter++)
	Ensure (niter->second.size() == dim);
}

//---------------------------------------------------------------------------//
// constant overloaded ()

template<class K, class T>
inline const T & CAR_CU_Mesh::LNCVF<K,T>::operator()(int dim, int node) const 
{
    vector<T> local_vector = * data.find(node);
    return local_vector[dim]; 
}

//---------------------------------------------------------------------------//
// assignment overloaded ()

template<class K, class T>
inline T & CAR_CU_Mesh::LNCVF<K,T>::operator()(int dim, int node)
{
    vector<T> local_vector = * data.find(node);
    return local_vector[dim]; 
}

//---------------------------------------------------------------------------//
// vector return overload()

template<class K, class T>
inline vector<T> CAR_CU_Mesh::LNCVF<K,T>::operator()(int node) const
{
    return data.find(node)->second;
}

//---------------------------------------------------------------------------//
// CAR_CU_Mesh inline functions
//---------------------------------------------------------------------------//

inline double CAR_CU_Mesh::get_mesh_min_coord(int d) const 
{
    // find the minimum surface for d over the whole mesh
    return node_coords[0][d]; 
}

//---------------------------------------------------------------------------//

inline double CAR_CU_Mesh::get_mesh_max_coord(int d) const 
{
    // find the maximum surface for d over the whole mesh
    return node_coords[get_num_cells() - 1][d]; 
}

//---------------------------------------------------------------------------//

inline double CAR_CU_Mesh::get_cell_center_coord(int d, int cell) const
{
    // find center position of cell

    // set return value
    double return_pos = 0.0;

    // loop over all vertices and take average value to get the center
    // point 
    for (int i = 0; i < pow(2.0,get_ndim()); i++)
	return_pos += node_coords[cell_nodes[cell][i]][d];

    // return value
    return return_pos /  pow(2.0,get_ndim());     
}

//---------------------------------------------------------------------------//

inline double CAR_CU_Mesh::get_cell_min_coord(int d, int cell) const 
{	
    // find minimum dimension along d of cell

    // loop over all vertices and find the minimum
    double minimum = node_coords[cell_nodes[cell][0]][d];
    for (int i = 1; i < pow(2.0,get_ndim()); i++)
    {
	double point = node_coords[cell_nodes[cell][i]][d];

	// update the minimum value point
	if (point < minimum)
	    minimum = point;
    }
	
    // return minimum dimension
    return minimum;
}

//---------------------------------------------------------------------------//

inline double CAR_CU_Mesh::get_cell_max_coord(int d, int cell) const
{
    // find maximum dimension of cell

    // loop over all vertices and find the maximum
    double maximum = node_coords[cell_nodes[cell][0]][d];
    for (int i = 1; i < pow(2.0,get_ndim()); i++)
    {
	double point = node_coords[cell_nodes[cell][i]][d];

	// update the maximum value point
	if (point > maximum)
	    maximum = point;
    }

    // return maximum dimension
    return maximum;
}

//---------------------------------------------------------------------------//

inline double CAR_CU_Mesh::get_cell_volume(int cell) const 
{
    // calculate volume of cell

    // loop through dimensions and get volume
    double volume = 1.0;
    for (int d = 1; d <= get_ndim(); d++)
	volume *= get_cell_width(d, cell);

    // return volume
    return volume;
}

//---------------------------------------------------------------------------//

inline double CAR_CU_Mesh::get_cell_face_area(int cell, int face) const 
{
    // calculate area of face on cell

    // loop through dimensions and multiply off-dimension widths
    double face_area = 1.0;
    int dim_face_on = abs(face - get_ndim()) + 1 - 
                      face/(get_ndim() + 1);

    for (int d = 1; d <= get_ndim(); d++)
    {
	if (d != dim_face_on)
	    face_area *= get_cell_width(d, cell);
    }

    // return face_area
    return face_area;
}

//---------------------------------------------------------------------------//

inline vector<int> CAR_CU_Mesh::get_cell_nodes(int cell) const
{
    // Return the set of nodes that make up a cell, including both the corner
    // nodes and the face-centered nodes.
    int nnodes = pow(2.0,get_ndim()) + 2.0 * get_ndim();
    vector<int> node_set(nnodes);
    for (int node = 0; node < nnodes; node++)
        node_set[node] = cell_nodes[cell][node];

    return node_set;
}

//---------------------------------------------------------------------------//

inline vector<int> CAR_CU_Mesh::get_cell_corner_nodes(int cell) const
{
    // Return the set of corner nodes that make up a cell.
    int nnodes = pow(2.0,get_ndim());
    vector<int> node_set(nnodes);
    for (int node = 0; node < nnodes; node++)
        node_set[node] = cell_nodes[cell][node];

    return node_set;
}

//---------------------------------------------------------------------------//

inline vector<int> CAR_CU_Mesh::get_cell_face_centered_nodes(int cell) const
{
    // Return the set of face-centered nodes for a cell.
    int nnodes = pow(2.0,get_ndim()) + 2.0 * get_ndim();
    vector<int> node_set(2.0 * get_ndim());
    for (int node = pow(2.0,get_ndim()); node < nnodes; node++)
        node_set[node - pow(2.0,get_ndim())] = 
	    cell_nodes[cell][node];

    return node_set;
}

//---------------------------------------------------------------------------//

inline vector<int> CAR_CU_Mesh::get_cell_face_nodes(int cell, int face) const
{
    // Return the set of nodes that make up a cell face
    vector<int> node_set(pow(2.0,get_ndim() - 1));

    int node = face/3 + face/4 + 2 * (face/5);
    node_set[0] = cell_nodes[cell][node];

    node = ((face+1)%2) * (1 + face/2) + (face%2) * (4 + face/2);
    node_set[1] = cell_nodes[cell][node];

    if (get_ndim() == 3)
    {
        node = ((face+1)%2) * (3 * (1 + face/2) - 2 * (face/4)) +
	       (face%2) * (5 + 2 * (face/3));
        node_set[2] = cell_nodes[cell][node];

        node = ((face+1)%2) * (face + 2) + (face%2) * face;
        node_set[3] = cell_nodes[cell][node];
    }

    return node_set;
}

//---------------------------------------------------------------------------//

inline vector<double> CAR_CU_Mesh::get_cell_face_normal(int cell, int face) 
    const
{
    // normal always has 3 components
    vector<double> normal(get_ndim(), 0.0);
	
    // calculate normal based on face, (-z, -y, -x, +x, +y, +z), only
    // one coordinate is non-zero    
    int index = abs(face + 1 - get_ndim()) - face/get_ndim();
    if (face < get_ndim())
        normal[index] = -1.0;
    else
        normal[index] = 1.0;

    // return the normal
    return normal;
}

//---------------------------------------------------------------------------//

inline vector<double> CAR_CU_Mesh::get_cell_face_normal_in(int cell, int face)
    const
{
    // normal always has 3 components
    vector<double> normal(get_ndim(), 0.0);
	
    // calculate normal based on face, (-z, -y, -x, +x, +y, +z), only
    // one coordinate is non-zero    
    int index = abs(face + 1 - get_ndim()) - face/get_ndim();
    if (face >= get_ndim())
        normal[index] = -1.0;
    else
        normal[index] = 1.0;

    // return the normal
    return normal;
}

//---------------------------------------------------------------------------//
// Return the node vertex

inline vector<double> CAR_CU_Mesh::get_node_coords(int node) const
{
    // determine the vertices along a cell-face

    // return vertices
    vector<double> ret_vert(get_ndim());

    for (int d = 0; d < get_ndim(); d++)
        ret_vert[d] = node_coords[node][d];

    // return vector of vertices
    return ret_vert;
}

//---------------------------------------------------------------------------//
// calculate the vertices bounding a cell face

inline CAR_CU_Mesh::ncvf_d CAR_CU_Mesh::get_cell_nodes_coords(int cell, 
							      int face) const
{
    // determine the vertices along a cell-face

    // return vertices
    ncvf_d ret_vert(get_ndim());

    // determine axis dimension of surface (x=0, y=1, z=2)
    int axis = abs(face + 1 - get_ndim()) - face/get_ndim();
    double plane;
    if (face < get_ndim())
	plane = get_cell_min_coord(axis, cell);
    else
	plane = get_cell_max_coord(axis, cell);

    // loop over vertices in cell and get the vertices that are in the plane
    for (int i = 0; i < pow(2.0,get_ndim()); i++)
	if (plane == node_coords[cell_nodes[cell][i]][axis])
	    for (int d = 0; d < get_ndim(); d++)
		ret_vert[d].push_back(node_coords[cell_nodes[cell][i]][d]);

    // return vector of vertices
    return ret_vert;
}

//---------------------------------------------------------------------------//
// calculate the vertices bounding a cell

inline CAR_CU_Mesh::ncvf_d CAR_CU_Mesh::get_cell_nodes_coords(int cell) const
{
    // determine the vertices bounding a cell
    
    // return vertices
    ncvf_d ret_vert(get_ndim());

    // loop over cell vertices and build the cell vertices
    for (int i = 0; i < pow(2.0,get_ndim()); i++)
	for (int d = 0; d < get_ndim(); d++)
	    ret_vert[d].push_back(node_coords[cell_nodes[cell][i]][d]);

    // return vector of vertices
    return ret_vert;
}

//---------------------------------------------------------------------------//
// determine the generation level of a cell

inline int CAR_CU_Mesh::get_generation(int cell) const
{
    return generation[cell];
}

//---------------------------------------------------------------------------//
// check if a point lies on a face

inline bool CAR_CU_Mesh::check_on_face(vector<double> & pos, int & cell, 
				       int & face) const
{
    int off_face = abs(face + 1 - get_ndim()) - face/get_ndim();
    int on_face = 1;

    for (int d = 0; d < get_ndim(); d++)
    {
	if (d != off_face)
	{
	    double side_max = get_cell_max_coord(d, cell);
	    double side_min = get_cell_min_coord(d, cell);
	    // multiply integer by zero if this coordinate is not on the face
	    on_face *= (pos[d] >= side_min && pos[d] <= side_max);
	}
    }
    return (on_face != 0);
}

} // end namespace rtt_amr

#endif                          // __amr_CAR_CU_Mesh_hh__

//---------------------------------------------------------------------------//
//                              end of amr_mesh/Mesh.hh
//---------------------------------------------------------------------------//
