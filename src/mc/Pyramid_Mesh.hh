//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/Pyramid_Mesh.hh
 * \author Jeffery Densmore (Stolen from RZWedge_Mesh.hh)
 * \date   Mon Oct  6 09:15:12 2003
 * \brief  Pyramid_Mesh header file.
 * \note   Copyright © 2003 The Regents of the University of California.
 *
 *
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef rtt_mc_Pyramid_Mesh_hh
#define rtt_mc_Pyramid_Mesh_hh

#include "Coord_sys.hh"
#include "AMR_Layout.hh"
//#include "rng/Sprng.hh"
#include "ds++/SP.hh"
//#include "ds++/Assert.hh"
//#include "Math.hh"
//#include "Constants.hh"
//#include "Sampler.hh"
#include <vector>
//#include <iostream>
#include <string>
//#include <utility>

namespace rtt_mc
{

//===========================================================================//
/*!
 * \class Pyramid_Mesh

 * \brief An XYZ Pyramid mesh constructed from an R mesh.
 *
 * I'll put a description in later

 */

// revision history:
// -----------------
// 0) (Mon Oct  6 09:15:12 2003) Jeffery Densmore: original
// 
//===========================================================================//

class Pyramid_Mesh 
{

  public:
    // Forward declaration of pack class
    //    struct Pack;

    // Forward declarations of cell-centered fields
    //template<class T> class CCSF;
    //template<class T> class CCVF;

  public:
    // Typedefs used throughout Pyramid_Mesh class.
    // typedef rtt_dsxx::SP<Pyramid_Mesh> SP_Mesh;
    typedef rtt_dsxx::SP<Coord_sys>            SP_Coord;
    //    typedef rtt_dsxx::SP<Pyramid_Mesh::Pack>  SP_Pack;
    // typedef rtt_rng::Sprng rng_Sprng;
    typedef std::vector<int>                   sf_int;
    //    typedef std::vector<std::vector<int>>     vf_int;
    typedef std::vector<double>                sf_double;
    typedef std::vector<std::vector<double> >  vf_double;
    typedef std::string                        std_string;
    //    typedef std::pair<sf_double, sf_double>   pair_sf_double;

    // Handy typedefs to CC fields (not formally needed in KCC3.3+).
    // typedef CCSF<double> CCSF_double;
    // typedef CCSF<int> CCSF_int;
    // typedef CCVF<int> CCVF_double;
    // typedef CCVF<int> CCVF_int;
    // typedef CCSF<std_string> CCSF_string;
  
  private:
    // Base class reference to a derived coordinate system class
    // (the Pyramid_Mesh is always three-dimensional, Cartesian)
    SP_Coord coord;

    // Layout (cell connectivity of mesh
    AMR_Layout layout;
    
    // vector<vector> of x-extents of each cell
    vf_double cell_x_extents;

    // Pyramid angle data; precalculated
    double beta_radians;
    double tan_beta;
    double sin_beta;
    double cos_beta;

    // Total, processor-local volume
    double total_volume;

    // >>> Private implementations

    // Pack up the cell extents
    //    void pack_extents(const sf_int &, char *, int, int) const;
    
    // Function to calculate frequently used wedge data
    void calc_angle_data(const double beta_radians);

    // Function to calculate and set the total, on-processor volume
    void calc_total_volume();

    // Private copy assignment operators (can't copy or assign a mesh).
    Pyramid_Mesh(const Pyramid_Mesh &);
    Pyramid_Mesh& operator=(const Pyramid_Mesh &);
  
  public:
    // Constructor
    Pyramid_Mesh(SP_Coord coord_, AMR_Layout &layout_, 
		 vf_double &cell_x_extents_, double beta_radians_);

    // >>> Member functions used by Pyramid_Mesh dependent classes
    
    // Return the number of cells.
    int num_cells() const { return layout.num_cells(); }

    // Return the x-dimension cell extents
    double get_low_x(int cell) const {return cell_x_extents[cell-1][0]; }
    double get_high_x(int cell) const {return cell_x_extents[cell-1][1]; }

    // Get the midpoint of a cell for a give dimension
    inline double get_x_midpoint(int cell) const;
    inline double get_y_midpoint(int cell) const;
    inline double get_z_midpoint(int cell) const;

    // get the dimension of a cell for a given coordinate
    // (y and z-coord => dim at midpoint)
    inline double dim(const int coordinate, const int cell) const;


    // Determine if a position is in a cell
    bool in_cell(int cell, const sf_double & r) const;

    // Diagnostic functions.
    //    void print(std::ostream &) const;
    //    void print(std::ostream &, int) const;

    // References to embedded objects.
    // const AMR_Layout& get_Layout() const {return layout; }
    const Coord_sys & get_Coord() const {return *coord;}
    SP_Coord get_SPCoord() const { return coord;}

    // Access total, on-procceser Pyramid volume
    inline double get_total_volume() const;

    // Services required for graphics dumps.
    sf_int get_cell_types() const;
    vf_double get_point_coord() const;
    //    vf_int get_cell_pair() const;

    // Required services for transport and source
    // get_spatial_dimension() const {return 3;}
    inline int next_cell(int cell, int face) const;
    int get_cell(const sf_double &) const;
    //    double get_db(const sf_double &, const sf_double &, int, int &)
    //    const;
    // inline double get_random_walk_sphere_radius(const sf_double &,int) const;
    inline sf_double get_normal(int cell, int face) const;
    inline sf_double get_normal_in(int cell, int face) const;
    inline double volume(int cell) const;
    inline double face_area(int cell,int face) const;
    vf_double get_vertices(int cell) const;
    vf_double get_vertices(int cell, int face) const;
    // inline sf_double sample_pos(int, rng_Sprng &) const;
    // inline sf_double sample_pos(int, rng_Sprng &, sf_double, double)
    // const;
    // inline sf_double sample_pos_on_face(int,int, rng_Sprng &) const;
    int get_bndface(std_string boundary, int cell) const;
    sf_int get_surcells(std_string boundary) const;
    bool check_defined_surcells(const std_string ss_face, 
				const sf_int &ss_list) const;
    inline sf_int get_neighbors(int cell) const;
    //    pair_sf_double sample_random_walk_sphere(int, const sf_double &, 
    //					     double, rng_Sprng &) const;

    //! Determine if this is a full mesh or partioned mesh (always full).
    bool full_Mesh() const { return 1;}

    // Pack function.
//    SP_Pack pack(const sf_int & = sf_int()) const;

    // Overloaded Operators.
//  bool operator==(const Pyramid_Mesh &) const;
    // bool operator!=(const Pyramid_Mes &rhs) const {return !(*this == rhs);
    // }



};
//---------------------------------------------------------------------------//
// PYRAMID_MESH INLINE FUNCTIONS
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
// OVERLOADED OPERATORS
//---------------------------------------------------------------------------//

//std::ostream & operator<<(std:ostream &output, const Pyramid_Mesh &object);

//---------------------------------------------------------------------------//
// PYRAMID_MESH SERVICES FOR PYRAMID_MESH-DEPENDENT CLASSES
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
/*! 
 * \brief Calculate the midpoint of a cell for the x-dimension.
 * 
 * \param cell Pyramid_Mesh cell.
 * \return midpoint of x-dimension
 */
double Pyramid_Mesh::get_x_midpoint(int cell) const
{
    Check ( ( cell >0) && (cell <=num_cells() ));
    return 0.5*(get_low_x(cell)+get_high_x(cell));
}
//---------------------------------------------------------------------------//
/*! 
 * \brief Calculate the midpoint of a cell for the y-dimension
 * 
 * \param cell Pyramid_Mesh cell
 *
 * \return midpoint of y-dimension
 */
double Pyramid_Mesh::get_y_midpoint(int cell) const
{
    Check ((cell>0) && (cell <= num_cells()));

    return 0.0;
}
//---------------------------------------------------------------------------//
/*! 
 * \brief Calculate the midpoint of a cell for the z-dimension
 * 
 * \param cell Pyramid_Mesh cell
 *
 * \return midpoint of z-dimension
 */
double Pyramid_Mesh::get_z_midpoint(int cell) const
{
    Check ((cell>0) && (cell <= num_cells()));

    return 0.0;
}
//---------------------------------------------------------------------------//
/*! 
 * \brief Get the dimension of the cell in the requested coordinate
 * 
 * \param cell Pryamid_Mesh cell.
 * \param coordinate coordinate
 *
 * \return dim the dimension (length) of the cell in the requested coordinate
 */
double Pyramid_Mesh::dim(const int coordinate, const int cell) const
{
    Check ( (cell>0) && (cell<= num_cells()) );
    Check ( (coordinate>0) && (coordinate <=3) );

    // return value
    double dimension = 0.0;

    // x-coordinate
    if (coordinate==1)
	dimension = get_high_x(cell)-get_low_x(cell);

    // y-coordinate
    else if (coordinate==2)
	dimension = 2.0 *get_x_midpoint(cell)*tan_beta;

    // z-coordinate
    else if (coordinate==3)
	dimension=2.0*get_x_midpoint(cell)*tan_beta;
    else
	Insist (0,"Requested coordinate in Pyramid_Mesh's dim no valid!");

    return dimension;
}


//---------------------------------------------------------------------------//
// RZWEDGE_MESH GENERALIZED MT SERVICES REQUIRED BY IMC
//---------------------------------------------------------------------------//

/*! 
 * \brief Calculate the cell across a face.
 * 
 *
 * \param cell current cell index
 * \param face face index
 * \return cell_across the cell across the face
 */
int Pyramid_Mesh::next_cell(int cell, int face) const
{
    Require (cell>0 && cell <= layout.num_cells());
    Require (face>0 && face<=6);

    // declare return cell
    int cell_across;

    if (layout.num_cells_across(cell,face) ==1)
    {
	cell_across = layout(cell,face,1);

	Check(face ==3 || face == 4 || face == 5 || face == 6 
	      ? cell_across == cell:cell_across==layout(cell,face,1));
    }
    else 
	Insist(0,"Must have exactly one cell across face in Pyramid Mesh!");

    return cell_across;
}


//---------------------------------------------------------------------------//
/*!
 * \brief Calculate the outward normal for the particular face of a cell.
 *
 * \param cell cell -- not used -- each Pyramid_Mesh cell has the same 
 * normals.
 * \param face face of cell for which to return outward normal
 *
 * \return outward normal
 */
Pyramid_Mesh::sf_double Pyramid_Mesh::get_normal(int cell, int face) const
{
    Check (coord->get_dim() == 3);
    Check ((face >=1) && (face<=6));

    sf_double normal(coord->get_dim(),0.0);

    // low x face
    if (face==1)
	normal[0]=-1.0;

    // high x face
    else if (face==2)
	normal[0]=1.0;

    // low y face
    else if (face==3)
    {
	normal[0]= -sin_beta;
	normal[1]= -cos_beta;
    }

    // high y face
    else if (face ==4)
    {
	normal[0]=-sin_beta;
	normal[1]=cos_beta;
    }
 

    // low z face
    else if (face==5)
    {
	normal[0]= -sin_beta;
	normal[2]= -cos_beta;
    }

    // high z face
    else if (face ==6)
    {
	normal[0]=-sin_beta;
	normal[2]=cos_beta;
    }
    
    // return outward normal
    return normal;
}
//---------------------------------------------------------------------------//
/*! 
 * \brief Calculate the inward normal for the particular face of a cell.
 * 
 * \param cell cell -- not used -- each Pyramid_Mesh cell has same normals
 * \param face face of cell for which to return inward normal.
 *
 * \return inward normal
 */
Pyramid_Mesh::sf_double Pyramid_Mesh::get_normal_in(int cell, int face) const
{
    Check (coord->get_dim() == 3);
    Check ((face>=1) && (face <= 6));

    // initialize inward normal
    sf_double normal_in(coord->get_dim(), 0.0);

    // get outward normal first
    sf_double normal = get_normal(cell, face);
    Check (normal.size() == coord->get_dim());
    
    // reverse direction
    for (int dir = 0; dir <coord->get_dim(); dir++)
	normal_in[dir]=-normal[dir];

    // return inward normal
    return normal_in;
}
//---------------------------------------------------------------------------//
/*!
 * \brief Return the volume of a Pyramid_Mesh cell.
 *
 * \param cell Pyramid_Mesh cell
 */
double Pyramid_Mesh::volume(int cell) const
{
    Require (cell>0 && cell <=num_cells());
    
    double lox = get_low_x(cell);
    double hix = get_high_x(cell);

    double vol = (4./3.)*tan_beta*tan_beta*(hix*hix*hix-lox*lox*lox);

    Ensure (vol>0.0);

    return vol;

}
//---------------------------------------------------------------------------//
/*! 
 * \brief Access the total (on-processor) volume of a Pyramid_Mesh
 * 
 * \return total volume of on-processor Pyramid cells
 */
double Pyramid_Mesh::get_total_volume() const
{
    Ensure (total_volume>0.0);

    return total_volume;
}
//---------------------------------------------------------------------------//
/*! 
 * \brief Return the area of a face in a Pyramid_Mesh cell.
 * 
 * \param cell Pyramid_Mesh cell
 * \param face face of cell.
 *
 * \return face area
 */
double Pyramid_Mesh::face_area(int cell, int face) const
{
    Require (face>0 && face <=6);
    Require (cell>0 && cell<=num_cells());

    // initialize face area
    double face_area =0.0;

    // low x face
    if(face ==1)
	face_area=4.0*get_low_x(cell)*get_low_x(cell)*tan_beta*tan_beta;

    // high x face
    if(face ==2)
	face_area=4.0*get_high_x(cell)*get_high_x(cell)*tan_beta*tan_beta;

    // all other faces
    else if (face ==3 || face ==4 || face == 5 || face ==6)
    {
	face_area=get_high_x(cell)*get_high_x(cell);
	face_area-=get_low_x(cell)*get_low_x(cell);
	face_area/=cos_beta;
    }

    // return face_area;
    return face_area;
}
//---------------------------------------------------------------------------//
/*! 
 * \brief Calculate the neighbors around a cell.
 * 
 * \param cell cell index
 * \return vector containing list of neighboring cells
 */
Pyramid_Mesh::sf_int Pyramid_Mesh::get_neighbors(int cell) const
{
    Require (layout.num_faces(cell) == 6);

    sf_int neighbors;

    for( int face =1; face<= layout.num_faces(cell); face++)
    {

	// loop over number of cells across the face and add the cells
	// across to the neighbor list
	Check(layout.num_cells_across(cell,face)==1);
	for (int i =1; i<= layout.num_cells_across(cell,face); i++)
	    neighbors.push_back(layout(cell,face,i));
    }
    Check(neighbors.size()==layout.num_faces(cell));

    return neighbors;
}

//===========================================================================//
/*!
 * \struct  Pyramid_Mesh::Pack
 * \brief  Pack and unpack a Pyramid_Mesh instance into raw c-style data
 * arrays.
 */
//===========================================================================//
//struct Pyramid_Mesh::Pack
//{
//  private:
    // Data contained in the mesh.
//    char *data;
//    int size;

    //disallow assignment.
//    const Pack& operator=(const Pack &);

//  public:
    // constructor
//    Pack(int, char *);

    // copy constructor.
//    Pack(const Pack &);

    // Destructor
//    ~Pack();

    // >>> Accessors

    //! Get pointer to beginning of char data stream.
//    const char* begin() const {return & data[0];}

    //! Get pointer to end of char data stream.
//    const char* end() const { return &data[size];}

    //! Return the number of cells in the packed mesh.
//    int get_num_packed_cells() const;

    //! Get size of data stream.
//    int get_size() const {return size; }

    // Unpack function.
//    SP_Mesh unpack() const;
//}


} // end namespace rtt_mc

#endif // rtt_mc_Pyramid_Mesh_hh

//---------------------------------------------------------------------------//
//              end of mc/Pyramid_Mesh.hh
//---------------------------------------------------------------------------//
