//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/Sphyramid_Mesh.hh
 * \author Jeffery Densmore (Stolen from RZWedge_Mesh.hh)
 * \date   Mon Nov 10 12:53 2003
 * \brief  Sphyramid_Mesh header file.
 * \note   Copyright © 2003 The Regents of the University of California.
 *
 *
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef rtt_mc_Sphyramid_Mesh_hh
#define rtt_mc_Sphyramid_Mesh_hh

#include "Coord_sys.hh"
#include "Layout.hh"
#include "rng/Sprng.hh"
#include "ds++/SP.hh"
#include "ds++/Assert.hh"
#include "Math.hh"
//#include "Constants.hh"
#include "Sampler.hh"
#include <vector>
//#include <iostream>
#include <string>
//#include <utility>

namespace rtt_mc
{

//===========================================================================//
/*!
 * \class Sphyramid_Mesh

 * \brief An XYZ Pyramid mesh constructed from an R Mesh
 *
 * We model an R mesh with an XYZ "Sphyramid Mesh".  On this mesh, all
 * cells are hexes, except at r=0, where the inner degenerate face is a
 * point.  In the construction of this mesh, a spherical cone (the original
 * mesh) is modeled with a pyramid, and the x-direction is aligned with the 
 * r-direction.  
 * The construction of this mesh requires an R mesh and an unfolding angle
 * alpha (the maximum polar angle of the spherical cone.  The Sphyramid mesh
 * then generates an x mesh and an unfolding angle beta (the angle between
 * the x-axis and sides of the mesh).  Please see memos CCS-4:03-56(U) and 
 * CCS-4:03-??(U)
 *
 * 

 */

// revision history:
// -----------------
// 0) (Mon Oct  6 09:15:12 2003) Jeffery Densmore: original
// 
//===========================================================================//

class Sphyramid_Mesh 
{

  public:
    // Forward declaration of pack class
    //    struct Pack;

    // Forward declarations of cell-centered fields
    //template<class T> class CCSF;
    //template<class T> class CCVF;

  public:
    // Typedefs used throughout Sphyramid_Mesh class.
    // typedef rtt_dsxx::SP<Sphyramid_Mesh> SP_Mesh;
    typedef rtt_dsxx::SP<Coord_sys>            SP_Coord;
    // typedef rtt_dsxx::SP<Sphyramid_Mesh::Pack>  SP_Pack;
    typedef rtt_rng::Sprng                     rng_Sprng;
    typedef std::vector<int>                   sf_int;
    // typedef std::vector<std::vector<int>>     vf_int;
    typedef std::vector<double>                sf_double;
    typedef std::vector<std::vector<double> >  vf_double;
    typedef std::string                        std_string;
    // typedef std::pair<sf_double, sf_double>   pair_sf_double;

    // Handy typedefs to CC fields (not formally needed in KCC3.3+).
    // typedef CCSF<double> CCSF_double;
    // typedef CCSF<int> CCSF_int;
    // typedef CCVF<int> CCVF_double;
    // typedef CCVF<int> CCVF_int;
    // typedef CCSF<std_string> CCSF_string;
  
  private:
    // Base class reference to a derived coordinate system class
    // (the Sphyramid_Mesh is always three-dimensional, Cartesian)
    SP_Coord coord;

    // Layout (cell connectivity) of mesh
    Layout layout;
    
    // vector<vector> of x-extents of each cell
    vf_double cell_x_extents;

    // Sphyramid angle data; precalculated
    double beta_radians;
    double tan_beta;
    double sin_beta;
    double cos_beta;

    // Total, processor-local volume
    double total_volume;

    // >>>> Private implementations <<<<

    // Pack up the cell extents
    // void pack_extents(const sf_int &, char *, int, int) const;
    
    // Function to calculate frequently used wedge data
    void calc_angle_data();

    // Function to calculate and set the total, on-processor volume
    void calc_total_volume();

    // Private copy assignment operators (can't copy or assign a mesh).
    Sphyramid_Mesh(const Sphyramid_Mesh &);
    Sphyramid_Mesh& operator=(const Sphyramid_Mesh &);
  
  public:
    // Constructor
    Sphyramid_Mesh(SP_Coord coord_, Layout &layout_, 
		   vf_double &cell_x_extents_, double beta_radians_);

    // Return the number of cells.
    int num_cells() const { return layout.num_cells(); }

    // Return the x-dimension cell extents
    double get_low_x(int cell) const { return cell_x_extents[cell-1][0]; }
    double get_high_x(int cell) const { return cell_x_extents[cell-1][1]; }

    // Get the midpoint of a cell for a given dimension
    inline double get_x_midpoint(int cell) const;
    inline double get_y_midpoint(int cell) const;
    inline double get_z_midpoint(int cell) const;

    // get the dimension of a cell for a given coordinate
    // (y and z-coord => dim at midpoint)
    inline double dim(const int coordinate, const int cell) const;


    // Determine if a position is in a cell
    bool in_cell(int cell, const sf_double & r) const;

    // Diagnostic functions
    // void print(std::ostream &) const;
    // void print(std::ostream &, int) const;

    // References to embedded objects
    // const AMR_Layout& get_Layout() const {return layout; }
    const Coord_sys & get_Coord() const { return *coord; }
    SP_Coord get_SPCoord() const { return coord; }

    // Access total, on-procceser Sphyramid volume
    inline double get_total_volume() const;

    // Services required for graphics dumps.
    sf_int get_cell_types() const;
    vf_double get_point_coord() const;
    // vf_int get_cell_pair() const;

    // Services for transport and source
    // get_spatial_dimension() const {return 3;}
    inline int next_cell(int cell, int face) const;
    int get_cell(const sf_double &) const;
    // double get_db(const sf_double &, const sf_double &, int, int &) const;
    // inline double get_random_walk_sphere_radius(const sf_double &,int) const;
    inline sf_double get_normal(int cell, int face) const;
    inline sf_double get_normal_in(int cell, int face) const;
    inline double volume(int cell) const;
    inline double face_area(int cell,int face) const;
    vf_double get_vertices(int cell) const;
    vf_double get_vertices(int cell, int face) const;
    inline sf_double sample_pos(int cell, rng_Sprng &random) const;
    inline sf_double sample_pos(int cell, rng_Sprng &random, sf_double slope, 
				double center_value) const;
    inline sf_double sample_pos_on_face(int cell, int face, rng_Sprng &random) const;
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
    // SP_Pack pack(const sf_int & = sf_int()) const;

    // Overloaded Operators.
    bool operator==(const Sphyramid_Mesh &rhs) const;
    bool operator!=(const Sphyramid_Mesh &rhs) const {return !(*this == rhs); }



};
//---------------------------------------------------------------------------//
// SPHYRAMID_MESH INLINE FUNCTIONS
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
/*! 
 * \brief Calculate the midpoint of a cell for the x-dimension.
 * 
 * \param cell cell number.
 * \return midpoint of x-dimension
 */
double Sphyramid_Mesh::get_x_midpoint(int cell) const
{
    Require (cell > 0);
    Require (cell <= num_cells());
    return 0.5*(get_low_x(cell)+get_high_x(cell));
}

//---------------------------------------------------------------------------//
/*! 
 * \brief Calculate the midpoint of a cell for the y-dimension
 * 
 * \param cell cell number
 *
 * \return midpoint of y-dimension
 */
double Sphyramid_Mesh::get_y_midpoint(int cell) const
{
    Require (cell > 0); 
    Require (cell <= num_cells());

    return 0.0;
}
//---------------------------------------------------------------------------//
/*! 
 * \brief Calculate the midpoint of a cell for the z-dimension
 * 
 * \param cell cell number
 *
 * \return midpoint of z-dimension
 */
double Sphyramid_Mesh::get_z_midpoint(int cell) const
{
    Require (cell > 0);
    Require (cell <= num_cells());

    return 0.0;
}
//---------------------------------------------------------------------------//
/*! 
 * \brief Get the dimension of the cell in the requested coordinate
 * 
 * \param cell cell number.
 * \param coordinate coordinate
 *
 * \return the dimension (length) of the cell in the requested coordinate
 */
double Sphyramid_Mesh::dim(const int coordinate, const int cell) const
{
    Require (cell > 0); 
    Require (cell <= num_cells() );
    Require (coordinate > 0);
    Require (coordinate <= 3);

    // return value
    double dimension = 0.0;

    // x-coordinate
    if (coordinate == 1)
    {
	dimension = get_high_x(cell)-get_low_x(cell);
    }

    // y-coordinate
    else if (coordinate==2)
    {
	dimension = 2.0 *get_x_midpoint(cell)*tan_beta;
    }

    // z-coordinate
    else if (coordinate==3)
    {
	dimension = 2.0*get_x_midpoint(cell)*tan_beta;
    }
    else
	Insist (0,"Requested coordinate in Sphyramid_Mesh's dim not valid!");

    return dimension;
}
//---------------------------------------------------------------------------//
/*! 
 * \brief Access the total (on-processor) volume of a Sphyramid_Mesh
 * 
 * \return total volume of on-processor Sphyramid cells
 */
double Sphyramid_Mesh::get_total_volume() const
{
    Ensure (total_volume > 0.0);

    return total_volume;
}
//---------------------------------------------------------------------------//
/*! 
 * \brief Calculate the cell across a face.
 * 
 *
 * \param cell current cell index
 * \param face face index
 * \return the cell across the face
 */
int Sphyramid_Mesh::next_cell(int cell, int face) const
{
    Require (cell > 0); 
    Require (cell <= this->layout.num_cells());
    Require (face > 0);
    Require (face <= 6);

    // declare return cell
    int cell_across;

    cell_across = this->layout(cell,face);
    
    Check (face == 3 || face == 4 || face == 5 || face == 6 
	      ? cell_across == cell : cell_across == layout(cell,face));

    return cell_across;
}
//---------------------------------------------------------------------------//
/*!
 * \brief Calculate the outward normal for the particular face of a cell.
 *
 * \param cell cell -- not used -- each Sphyramid_Mesh cell has the same 
 * normals.
 * \param face face of cell for which to return outward normal
 *
 * \return outward normal
 */
Sphyramid_Mesh::sf_double Sphyramid_Mesh::get_normal(int cell, int face) const
{
    Check (this->coord->get_dim() == 3);
    Require (face >= 1);
    Require (face <= 6);

    sf_double normal(this->coord->get_dim(),0.0);

    // low x face
    if (face == 1)
    {
	normal[0] = -1.0;
    }
    
    // high x face
    else if (face == 2)
    {
	normal[0] = 1.0;
    }

    // low y face
    else if (face == 3)
    {
	normal[0] = -this->sin_beta;
	normal[1] = -this->cos_beta;
    }

    // high y face
    else if (face == 4)
    {
	normal[0] = -this->sin_beta;
	normal[1] =  this->cos_beta;
    }

    // low z face
    else if (face == 5)
    {
	normal[0] = -this->sin_beta;
	normal[2] = -this->cos_beta;
    }

    // high z face
    else if (face == 6)
    {
	normal[0] = -this->sin_beta;
	normal[2] =  this->cos_beta;
    }
    
    // return outward normal
    return normal;
}
//---------------------------------------------------------------------------//
/*! 
 * \brief Calculate the inward normal for the particular face of a cell.
 * 
 * \param cell cell -- not used -- each Sphyramid_Mesh cell has same normals
 * \param face face of cell for which to return inward normal.
 *
 * \return inward normal
 */
Sphyramid_Mesh::sf_double Sphyramid_Mesh::get_normal_in(int cell, int face) 
    const
{
    Check (this->coord->get_dim() == 3);
    Require (face >= 1);
    Require (face <= 6);

    // initialize inward normal
    sf_double normal_in(this->coord->get_dim(), 0.0);

    // get outward normal first
    sf_double normal = get_normal(cell, face);
    Check (normal.size() == this->coord->get_dim());
    
    // reverse direction
    for (int dir = 0; dir < this->coord->get_dim(); dir++)
    {
	normal_in[dir] = -normal[dir];
    }

    // return inward normal
    return normal_in;
}
//---------------------------------------------------------------------------//
/*!
 * \brief Return the volume of a Sphyramid mesh cell.
 *
 * \param cell cell number
 */
double Sphyramid_Mesh::volume(int cell) const
{
    Require (cell > 0);
    Require (cell <= num_cells());
    
    double lox = get_low_x(cell);
    double hix = get_high_x(cell);

    Check (hix >= 0.);
    Check (lox >= 0.);

    double vol = (4./3.)*this->tan_beta*this->tan_beta*
	(hix*hix*hix-lox*lox*lox);

    Ensure (vol > 0.0);

    return vol;
}
//---------------------------------------------------------------------------//
/*! 
 * \brief Return the area of a face in a Sphyramid_Mesh cell.
 * 
 * \param cell cell number
 * \param face face of cell.
 *
 * \return face area
 */
double Sphyramid_Mesh::face_area(int cell, int face) const
{
    Require (face > 0);
    Require (face <= 6);
    Require (cell > 0);
    Require (cell <= num_cells());

    // initialize face_area and side
    double face_area;
    double side;

    // low x face
    if (face == 1)
    {
	side      = 2.0*get_low_x(cell)*this->tan_beta;
	face_area = side*side;
    }

    // high x face
    if (face == 2)
    {
	side      = 2.0*get_high_x(cell)*this->tan_beta;
	face_area = side*side;
    }

    // all other faces
    else if (face == 3 || face == 4 || face == 5 || face == 6)
    {
	face_area  = get_high_x(cell)*get_high_x(cell);
	face_area -= get_low_x(cell)*get_low_x(cell);
	face_area *= this->tan_beta/this->cos_beta;
    }

    // return face_area;
    Ensure (face_area >= 0.);
    return face_area;
}
//---------------------------------------------------------------------------//
/*! 
 * \brief Sample a position uniformly in an Sphyramid_Mesh cell
 * the form of the (unormalized) pdf is
 * f(x)=x^2 x_low<x<x_high
 * 
 * \param cell cell number
 * \param random random number object
 * \return location in cell
 */
Sphyramid_Mesh::sf_double Sphyramid_Mesh::sample_pos(int cell, 
						     rng_Sprng &random) const
{
    Require (cell > 0);
    Require (cell <= num_cells());
    Check   (this->coord->get_dim() == 3);

    // initialize location vector
    sf_double position(this->coord->get_dim(), 0.0);

    // get x-dimension cell extents
    double lox = get_low_x(cell);
    double hix = get_high_x(cell);
    Check (lox >= 0.0);
    Check (hix >= 0.0);

    // calculate corresponding y (and z) cell extents
    double loy = lox*this->tan_beta;
    double hiy = hix*this->tan_beta;
    Check (loy >= 0.0);
    Check (hiy >= 0.0);
    
    // sample x uniformly in cell
    position[0] = rtt_mc::sampler::sample_xsquared(random, lox, hix);
    
    // sample y and z value
    double pos_y = position[0]*this->tan_beta;
    double neg_y =-pos_y;
    Check (pos_y >= loy);
    Check (pos_y <= hiy);

    position[1] = neg_y+random.ran()*(pos_y-neg_y);
    position[2] = neg_y+random.ran()*(pos_y-neg_y);

    // check that sampled position is within cell
    Ensure (position[0] >= lox);
    Ensure (position[0] <= hix);
    Ensure (position[1] >=-tan_beta*position[0]);
    Ensure (position[1] <= tan_beta*position[0]);
    Ensure (position[2] >=-tan_beta*position[0]);
    Ensure (position[2] <= tan_beta*position[0]);

    // return position
    return position;

}
//---------------------------------------------------------------------------//
/*! 
 * \brief  sample a position in a Sphyramid_Mesh cell according to a slope
 *
 * This function samples a position in a Sphyramid_Mesh cell according to a
 * weighting function (i.e. T^4).  The weighting function is assumed linear
 * in each cell, with the slope and cell-centered value provided (the slope
 * in the y and z direction are assumed to be zero).
 * 
 * \param cell cell number
 * \param random random number object
 * \param slope slope sampling weight function (i.e. T^4)
 * \param center_value cell-centered value of sampling weight function
 * \return location in cell
 */
Sphyramid_Mesh::sf_double Sphyramid_Mesh::sample_pos(int cell, 
						     rng_Sprng &random,
						     sf_double slope,
						     double center_value) const
{
    using global::soft_equiv

    Require (cell >  0);
    Require (cell <= num_cells());
    Check   (this->coord->get_dim() == 3);

    // make sure that there is no y or z slopes
    Require (soft_equiv(slope[1], 0.0));
    Require (soft_equiv(slope[2], 0.0));
	

    // initialize location vector
    sf_double position(this->coord->get_dim(), 0.0);

    // get x-dimension cell extents
    double lox = get_low_x(cell);
    double hix = get_high_x(cell);
    double delta_x = hix-lox;
    Check (lox >= 0.0);
    Check (hix >= 0.0);

    // calculate corresponding y (and z) cell extents
    double loy = lox*this->tan_beta;
    double hiy = hix*this->tan_beta;
    Check (loy >= 0.0);
    Check (hiy >= 0.0);

    // calculate x-value of centroid
    // in the future, a Sphyramid_Mesh member function may perfom this
    double x_center = 4.*lox*lox+8.*lox*hix+12*hix*hix;
    x_center       /= 4.*lox*lox+4.*lox*hix+4.*hix*hix;
    x_center        = x_center*delta_x/4.+lox;
    Check (x_center >= lox);
    Check (x_center <= hix);

    // calculate low and high values of the weighting function
    double low_w  = center_value+slope[0]*(lox-x_center);
    double high_w = center_value+slope[0]*(hix-x_center);

    // check that weighting function is positive throughout the cell
    Insist (low_w  >= 0, "Weighting function is negative!");
    Insist (high_w >= 0, "Weighting function is negative!");

    // >>> sample x position <<<

    // calculate coefficients
    double C1 = slope[0];
    double C2 = low_w-C1*lox;
    double C3 = 0.0;
    double C4 = 0.0;

    // calculate maximum value
    double fmax = 0.0;
    if (slope[0] > -low_w/(2.*lox+3.*delta_x))
    {
	// maximum is at hix
	fmax=C1*pow(hix,3)+C2*pow(hix,2)+C3*hix+C4;
    }
    else if (lox <= 0.0 || slope[0] > -low_w/(2.*lox))
    {
	if (slope[0] < -2.*low_w/(lox+3.*delta_x))
	{
	    // local maximum
	    double xmax = -(2./3)*C2/C1;
	    fmax=C1*pow(xmax,3)+C2*pow(xmax,2)+C3*xmax+C4;
	}
	else
	{
	    // maximum at hix
	    fmax=C1*pow(hix,3)+C2*pow(hix,2)+C3*hix+C4;
	}
    }
    else if (slope[0] > -2.*low_w/lox)
    {
	if (slope[0] < -2.*low_w/(lox+3.*delta_x))
	{
	    // local maximum
	    double xmax = -(2./3)*C2/C1;
	    fmax=C1*pow(xmax,3)+C2*pow(xmax,2)+C3*xmax+C4;
	}
	else
	{
	    // maximum at hix
	    fmax=C1*pow(hix,3)+C2*pow(hix,2)+C3*hix+C4;
	}
    }
    else
    {
	// maximum at lox
	fmax=C1*pow(lox,3)+C2*pow(lox,2)+C3*lox+C4;
    }
    Check (fmax > 0.0);

    // sample x
    position[0] = rtt_mc::sampler::sample_cubic(random, lox, hix, C1, C2, C3,
						C3, fmax);

    // sample y and z value
    double pos_y =  position[0]*this->tan_beta;
    double neg_y = -pos_y;
    Check (pos_y >= loy);
    Check (pos_y <= hiy);

    position[1] = neg_y+random.ran()*(pos_y-neg_y);
    position[2] = neg_y+random.ran()*(pos_y-neg_y);

    // check that sampled position is within cell
    Ensure (position[0] >= lox);
    Ensure (position[0] <= hix);
    Ensure (position[1] >=-tan_beta*position[0]);
    Ensure (position[1] <= tan_beta*position[0]);
    Ensure (position[2] >=-tan_beta*position[0]);
    Ensure (position[2] <= tan_beta*position[0]);


    return position;
}
//---------------------------------------------------------------------------//
/*! 
 * \brief Sample a position on a Sphyramid_Mesh cell face
 * 
 * \param cell cell number
 * \param face face number
 * \param random random number object
 *
 * \return position on face
 */
Sphyramid_Mesh::sf_double Sphyramid_Mesh::sample_pos_on_face(int cell, int face,
							     rng_Sprng &random) const
{
    Require (cell >  0);
    Require (cell <= num_cells());

    Check (this->coord->get_dim() == 3);

    // initialize position vector
    sf_double position(this->coord->get_dim(), 0.0);

    // high x face
    if (face == 2)
    {
      double hix = get_high_x(cell);
      double hiy = this->tan_beta*hix;
      Check (hix >= 0.0);
      
      position[0] = hix;
      position[1] = -hiy+random.ran()*2.*hiy;
      position[2] = -hiy+random.ran()*2.*hiy;

      Ensure (position[1] >= -hiy);
      Ensure (position[1] <=  hiy);
      Ensure (position[2] >= -hiy);
      Ensure (position[2] <=  hiy);
    } 
    else{
	Insist (0, "Can only sample on hir/hix face in Sphyramid_Mesh!");
    }

    return position;
}

    

//---------------------------------------------------------------------------//
/*! 
 * \brief Calculate the neighbors around a cell.
 * 
 * \param cell cell index
 * \return vector containing list of neighboring cells
 */
Sphyramid_Mesh::sf_int Sphyramid_Mesh::get_neighbors(int cell) const
{
    Check (this->layout.num_faces(cell) == 6);

    sf_int neighbors;

    for (int face = 1; face <= this->layout.num_faces(cell); face++)
    {
	neighbors.push_back(this->layout(cell,face));
    }
    Check (neighbors.size() == this->layout.num_faces(cell));

    return neighbors;
}

} // end namespace rtt_mc

#endif // rtt_mc_Sphyramid_Mesh_hh

//---------------------------------------------------------------------------//
//              end of mc/Sphyramid_Mesh.hh
//---------------------------------------------------------------------------//
