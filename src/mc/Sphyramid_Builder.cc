//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/Sphyramid_Builder.cc
 * \author Jeffery Densmore (stolen from RZWedge_builder.cc)
 * \date   Mon Nov  10 7:46:00 2003
 * \brief  Sphyramid_Builder implementation file.
 * \note   Copyright © 2003 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Sphyramid_Builder.hh"
#include "XYZCoord_sys.hh"
#include "Math.hh"

#include <algorithm>
#include <cmath>

namespace rtt_mc
{

//---------------------------------------------------------------------------//
// PUBLIC INTERFACE
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
/*! 
 * \brief Build a Sphyramid mesh from data specified in the OS_Mesh format file.
 * 
 * The builder takes data from the OS_Mesh input format that was parsed in
 * the constructor and builds a Sphyramid mesh.  The builder checks for the
 * existence for the built mesh; thus, a builder can only build one mesh.  To
 * get extra copies of the mesh call the get_Mesh() accessor function.
 *
 *
 */
Sphyramid_Builder::SP_Mesh Sphyramid_Builder::build_Mesh()
{
    Require(!mesh);
    Insist(((coord_system =="r") || (coord_system=="R")),
	   "You are using Sphyramid_Builder, but coord_system is not R!");

    // declare smart pointers
    SP_Coord_sys coord;
    SP_Layout layout;

    //build coordinate and layout objects
    coord  = build_Coord();
    layout = build_Sphyramid_Layout(*coord);

    // build the Sphyramid_Mesh
    mesh = build_Sphyramid_Mesh(coord, *layout);

    // calculate defined surface cells;
    calc_defined_surcells();

    // do some check and return completed mesh
    Ensure (mesh->full_Mesh());
    Ensure (mesh->get_Coord().get_dim() ==3);
    Ensure (mesh->num_cells() == zone.size());

    return mesh;
}

//---------------------------------------------------------------------------//
/*! 
 * \brief Return the cell regions for graphics dumping 
 *
 * 
 * \return cell-sized vector with each element the corresponding region
 */

Sphyramid_Builder::sf_int Sphyramid_Builder::get_regions() const
{
    using std::fill;
    using std::vector;

    sf_int return_regions(zone.size());
    vector<int>::const_iterator cell_itr;

    // if no user defined regions given, simply return region 1 in each cell
    if (regions.size() ==1)
    {
	// make all regions 1
	fill(return_regions.begin(), return_regions.end(), 1);
    }
    else if(regions.size()>1)
    {

	// preliminary declarations
	int zone;
	for(int i=0; i<regions.size();i++)
	    for(int j=0;j<regions[i].size(); j++)
	    {
		// get a zone in this region
		zone=regions[i][j];

		// find the cells in this zone and add them to the region
		for (int k=0; k < cell_zone[zone-1].size();k++)
		    return_regions[cell_zone[zone-1][k]-1]=i+1;
	    }
    }
    else 
	Insist (0, "something wrong with the regions setting!");

    // make sure there are no undefined regions
    cell_itr=std::find(return_regions.begin(),return_regions.end(), 0);
    Insist(cell_itr == return_regions.end(), "Incomplete regions defined!");

    return return_regions;
}
//---------------------------------------------------------------------------//
/*! 
 * \brief Return the list of defined surface cells.
 * 
 */
Sphyramid_Builder::vf_int Sphyramid_Builder::get_defined_surcells() const
{
    Require (mesh);
    return defined_surcells;
}

//---------------------------------------------------------------------------//
// PRIVATE MESH FILE PARSING IMPLEMENTATION
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
/*! 
 * \brief parse input file
 * 
 *
 */

void Sphyramid_Builder::parser()
{

    using std::pow;
    using std::fabs;

    // first open the input file
    const char *file=mesh_file.c_str();
    std_ifstream input(file);

    // make sure the file exists
   Insist (input, "You must supply a mesh input file!");

    // determine coord_sys 
    std_string keyword; 
    while (keyword != "end-title")
    {
	//test that we have not reached end-of-file
	Insist (!input.eof(), "You did not supply an end-title!");

	//do title block input
	input>>keyword;
	if(keyword == "coord:")
	    input >> coord_system;
    }

    // check to make sure an appropriate coord_sys is called
    Insist (coord_system == "xy"  || coord_system == "XY"  || 
	    coord_system == "rz"  || coord_system == "RZ"  ||
	    coord_system == "XYZ" || coord_system == "xyz" ||
	    coord_system == "R"   || coord_system == "r"   ,
	    "Invalid coordinate system choosen!");

    // call sub-parser
   if (coord_system == "r" || coord_system == "R")
   {
       // the boundary cond is always reflecting at r=0 for a Sphyramid mesh
       bnd_cond[0] = "reflect";

       //parse mesh particulars
       parser1D(input);
   }
   else if (coord_system == "xyz" || coord_system == "XYZ" ||
	    coord_system == "xy" || coord_system == "XY" ||
	    coord_system == "RZ" || coord_system == "rz")
       Insist(0,"On input, Sphyramid needs r or R coord_system!");

   //parse the source block that contains surface source information
   source_parser(input);

   // size regions data if not already sized
   if (regions.empty())
       regions.resize(1);

   // >>> calculate fine_edge array <<<

   // determine size of fine_edge array
   int nfine = 0;
   for (int n = 0; n<fine_cells.size(); n++)
       nfine+=fine_cells[n];
   fine_edge.resize(nfine+1);

   // assign edge data to array

   // initialize fine cell counter
   int ifine=0;

   // initalize local-use widths
   double coarse_cell_width  =0.0;
   double first_fine_cell_width =0.0;

   // loop over coarse cells -- i-th coarse cell bracketed by the i-th and
   // i+1st coarse cell edges
   for (int i =0; i < coarse_edge.size() - 1; i++)
   {
       // set the first fine edge -- equivalent to first coarse edge
       fine_edge[ifine] =coarse_edge[i];

       // calculate the width of this coarse cell
       coarse_cell_width=coarse_edge[i+1]-coarse_edge[i];

       // set the width of the first fine cell in this coarse cell
       if (fine_ratio[i] ==1.0)
	   first_fine_cell_width=coarse_cell_width/fine_cells[i];
       else
	   first_fine_cell_width=(1.0-fine_ratio[i]) * coarse_cell_width / 
	       (1.0-pow(fine_ratio[i],fine_cells[i]));

       // initialize the offset from the low edge of the coarse cell
       double coarse_offset=0.0;

       // set each fine edge in this coarse cell
       for(int j=1; j<=fine_cells[i];j++)
       {
	   // increment fine edge counter (for entire direction d)
	   ifine++;

	   // set offset from low edge of this coarse cell; set fine edge
	   if (fine_ratio[i] == 1.0)
	       coarse_offset=j* first_fine_cell_width;
	   else
	       coarse_offset += pow(fine_ratio[i],j-1)*first_fine_cell_width;
	   fine_edge[ifine] = coarse_edge[i]+coarse_offset;
       }

       
       // check last fine cell in this coarse cell; it should be equal to the
       // "first fine cell width" if we use the inverse of the ratio.
       double last_width=fine_edge[ifine]-fine_edge[ifine-1];
       double inv_ratio=1.0 / fine_ratio[i];
       double check_last_w;
       if (inv_ratio == 1.0)
	   check_last_w= coarse_cell_width/ fine_cells[i];
       else
	   check_last_w= (1.-inv_ratio) * coarse_cell_width / 
	       (1.0-pow(inv_ratio,fine_cells[i]));
       Check (fabs(last_width -check_last_w) <1.0e-5 *last_width);
   }

   //reset the last fine edge -- equiv to the last coarse edge.
   fine_edge[ifine] = coarse_edge.back();
}

//---------------------------------------------------------------------------//
/*! 
 * \brief parser for geometry/mesh information
 * 
 *
 */

void Sphyramid_Builder::parser1D(std_ifstream &in)
{

    using std::fill;
    using global::soft_equiv;

    std_string keyword;
    int data;

    // initialization block input
    while (keyword != "end-init")
    {
	// test that we have not reached end-of-file
	Insist (!in.eof(), "You did not specify a proper end-init!");

	// do input
	in >> keyword;
	if (keyword == "num_xcoarse:" || keyword == "num_rcoarse:")
	{
	    in>> data;
	    Require(data>0);
	    fine_cells.resize(data);
	    fine_ratio.resize(data);
	    fill(fine_ratio.begin(), fine_ratio.end(), 1.0);
	    coarse_edge.resize(data+1);
	}
	if (keyword == "lox_bnd:" || keyword == "lor_bnd:")
	{
	    in >> bnd_cond[0];
	    Insist(((coord_system == "r") || (coord_system=="R")) ?
		   bnd_cond[0] == "reflect" : true, "Invalid B.C. for R!");
	}
	if (keyword == "hix_bnd:" || keyword== "hir_bnd:")
	{
	    in>>bnd_cond[1];
	    Require(bnd_cond[1]=="reflect"||bnd_cond[1]=="vacuum");
	}
	if (keyword == "num_regions:")
	{
	    in >> data;
	    Require(data>0);
	    regions.resize(data);
	}
	if (keyword == "num_zones_per_region:")
	{
	    Insist (regions.size()>0, "region size not set!");
	    for (int i = 0; i<regions.size(); i++)
	    {
		in >>data;
		Require(data>0);
		regions[i].resize(data);
	    }
	}
	if (keyword == "regions:")
	{
	    Insist (regions.size()>0, "Region size not set!");
	    in >> data;
	    Require(data>0);
	    Insist (regions[data-1].size()>0, "Zones per region not set!");
	    for (int i=0; i< regions[data-1].size(); i ++){
		in >> regions[data-1][i];
		Require(regions[data-1][i]>0);
	    }
	}
    }

    // mesh block input
    while (keyword != "end-mesh")
    {

	// test that we have not reached end-of-file
	Insist(!in.eof(), "You did not specify a proper end-mesh!");

	// do input
	in >> keyword;
	if (keyword == "xcoarse:" || keyword == "rcoarse:")
	{
	    for (int i=0; i<coarse_edge.size(); i++)
		in >> coarse_edge[i];
	    Insist(soft_equiv(coarse_edge[0],0.0),
		   "The first radial value must be zero!");
	}
	if (keyword == "num_xfine:" || keyword == "num_rfine:"){
	    for(int i = 0; i<fine_cells.size(); i++){
		in >> fine_cells[i];
		Require(fine_cells[i]>0);
	    }
	}
	if (keyword == "xfine_ratio:" || keyword == "rfine_ratio:")
	    for(int i=0;i<fine_cells.size(); i++)
	    {
		in >> fine_ratio[i];
		Require(fine_ratio[i] > 0.0);
	    }
	// unfolding angle for Sphyramid_Mesh
	if (keyword == "cone_angle_degrees:")
	{
	    Insist ((coord_system == "r") || (coord_system =="R"),
		    "Cone angle specified, but coord_sys not equal to R!");

	    in >> alpha_degrees;
	    Insist (alpha_degrees >0.0, "Cone Angle <= 0!");
	    Insist (alpha_degrees <=45.0, "Cone Angle > 45!");
	}
    }
}
//---------------------------------------------------------------------------//
/*! 
 * \brief  Parse source descritptions that are part of the mesh format
 * 
 */
void Sphyramid_Builder::source_parser(std_ifstream &in)
{
    using std::fill;

    //input keywords
    std_string keyword;
    int data;

    //read source info, at present this information is exclusively surface
    //source position information
    while (keyword != "end-source")
    {
	// test that we have not reached end-of-file
	Insist (!in.eof(), "You did not specify a proper end-source!");

	// do input
	in >> keyword;

	// check keyword and read appropriate data

	// get the number of surface sources
	if (keyword == "num_ss:")
	{
	    in >> data;
	    Insist(data==1,"Can only have one surface source in R geometry!");

	    // size surface source data
	    ss_pos.resize(data);
	    num_defined_surcells.resize(data);
	    defined_surcells.resize(data);

	    // initalize to zero
	    fill(num_defined_surcells.begin(),
		 num_defined_surcells.end(), 0);
	}
	// read in the surface source positions
	if (keyword == "sur_source:")
	{

	    // make sure that surface source has been sized
	    Insist (ss_pos.size() !=0, "Number of surface sources =0!");

	    // input the surface source position
	    for (int i=0; i< ss_pos.size(); i++){
		in >> ss_pos[i];
		Insist(ss_pos[i]=="hir",
		       "Surface source can only be on exterior!");
	    }
	}

	// read in the number of defined surface cells for each position
	if (keyword == "num_defined_surcells:")
	{
	    // make sure number of user/host-defined surcells has been sized
	    Insist (num_defined_surcells.size() !=0,
		    "Number of surface sources = 0!");

	    // input the number of user/host defined surcells.
	    // input file requires any preceding zeros
	    for (int i=0; i < num_defined_surcells.size(); i++){
		in >> num_defined_surcells[i];
		Require(num_defined_surcells[i]==0
			||num_defined_surcells[i]==1);
	    }

	    // resize defined_surcells
	    for (int i=0; i< num_defined_surcells.size(); i++)
		defined_surcells[i].resize(num_defined_surcells[i]);
	}
   
	// read in the defined surface cells list; keyword repeated for each
	// surface source with defined surface cells
	if (keyword == "defined_surcells:")
	{
       
	    	
	    // make sure user/host-defined surcells has been sized
	    Insist (defined_surcells.size()!=0,
		"Number of surface sources = 0!");

	    // surface source index (first entry after keyword)
	    in >> data;
	    Insist (data <= defined_surcells.size(),
		    "Defined surface source out of range!");

	    Insist (defined_surcells[data-1].size() == 
		    num_defined_surcells[data-1],
		    "Range error in defined surface source cells!");

	    // input user/host-defined surcells
	    for (int i = 0; i< num_defined_surcells[data-1]; i++)
		in >> defined_surcells[data-1][i];
	}
    }
}

//---------------------------------------------------------------------------//
// PRIVATE MESH BUILDING IMPLEMENTATION
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
/*! 
 * \brief build mesh specifically for Sphyramid_Mesh 
 * 
 * \param coord smart pointer to Coord_sys object
 * \param layout reference to AMR_layout object
 * \return smart pointer to built mesh
 */

Sphyramid_Builder::SP_Mesh
Sphyramid_Builder::build_Sphyramid_Mesh(SP_Coord_sys coord, AMR_Layout &layout)
{
    using global::pi;
    using std::pow;
    using std::tan;
    using std::cos;
    using std::sqrt;
    using std::atan;

    //consistency checks
    Require(alpha_degrees>0.0);
    Require(alpha_degrees<=45.0);

    // variable declarations
    int num_xsur = fine_edge.size();
    int num_xcells=num_xsur-1;
    int num_cells=num_xcells;
    int dimension = coord->get_dim();

    // check some assertions
    Check(layout.num_cells() == num_cells);
    Check(dimension == 3);

    // radius to x-value conversion
    double alpha_radians =pi/180.0*alpha_degrees;
    double r_to_x=pow((2.*(1.-cos(alpha_radians))
		       /(tan(alpha_radians)*tan(alpha_radians))),1./3.);

    // create the x cell extents vector
    vf_double cell_x_extents(num_cells);
    for (int xsur = 0; xsur<num_xcells;xsur++)
    {
	int cell=xsur;
	cell_x_extents[cell].resize(2);
	cell_x_extents[cell][0] = fine_edge[xsur]*r_to_x;
	cell_x_extents[cell][1]=fine_edge[xsur+1]*r_to_x;
    }

    // calculate beta angle of pyramid
    double beta_radians=atan(sqrt(pi)/2.*tan(alpha_radians));

    // create mesh
    SP_Mesh mesh_return(new Pyramid_Mesh(coord,layout,cell_x_extents,
					 beta_radians));

    // return mesh to builder
    return mesh_return;
}
   
//---------------------------------------------------------------------------//
// PRIVATE COORD_SYS BUILDER IMPLEMENTATION
//---------------------------------------------------------------------------//
//---------------------------------------------------------------------------//
/*! 
 * \brief constructs Coord_sys object
 * 
 *
 * \return smart pointer to new Coord_sys object
 */

Sphyramid_Builder::SP_Coord_sys Sphyramid_Builder::build_Coord() const
{
    using rtt_dsxx::SP

    //build coordinate system
    SP_Coord_sys coord;
    Check (coord_system !="xy"  || coord_system !="XY"  ||
	   coord_system !="xyz" || coord_system !="XYZ" ||
	   coord_system !="rz"  || coord_system !="RZ");
    Check (coord_system == "r"  || coord_system =="R");

    // the Sphyramid mesh uses a 3D XYZ coordinate system
    SP<XYZCoord_sys> xyzcoord(new XYZCoord_sys);
    coord=xyzcoord;
    
    // return base class SP to a derived Coord_sys
    return coord;
}

//---------------------------------------------------------------------------//
// PRIVATE LAYOUT BUILDER IMPLEMENTATION
//---------------------------------------------------------------------------//
//---------------------------------------------------------------------------//
/*! 
 * \brief Construct AMR Layout
 * 
 * \param coord  Coord_sys object
 * \return smart pointer to new layout object
 */

Sphyramid_Builder::SP_Layout
Sphyramid_Builder::build_Sphyramid_Layout(const Coord_sys &coord) const
{
    // set size of new Layout
    int size = fine_edge.size()-1;

    //build layout object
    SP_Layout layout(new AMR_Layout(size));

    // set number of faces for each cell in layout.  For the (OS)
    // Sphyramid mesh, the layout must be 3D XYZ, so there are 6 faces/cell.
    for(int i=1;i<=size;i++)
	layout->set_size(i,6);

    // assign cells and faces to Layout;
    assignSphyramid_Layout(*layout);

    // return built Layout
    return layout;
}		   
//---------------------------------------------------------------------------//
/*! 
 * \brief Construct layout specifically
 * 
 * \param layout reference to AMR_layout object
 */

void Sphyramid_Builder::assignSphyramid_Layout(AMR_Layout &layout) const
{

    // 3D map of Mesh
    int num_xcells = fine_edge.size()-1;

    // loop over num_cells and assign cell acros faces
    // 1:x(-), 2:x(+), 3:y(-), 4:y(+), 5:z(-), 6:z(+)
    // the high and low y and z faces are always reflecting.
    // the x bc's will be taken care of later
    for(int cell=1;cell<=layout.num_cells(); cell++)
    {
	layout(cell,1,1)=cell-1;
	layout(cell,2,1)=cell+1;
	layout(cell,3,1)=cell;
	layout(cell,4,1)=cell;
	layout(cell,5,1)=cell;
	layout(cell,6,1)=cell;
    }

    // take care of boundary conditions

    // low x boundary (always reflecting)
    Insist(bnd_cond[0] =="reflect",
	   "Sphyramid Mesh must be reflecting at the center!");
    layout(1,1,1)=1;

    // high x boundary
    Insist(bnd_cond[1]=="vacuum" || bnd_cond[1]=="reflect",
	   "Invalid boundary condition!");
    if(bnd_cond[1] == "vacuum")
	layout(num_xcells,2,1) = 0;
    else if(bnd_cond[1] == "reflect")
	layout(num_xcells,2,1)= num_xcells;
}

//---------------------------------------------------------------------------//
// CELL ZONING FUNCTIONS (PRIVATE IMPLEMENTATION)
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
/*! 
 * \brief Map mesh zone indices to cells
 * 
 */

void Sphyramid_Builder::zone_mapper()
{
    // determine number of zones and accumulated fine_cell data
    int num_cells =1;
    int num_zones =1;

    // calculate number of zones
    num_cells *= fine_edge.size() -1 ;
    num_zones *= coarse_edge.size()-1;

    // calculate accumulated fine_cell data
    accum_cells.resize(fine_cells.size()+1);
    accum_cells[0] = 0;
    for(int i = 0 ; i < fine_cells.size(); i++)
	accum_cells[i+1] = accum_cells[i]+fine_cells[i];

    // size zone arrays
    zone.resize(num_cells);
    cell_zone.resize(num_zones);

    // determine cell-zone map, i.e. map a fine-cell to a zone
    for (int i=1; i<=fine_cells.size(); i++)
	cell_zoner(i);
}
//---------------------------------------------------------------------------//
/*! 
 * \brief  Map cells to zones in 1D mesh.
 * 
 * \param iz zone number
 *
 */

void Sphyramid_Builder::cell_zoner(int iz)
{
    // match a fine-cell to a zone for 1D meshes

    // descriptive variables
    int num_xzones = coarse_edge.size() -1 ;
    int zone_index = iz;
    int num_xcells = fine_edge.size() - 1;

    // loop boundaries
    int starti = 1+accum_cells[iz-1];
    int endi = accum_cells[iz-1]+fine_cells[iz-1];

    // loop over zone and assign zone_index to those fine-cells residing in
    // the zone
    for (int i = starti; i<=endi; i++)
    {
	int cell =i;
	zone[cell-1]=zone_index;
	cell_zone[zone_index-1].push_back(cell);
    }
}
//---------------------------------------------------------------------------//
// CALCULATE DEFINED SURFACE CELLS
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
/*! 
 * \brief Use built global mesh to calculate a list of defined surface cells.
 *
 * This function uses the global mesh to convert a global boundary surface
 * source into a list of explicit, global cells along the boundary.  For
 * example, a "hix" boundary could be defined as a surface source.  This
 * description will be converted into a list of all global cells along the
 * hi-x boundary.
 *
 * Surface source cell lists are only calculated for surface sources that do
 * not already have user-defined lists of cells.
 */
void Sphyramid_Builder::calc_defined_surcells()
{
    Require (mesh);
    Require (mesh->full_Mesh());

    // define number of surface sources
    int num_ss= ss_pos.size();

    //exit if no surface sources
    if (num_ss == 0)
	return ;

    // check the defined surface source list
    Check(defined_surcells.size() == num_ss);

    // check to see if we need to do any work
    for(int i=0;i<num_ss;i++)
    {
	if(defined_surcells[i].size() == 0)
	{
	    //get a vector of surface cells along the boundary
	    defined_surcells[i] = mesh->get_surcells(ss_pos[i]);
	    Check (defined_surcells[i].size() >0 );
	}
	else if (defined_surcells[i].size() >0)
	{
	    Check(mesh->check_defined_surcells(ss_pos[i],defined_surcells[i]));
	    continue;
	}
	else
	{
	    Insist(0,"Number of defined surface source cells < 0!");
	}
    }
}


} // end namespace rtt_mc

//---------------------------------------------------------------------------//
//                 end of Sphyramid_Builder.cc
//---------------------------------------------------------------------------//



