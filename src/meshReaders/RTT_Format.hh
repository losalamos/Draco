//----------------------------------*-C++-*----------------------------------//
// RTT_Format.hh
// Shawn Pautz (TychoMesh.cc original) / B.T. Adams (Extended to RTT_Format.hh)
// 7 June 99
/*! 
 * \file   meshReaders/RTT_Format.hh
 * \author Shawn Pautz/B.T. Adams
 * \date   Mon Jun 7 10:33:26 1999
 * \brief  Header file for RTT_Format library.
 */
//---------------------------------------------------------------------------//
// @> 
//---------------------------------------------------------------------------//

#ifndef __RTT_Format_hh__
#define __RTT_Format_hh__

#include "ds++/Mat.hh"
#include "ds++/SP.hh"
#include <string>
#include <algorithm>
#include <set>
#include <fstream>
#include <iostream>
#include <vector>
#include <map>

using std::string;
using std::set;
using std::ifstream;
using std::vector;
using std::multimap;
using rtt_dsxx::SP;

/*!
 * \brief RTT format namespace.
 *
 * Provides namespace protection for the Draco RTT format utilities used to
 * parse mesh files in the \ref rtt_format_defined and connect the mesh.
 *
 *\sa The \ref rtt_format_overview page presents a summary of the capabilities
 *    provided within the namespace.
 */
namespace rtt_format
{
 
//===========================================================================//
// class RTT_Format - 
//
/*!
 * \brief  A generalized input routine to parse an RTT Format mesh file and
 *         determine the associated connectivity.
 *
 *\sa The RTT_Format class constructor automatically instantiates and executes
 *    the readMesh and calculateConnectivity member classes used to parse the 
 *    mesh data and determine the mesh connectivity, respectively. Accessor 
 *    functions are provided for all of the remaining member classes to allow 
 *    data retrieval. The \ref rtt_format_overview page presents a summary of 
 *    the capabilities provided within the namespace.
 */ 
//
// revision history:
// -----------------
// 0) original (developed by extending the TychoMesh.hh file of Shawn Pautz)
// 
//===========================================================================//

class RTT_Format 
{
/*!
 * \brief RTT_Format nested class member that controls parsing, storing, and 
 *        accessing the data contained in the header section of the mesh file.
 */
    // NESTED CLASSES AND TYPEDEFS

  public:

    class Header
    {
	string version;
	string title;
	string date;
	int cycle;
	double time;
        int ncomments;
	rtt_dsxx::Mat1<string> comments;

      public:
	Header() {}
	~Header() {}

	void readHeader(ifstream & meshfile);

      private:
	void readKeyword(ifstream & meshfile);
	void readData(ifstream & meshfile);
	void readEndKeyword(ifstream & meshfile);

      public:
        // header data access
        string get_version() const { return version; }
        string get_title() const { return title; }
        string get_date() const { return date; }
	int get_cycle() const { return cycle; }
	double get_time() const { return time; }
	int get_ncomments() const { return ncomments; }	
        string get_comments(int i) const { return comments(i); }

    };

/*!
 * \brief RTT_Format nested class member that controls parsing, storing, and 
 *        accessing the data contained in the dimensions section of the mesh 
 *        file.
 */
    class Dims
    {
	string coor_units;
	string prob_time_units;

	int ncell_defs;
	int nnodes_max;
	int nsides_max;
	int nnodes_side_max;

	int ndim;
	int ndim_topo;

	int nnodes;
	int nnode_flag_types;
	rtt_dsxx::Mat1<int> nnode_flags;
	int nnode_data;

	int nsides;
	int nside_types;
	rtt_dsxx::Mat1<int> side_types;
	int nside_flag_types;
	rtt_dsxx::Mat1<int> nside_flags;
	int nside_data;

	int ncells;
	int ncell_types;
	rtt_dsxx::Mat1<int> cell_types;
	int ncell_flag_types;
	rtt_dsxx::Mat1<int> ncell_flags;
	int ncell_data;

        // flag to indicate node, side, and cell renumbering is performed.
        bool renumber;

      public:
	Dims() { ncells = 0; }
	~Dims() {}

	void readDims(ifstream & meshfile, const bool & renumber);

      private:
	void readKeyword(ifstream & meshfile);
	void readUnits(ifstream & meshfile);
	void readCellDefs(ifstream & meshfile);
	void readDimensions(ifstream & meshfile);
	void readNodes(ifstream & meshfile);
	void readSides(ifstream & meshfile);
	void readCells(ifstream & meshfile);
	void readEndKeyword(ifstream & meshfile);

      public:	
        // units and cell definition data access
        string get_coor_units() const { return coor_units; }
        string get_prob_time_units() const { return prob_time_units; }
	int get_ncell_defs() const { return ncell_defs; }
	int get_nnodes_max() const { return nnodes_max; }
	int get_nsides_max() const { return nsides_max; }
	int get_nnodes_side_max() const { return nnodes_side_max; }

        // node data access
	int get_ndim() const { return ndim; }
	int get_ndim_topo() const { return ndim_topo; }
	int get_nnodes() const { return nnodes; }
	int get_nnode_flag_types() const { return nnode_flag_types; }
	int get_nnode_flags(int i) const { return nnode_flags(i); }
	int get_nnode_data() const { return nnode_data; }

        // side data access
	int get_nsides() const { return nsides; }
	int get_nside_types() const { return nside_types; }
	int get_side_types(int i) const  { return  side_types(i); }
	int get_nside_flag_types() const { return nside_flag_types; }
	int get_nside_flags(int i) const { return nside_flags(i); }
	int get_nside_data() const { return nside_data; }

        // cell data access
	int get_ncells() const { return ncells; }
	int get_ncell_types() const { return ncell_types; }	
	int get_cell_types(int i)  const { return  cell_types(i); }	
	int get_ncell_flag_types() const { return ncell_flag_types; }
	int get_ncell_flags(int i) const { return ncell_flags(i); }
	int get_ncell_data() const { return ncell_data; }

        // renumbering flag
        bool get_renumber() const { return renumber; }

	bool allowed_side_type(int sidetype) const
	{ return side_types.end()
	      != std::find(side_types.begin(), side_types.end(), sidetype); }
	bool allowed_cell_type(int celltype) const
	{ return cell_types.end()
	      != std::find(cell_types.begin(), cell_types.end(), celltype); }
    };

/*!
 * \brief RTT_Format nested base class member that controls parsing, storing, 
 *        and accessing the data contained in the node, side, and cell flag 
 *        sections of the mesh file.
 */
    class Flags
    {
	int nflags;
	string name;
	rtt_dsxx::Mat1<int> flag_nums;
	rtt_dsxx::Mat1<string> flag_names;

      public:
	Flags(int nflags_, const string & name_)
	    : nflags(nflags_), name(name_), flag_nums(nflags), flag_names(nflags) {}
	~Flags() {}

	void readFlags(ifstream & meshfile);

      public:
	bool allowed_flag(int flag) const
	{ return flag_nums.end()
	      != std::find(flag_nums.begin(), flag_nums.end(), flag); }

        string getFlagType() const {return name;}
        int getFlagNumber(int flag) const {return flag_nums(flag);}
        string getFlagName(int flag) const {return flag_names(flag);}
        int getFlagSize() const {return nflags;}
    };

/*!
 * \brief RTT_Format nested class member that controls parsing, storing, and
 *        accessing the data specific to the node flags section of the mesh 
 *        file.
 */
    class NodeFlags
    {
	const Dims & dims;
	rtt_dsxx::Mat1<rtt_dsxx::SP<Flags> > flagTypes;

      public:
	NodeFlags(const Dims & dims_)
	    : dims(dims_), flagTypes(dims.get_nnode_flag_types()) {}
	~NodeFlags() {}

	void readNodeFlags(ifstream & meshfile);

      private:
	void readKeyword(ifstream & meshfile);
	void readFlagTypes(ifstream & meshfile);
	void readEndKeyword(ifstream & meshfile);

      public:
	bool allowed_flag(int flagtype, int flag) const
	{
	    Insist(flagtype <= dims.get_nnode_flag_types() - 1,
		   "Invalid node flag type number!");
	    return flagTypes(flagtype)->allowed_flag(flag); 
	}
        string get_flag_type(int flagtype) const 
	{ 
	    Insist(flagtype <= dims.get_nnode_flag_types() - 1,
		   "Invalid node flag type number!");
	    return flagTypes(flagtype)->getFlagType();
	}
        int get_flag_number(int flagtype, int flag_index) const 
	{
	    Insist(flagtype  <= dims.get_nnode_flag_types() - 1,
		   "Invalid node flag type number!");
	    Insist(flag_index  <= flagTypes(flagtype)->getFlagSize() - 1,
		   "Invalid node flag number index number!");
	    return flagTypes(flagtype)->getFlagNumber(flag_index);
	}
        int get_flag_size(int flagtype) const 
	{
	    Insist(flagtype  <= dims.get_nnode_flag_types() - 1,
		   "Invalid node flag type number!");
	    return flagTypes(flagtype)->getFlagSize();
	}
        string get_flag_name(int flagtype, int flag_index) const 
	{
	    Insist(flagtype  <= dims.get_nnode_flag_types() - 1,
		   "Invalid node flag type number!");
	    Insist(flag_index  <= flagTypes(flagtype)->getFlagSize() - 1,
		   "Invalid node flag name index number!");
	    return flagTypes(flagtype)->getFlagName(flag_index);
	}
    };

/*!
 * \brief RTT_Format nested class member that controls parsing, storing, and
 *        accessing the data specific to the side flags section of the mesh 
 *        file.
 */
    class SideFlags
    {
	const Dims & dims;
	rtt_dsxx::Mat1<rtt_dsxx::SP<Flags> > flagTypes;

      public:
	SideFlags(const Dims & dims_)
	    : dims(dims_), flagTypes(dims.get_nside_flag_types()) {}
	~SideFlags() {}

	void readSideFlags(ifstream & meshfile);

      private:
	void readKeyword(ifstream & meshfile);
	void readFlagTypes(ifstream & meshfile);
	void readEndKeyword(ifstream & meshfile);

      public:
	bool allowed_flag(int flagtype, int flag) const
	{
	    Insist(flagtype <= dims.get_nside_flag_types() - 1,
		   "Invalid side flag type number!");
	    return flagTypes(flagtype)->allowed_flag(flag); 
	}

        string get_flag_type(int flagtype) const 
	{ 
	    Insist(flagtype <= dims.get_nside_flag_types() - 1,
		   "Invalid side flag type number!");
	    return flagTypes(flagtype)->getFlagType();
	}
        int get_flag_number(int flagtype, int flag_index) const 
	{
	    Insist(flagtype  <= dims.get_nside_flag_types() - 1,
		   "Invalid side flag type number!");
	    Insist(flag_index  <= flagTypes(flagtype)->getFlagSize() - 1,
		   "Invalid side flag number index number!");
	    return flagTypes(flagtype)->getFlagNumber(flag_index);
	}
        int get_flag_size(int flagtype) const 
	{
	    Insist(flagtype  <= dims.get_nside_flag_types() - 1,
		   "Invalid side flag type number!");
	    return flagTypes(flagtype)->getFlagSize();
	}
        string get_flag_name(int flagtype, int flag_index) const 
	{
	    Insist(flagtype  <= dims.get_nside_flag_types() - 1,
		   "Invalid side flag type number!");
	    Insist(flag_index  <= flagTypes(flagtype)->getFlagSize() - 1,
		   "Invalid side flag name index number!");
	    return flagTypes(flagtype)->getFlagName(flag_index);
	}
        int get_boundary_flag_number() const 
	{
	    // Allow any combination of the phrase "boundary_conditions". 
	    string bc("boundarycditBOUNDARYCDIT_");
	    int boundary = -1;
	    int length = 0;
	    for (int f = 0; f < dims.get_nside_flag_types(); f++)
	    {
	        string flag = flagTypes(f)->getFlagType();
	        if ((flag[0] == 'b' || flag[0] == 'B') &&
		    flag.find_first_not_of(bc) == string::npos &&  
		    flag.find_first_not_of(bc) >= length)
		{
		    length = flag.size();
		    boundary = f;
		}
	    }
	    Insist(boundary >= 0, "Boundary conditions not found!");
	    return boundary;
	}
        int get_surface_src_flag_number() const 
	{
	    // Allow any combination of the phrase "surface_source". 
	    string surface("surfaceoSURFACEO_");
	    int source = -1;
	    int length = 0;
	    for (int f = 0; f < dims.get_nside_flag_types(); f++)
	    {
	        string flag = flagTypes(f)->getFlagType();
	        if ((flag[0] == 's' || flag[0] == 'S') &&
		    flag.find_first_not_of(surface) == string::npos &&  
		    flag.find_first_not_of(surface) >= length)
		{
		    length = flag.size();
		    source = f;
		}
	    }
	    return source;
	}
    };

/*!
 * \brief RTT_Format nested class member that controls parsing, storing, and
 *        accessing the data specific to the cell flags section of the mesh 
 *        file.
 */
    class CellFlags
    {
	const Dims & dims;
	rtt_dsxx::Mat1<rtt_dsxx::SP<Flags> > flagTypes;

      public:
	CellFlags(const Dims & dims_)
	    : dims(dims_), flagTypes(dims.get_ncell_flag_types()) {}
	~CellFlags() {}

	void readCellFlags(ifstream & meshfile);

      private:
	void readKeyword(ifstream & meshfile);
	void readFlagTypes(ifstream & meshfile);
	void readEndKeyword(ifstream & meshfile);

      public:
	bool allowed_flag(int flagtype, int flag) const
	{ 
	    Insist(flagtype <= dims.get_ncell_flag_types() - 1,
		   "Invalid cell flag type number!");
	    return flagTypes(flagtype)->allowed_flag(flag); 
	}

        string get_flag_type(int flagtype) const 
	{ 
	    Insist(flagtype <= dims.get_ncell_flag_types() - 1,
		   "Invalid cell flag type number!");
	    return flagTypes(flagtype)->getFlagType();
	}
        int get_flag_number(int flagtype, int flag_index) const 
	{
	    Insist(flagtype  <= dims.get_ncell_flag_types() - 1,
		   "Invalid cell flag type number!");
	    Insist(flag_index  <= flagTypes(flagtype)->getFlagSize() - 1,
		   "Invalid cell flag number index number!");
	    return flagTypes(flagtype)->getFlagNumber(flag_index);
	}
        int get_flag_size(int flagtype) const 
	{
	    Insist(flagtype  <= dims.get_ncell_flag_types() - 1,
		   "Invalid cell flag type number!");
	    return flagTypes(flagtype)->getFlagSize();
	}
        string get_flag_name(int flagtype, int flag_index) const 
	{
	    Insist(flagtype  <= dims.get_ncell_flag_types() - 1,
		   "Invalid cell flag type number!");
	    Insist(flag_index  <= flagTypes(flagtype)->getFlagSize() - 1,
		   "Invalid cell flag name index number!");
	    return flagTypes(flagtype)->getFlagName(flag_index);
	}
        int get_material_flag_number() const 
	{
	    // Allow any combination of the phrase "material". 
	    string matl("materilMATERIL");
	    int material = -1;
	    int length = 0;
	    for (int f = 0; f < dims.get_ncell_flag_types(); f++)
	    {
	        string flag = flagTypes(f)->getFlagType();
	        if ((flag[0] == 'm' || flag[0] == 'M') &&
		    flag.find_first_not_of(matl) == string::npos &&  
		    flag.find_first_not_of(matl) >= length)
		{
		    length = flag.size();
		    material = f;
		}
	    }
	    Insist(material >= 0, "Material conditions not found!");
	    return material;
	}    
        int get_volume_src_flag_number() const 
	{
	    // Allow any combination of the phrase "volume_source". 
	    string source("volumesrcVOLUMESRC_");
	    int vol_src = -1;
	    int length = 0;
	    for (int f = 0; f < dims.get_ncell_flag_types(); f++)
	    {
	        string flag = flagTypes(f)->getFlagType();
	        if ((flag[0] == 'v' || flag[0] == 'V') &&
		    flag.find_first_not_of(source) == string::npos &&  
		    flag.find_first_not_of(source) >= length)
		{
		    length = flag.size();
		    vol_src = f;
		}
	    }
	    return vol_src;
	}    
        int get_radiation_src_flag_number() const 
	{
	    // Allow any combination of the phrase "raditiation_source". 
	    string source("raditonsuceRADITONSUCE_");
	    int rad_src = -1;
	    int length = 0;
	    for (int f = 0; f < dims.get_ncell_flag_types(); f++)
	    {
	        string flag = flagTypes(f)->getFlagType();
	        if ((flag[0] == 'r' || flag[0] == 'R') &&
		    flag.find_first_not_of(source) == string::npos &&  
		    flag.find_first_not_of(source) >= length)
		{
		    length = flag.size();
		    rad_src = f;
		}
	    }
	    return rad_src;
	}    
    };
/*!
 * \brief RTT_Format nested class member that controls parsing, storing, and
 *        accessing the data specific to the node data ids section of the mesh 
 *        file.
 */
    class NodeDataIDs
    {
	const Dims & dims;
	rtt_dsxx::Mat1<string> names;
	rtt_dsxx::Mat1<string> units;

      public:
	NodeDataIDs(const Dims & dims_)
	    : dims(dims_), names(dims.get_nnode_data()),
	      units(dims.get_nnode_data()) {}
	~NodeDataIDs() {}

	void readDataIDs(ifstream & meshfile);

      private:
	void readKeyword(ifstream & meshfile);
	void readData(ifstream & meshfile);
	void readEndKeyword(ifstream & meshfile);

      public:
        string get_data_id_name(int id_numb) const 
	{ 
	    Insist(id_numb <= dims.get_nnode_data() - 1,
		   "Invalid node data id number!");
	    return names(id_numb);
	}
        string get_data_id_units(int id_numb) const 
	{ 
	    Insist(id_numb <= dims.get_nnode_data() - 1,
		   "Invalid node data id number!");
	    return units(id_numb);
	}
    };

/*!
 * \brief RTT_Format nested class member that controls parsing, storing, and
 *        accessing the data specific to the side data ids section of the mesh 
 *        file.
 */
    class SideDataIDs
    {
	const Dims & dims;
	rtt_dsxx::Mat1<string> names;
	rtt_dsxx::Mat1<string> units;

      public:
	SideDataIDs(const Dims & dims_)
	    : dims(dims_), names(dims.get_nside_data()),
	      units(dims.get_nside_data()) {}
	~SideDataIDs() {}

	void readDataIDs(ifstream & meshfile);

      private:
	void readKeyword(ifstream & meshfile);
	void readData(ifstream & meshfile);
	void readEndKeyword(ifstream & meshfile);

      public:
        string get_data_id_name(int id_numb) const 
	{ 
	    Insist(id_numb <= dims.get_nside_data() - 1,
		   "Invalid side data id number!");
	    return names(id_numb);
	}
        string get_data_id_units(int id_numb) const 
	{ 
	    Insist(id_numb <= dims.get_nside_data() - 1,
		   "Invalid side data id number!");
	    return units(id_numb);
	}
    };

/*!
 * \brief RTT_Format nested class member that controls parsing, storing, and
 *        accessing the data specific to the cell data ids section of the mesh 
 *        file.
 */
    class CellDataIDs
    {
	const Dims & dims;
	rtt_dsxx::Mat1<string> names;
	rtt_dsxx::Mat1<string> units;

      public:
	CellDataIDs(const Dims & dims_)
	    : dims(dims_), names(dims.get_ncell_data()),
	      units(dims.get_ncell_data()) {}
	~CellDataIDs() {}

	void readDataIDs(ifstream & meshfile);

      private:
	void readKeyword(ifstream & meshfile);
	void readData(ifstream & meshfile);
	void readEndKeyword(ifstream & meshfile);

      public:
        string get_data_id_name(int id_numb) const 
	{ 
	    Insist(id_numb <= dims.get_ncell_data() - 1,
		   "Invalid cell data id number!");
	    return names(id_numb);
	}
        string get_data_id_units(int id_numb) const 
	{ 
	    Insist(id_numb <= dims.get_ncell_data() - 1,
		   "Invalid cell data id number!");
	    return units(id_numb);
	}
    };

    class CellDefs;

/*!
 * \brief RTT_Format nested base class member that controls parsing, storing, 
 *        and accessing the data contained in the cell definition section of 
 *        the mesh file.
 */
    class CellDef
    {
	const CellDefs & cellDefs;
	const string name;
	int nnodes;
	int nsides;
	rtt_dsxx::Mat1<int> side_types;
	rtt_dsxx::Mat1<set<int> > sides;
        // Add the capability to maintain the sense of the outward normals.
        vector<vector<int> > ordered_sides;

      public:
	CellDef(const CellDefs & cellDefs_, const string & name_)
	    : name(name_), cellDefs(cellDefs_), ordered_sides(0) {}
	~CellDef() {}

	void readDef(ifstream & meshfile);
        void sortData();

        string get_name() const  { return name; }
	int get_nnodes() const { return nnodes; }
	int get_nsides() const { return nsides; }
        int get_side_types(int s) const { return side_types(s); }
	const set<int> & get_side(int s) const { return sides(s); }
	const vector<int> & get_ordered_side(int s) const 
	{ return ordered_sides[s]; }
    };

/*!
 * \brief RTT_Format nested class member that controls parsing, storing, and
 *        accessing the data contained in the cell definition section of the 
 *        mesh file.
 */
    class CellDefs
    {
	const Dims & dims;
	rtt_dsxx::Mat1<rtt_dsxx::SP<CellDef> > defs;

      public:
	CellDefs(const Dims & dims_)
	    : dims(dims_), defs(dims.get_ncell_defs()) {}
	~CellDefs() {}

	void readCellDefs(ifstream & meshfile);

      private:
	void readKeyword(ifstream & meshfile);
	void readDefs(ifstream & meshfile);
        void sortData();
	void readEndKeyword(ifstream & meshfile);

      public:
	string get_name(int i) const { return defs(i)->get_name(); }
	const CellDef & get_cell_def(int i) const { return *(defs(i)); }
	int get_nnodes(int i) const { return defs(i)->get_nnodes(); }
	int get_nsides(int i) const { return defs(i)->get_nsides(); }
        int get_side_types(int i, int s) const 
	{ return defs(i)->get_side_types(s); }
        const set<int> & get_side(int i, int s) const 
            { return defs(i)->get_side(s); }
        const vector<int> & get_ordered_side(int i, int s) const 
	    { return defs(i)->get_ordered_side(s); }
	int get_ncell_defs() const { return dims.get_ncell_defs(); }

    };

/*!
 * \brief RTT_Format nested class member that controls parsing, storing, and
 *        accessing the data specific to the nodes section of the mesh file.
 */
    class Nodes
    {
	const NodeFlags & nodeFlags;
	const Dims & dims;
	rtt_dsxx::Mat2<double> coords;
	rtt_dsxx::Mat1<int> parents;
	rtt_dsxx::Mat2<int> flags;
        // This vector is a map from the input node numbers (vector index) to
        // the sorted node number (stored value)
        vector<int> sort_map;

      public:
	Nodes(const NodeFlags & nodeFlags_, const Dims & dims_)
	    : nodeFlags(nodeFlags_), dims(dims_),
	      coords(dims.get_nnodes(),	dims.get_ndim()), 
	      parents(dims.get_nnodes()),sort_map(0),
	      flags(dims.get_nnodes(),dims.get_nnode_flag_types()) {}
	~Nodes() {}

	void readNodes(ifstream & meshfile);

      private:
	void readKeyword(ifstream & meshfile);
	void readData(ifstream & meshfile);
        static bool sortXYZ(const vector<double> & low_value, 
			    const vector<double> & high_value);
	void sortData();
	void readEndKeyword(ifstream & meshfile);


      public:
	double get_coords(int node_numb, int coord_index) const
	{ return coords(node_numb,coord_index); }

	vector<double> get_coords(int node_numb) const
	{ 
	    vector<double> local_coords(dims.get_ndim());
	    for (int d = 0; d < dims.get_ndim(); d++)
	        local_coords[d] = coords(node_numb,d);
	    return local_coords; 
	}

        int get_node(vector<double> node_coords) const;

	int get_parents(int cell_numb) const
	{ return parents(cell_numb); }

	int get_flags(int cell_numb, int flag_numb) const
	{ return flags(cell_numb,flag_numb); }

        int get_map(int cell_numb) const
	{ return sort_map[cell_numb];}
    };

/*!
 * \brief RTT_Format nested class member that controls parsing, storing, and
 *        accessing the data specific to the sides section of the mesh file.
 */
    class Sides
    {
	const SideFlags & sideFlags;
	const Dims & dims;
	const CellDefs & cellDefs;
        const Nodes & ptrNodes;
	rtt_dsxx::Mat1<int> sideType;
	rtt_dsxx::Mat2<int> nodes;
	rtt_dsxx::Mat2<int> flags;

      public:
	Sides(const SideFlags & sideFlags_, const Dims & dims_,
	      const CellDefs & cellDefs_, const Nodes & ptrNodes_)
	    : sideFlags(sideFlags_), dims(dims_), cellDefs(cellDefs_),
	      ptrNodes(ptrNodes_), sideType(dims.get_nsides()),
	      nodes(dims.get_nsides(), dims.get_nnodes_side_max()),
	      flags(dims.get_nsides(), dims.get_nside_flag_types()) {}
	~Sides() {}

	void readSides(ifstream & meshfile);

      private:
	void readKeyword(ifstream & meshfile);
	void readData(ifstream & meshfile);
        // Note that these arguements are passed by value for internal use
        // in a sorting routine.
        static bool sortXYZ(vector<int> low_value, vector<int> high_value);
	void sortData();
	void readEndKeyword(ifstream & meshfile);

      public:
	int get_type(int side_numb) const
	{ return sideType(side_numb); }

	int get_nodes(int side_numb,int node_numb) const
	{ return nodes(side_numb,node_numb); }

	int get_flags(int side_numb,int flag_numb) const
	{ return flags(side_numb,flag_numb); }

        int get_boundary_flag_number() const 
	{ return sideFlags.get_boundary_flag_number(); }

        int get_surface_src_flag_number() const 
	{ return sideFlags.get_surface_src_flag_number(); }
    };

/*!
 * \brief RTT_Format nested class member that controls parsing, storing, and
 *        accessing the data specific to the cells section of the mesh file.
 */
    class Cells
    {
	const CellFlags & cellFlags;
	const Dims & dims;
	const CellDefs & cellDefs;
        const Nodes & ptrNodes;
	rtt_dsxx::Mat1<int> cellType;
	rtt_dsxx::Mat2<int> nodes;
	rtt_dsxx::Mat2<int> flags;

      public:
	Cells(const CellFlags & cellFlags_, const Dims & dims_,
	      const CellDefs & cellDefs_, const Nodes & ptrNodes_) 
	    : cellFlags(cellFlags_), dims(dims_), cellDefs(cellDefs_),
	      ptrNodes(ptrNodes_), cellType(dims.get_ncells()),
	      nodes(dims.get_ncells(), dims.get_nnodes_max()),
	      flags(dims.get_ncells(), dims.get_ncell_flag_types()) {}
	~Cells() {}

	void readCells(ifstream & meshfile);

      private:
	void readKeyword(ifstream & meshfile);
	void readData(ifstream & meshfile);
	void sortData();
	void readEndKeyword(ifstream & meshfile);


      public:
	int get_type(int cell_numb) const
	{ return cellType(cell_numb); }

	int get_nodes(int cell_numb,int node_numb) const
	{ return nodes(cell_numb,node_numb); }

	int get_flags(int cell_numb,int flag_numb) const
	{ return flags(cell_numb,flag_numb); }
    };

/*!
 * \brief RTT_Format nested class member that controls parsing, storing, and
 *        accessing the data specific to the node data section of the mesh 
 *        file.
 */
    class NodeData
    {
	const Dims & dims;
	rtt_dsxx::Mat2<double> data;
        const Nodes & ptrNodes;

      public:
	NodeData(const Dims & dims_, const Nodes & ptrNodes_) 
	    : dims(dims_),ptrNodes(ptrNodes_),
	      data(dims.get_nnodes(),dims.get_nnode_data()) {}
	~NodeData() {}

	void readNodeData(ifstream & meshfile);

      private:
	void readKeyword(ifstream & meshfile);
	void readData(ifstream & meshfile);
	void sortData(ifstream & meshfile);
	void readEndKeyword(ifstream & meshfile);

      public:
	double get_data(int node_numb,int data_index) const
	{ return data(node_numb,data_index); }
    };

/*!
 * \brief RTT_Format nested class member that controls parsing, storing, and
 *        accessing the data specific to the side data section of the mesh 
 *        file.
 */
    class SideData
    {
	const Dims & dims;
	rtt_dsxx::Mat2<double> data; 
        const Nodes & ptrNodes;

      public:
	SideData(const Dims & dims_, const Nodes & ptrNodes_)
	    : dims(dims_),ptrNodes(ptrNodes_),
	      data(dims.get_nsides(), dims.get_nside_data()) {}
	~SideData() {}

	void readSideData(ifstream & meshfile);

      private:
	void readKeyword(ifstream & meshfile);
	void readData(ifstream & meshfile);
	void sortData(ifstream & meshfile);
	void readEndKeyword(ifstream & meshfile);

      public:
	double get_data(int side_numb,int data_index) const
	{ return data(side_numb,data_index); }
    };

/*!
 * \brief RTT_Format nested class member that controls parsing, storing, and
 *        accessing the data specific to the cell data section of the mesh 
 *        file.
 */
    class CellData
    {
	const Dims & dims;
	rtt_dsxx::Mat2<double> data; 
        const Nodes & ptrNodes;

      public:
	CellData(const Dims & dims_, const Nodes & ptrNodes_)
	    : dims(dims_),ptrNodes(ptrNodes_),
	      data(dims.get_ncells(), dims.get_ncell_data()) {}
	~CellData() {}

	void readCellData(ifstream & meshfile);

      private:
	void readKeyword(ifstream & meshfile);
	void readData(ifstream & meshfile);
	void sortData(ifstream & meshfile);
	void readEndKeyword(ifstream & meshfile);

      public:
	double get_data(int cell_numb,int data_index) const
	{ return data(cell_numb,data_index); }
    };

/*!
 * \brief RTT_Format nested class member that determines the mesh connectivity
 *        from the input mesh file data.
 */
    class Connectivity
    {
	const Dims & dims;
	const CellDefs & cellDefs;
	const Cells & cells;
	const Sides & sides;
        const Nodes & nodes;
        vector<vector<vector<int> > > adjCell;
        multimap<int, int> bndryFaces;
        multimap<int, vector<int> > Side_to_Cell_Face;

      public:
	Connectivity(const Dims & dims_, const CellDefs & cellDefs_,
		     const Cells & cells_, const Sides & sides_,
		     const Nodes & nodes_);
	~Connectivity() {}

	int get_adjCell(int cell, int face, int adjcell = 0) const
	{ return adjCell[cell][face][adjcell]; }
	int get_adjCell_size(int cell, int face) const
	{ return adjCell[cell][face].size(); }
	int get_bndryFaces_count(int face) const
	{ return bndryFaces.count(face); }
	set<int> get_bndryCells(int face) const
	{ 
	    set<int> bndryCells;
	    multimap<int, int>::const_iterator citer = bndryFaces.find(face);
	    int cell_index = 0;
	    int face_index = citer->first;
	    
	    while (citer !=  bndryFaces.end() && citer->first == face_index)
	    {
	        bndryCells.insert(citer->second);
		++cell_index;
		++citer;
	    }
	    return bndryCells;
	}
	bool check_bndryFace(int cell, int face) const
	{
	    multimap<int, int>::const_iterator citer = bndryFaces.find(face);
	    int face_index = citer->first;
	    bool bndryFace = false;
	    
	    while (citer !=  bndryFaces.end() && citer->first == face_index &&
		   !bndryFace)
	    {
	        if (citer->second == cell)
		    bndryFace = true;
		else
		    ++citer;
	    }
	    return bndryFace;
	}
        int get_Cell_from_Side(int side) const
	{
	    multimap<int, vector<int> >::const_iterator sitr = 
	        Side_to_Cell_Face.find(side);

	    if (sitr != Side_to_Cell_Face.end())
	    {
	        vector<int> Cell_Face = sitr->second;
		return Cell_Face[0];
	    }
	    else
	        return -1;
	}
        int get_Cell_Face_from_Side(int side) const
	{
	    multimap<int, vector<int> >::const_iterator sitr = 
	        Side_to_Cell_Face.find(side);

	    if (sitr != Side_to_Cell_Face.end())
	    {
	        vector<int> Cell_Face = sitr->second;
		return Cell_Face[1];
	    }
	    else
	        return -1;
	}

      private:
	void calcAdjacentCells();
    };

    // DATA

  private:
    Header header;
    Dims dims;
    SP<NodeFlags> spNodeFlags;
    SP<SideFlags> spSideFlags;
    SP<CellFlags> spCellFlags;
    SP<NodeDataIDs> spNodeDataIds;
    SP<SideDataIDs> spSideDataIds;
    SP<CellDataIDs> spCellDataIds;
    SP<CellDefs> spCellDefs;
    SP<Nodes> spNodes;
    SP<Sides> spSides;
    SP<Cells> spCells;
    SP<NodeData> spNodeData;
    SP<SideData> spSideData;
    SP<CellData> spCellData;
    SP<Connectivity> spConnectivity;

  public:

    // Constructors
    RTT_Format(const string & RTT_File, const bool & renumber);

/*!
 * \brief Destroys an RTT_Format class object
 */
    // Destructors
    ~RTT_Format() {}

    // ACCESSORS
    // header data access
/*!
 * \brief Returns the version number read from the mesh file header.
 * \return Mesh file version number.
 */
    string get_header_version() const { return header.get_version(); }
/*!
 * \brief Returns the title read from the mesh file header.
 * \return Mesh file title.
 */
    string get_header_title() const { return header.get_title(); }
/*!
 * \brief Returns the date read from the mesh file header.
 * \return Date the mesh file was generated.
 */
    string get_header_date() const { return header.get_date(); }
/*!
 * \brief Returns the cycle number read from the mesh file header.
 * \return Cycle number associated with the data in the mesh file.
 */
    int get_header_cycle() const { return header.get_cycle(); }
/*!
 * \brief Returns the problem time read from the mesh file header.
 * \return Problem time associated with the data in the mesh file.
 */
    double get_header_time() const { return header.get_time(); }
/*!
 * \brief Returns the number of comment lines read from the mesh file header.
 * \return The number of comment lines contained in the mesh file.
 */
    int get_header_ncomments() const { return header.get_ncomments(); }	
/*!
 * \brief Returns the specified comment line read from the mesh file header.
 * \param i Line number of the comment to be returned.
 * \return The specified comment line contained in the mesh file.
 */
    string get_header_comments(int i) const { return header.get_comments(i); }

    // dimensions units and cell definition data access
/*!
 * \brief Returns the problem coordinate units (e.g, cm) read from the mesh 
 *        file dimension data.
 * \return Mesh file coordinate units.
 */
    string get_dims_coor_units() const { return dims.get_coor_units(); }
/*!
 * \brief Returns the problem time units (e.g, shakes) read from the mesh 
 *        file dimension data.
 * \return Mesh file time units.
 */
    string get_dims_prob_time_units() const 
    { return dims.get_prob_time_units(); }
/*!
 * \brief Returns the total number of cell type definitions read from the mesh 
 *        file dimension data.
 * \return The number of unique cell definitions in the mesh file.
 */
    int get_dims_ncell_defs() const { return dims.get_ncell_defs(); }
/*!
 * \brief Returns the maximum number of nodes per cell type read from the mesh 
 *        file dimension data.
 * \return The maximum number of nodes per cell type.
 */
    int get_dims_nnodes_max() const { return dims.get_nnodes_max(); }
/*!
 * \brief Returns the maximum number of sides per cell type read from the mesh 
 *        file dimension data.
 * \return The maximum number of sides per cell type.
 */
    int get_dims_nsides_max() const { return dims.get_nsides_max(); }
/*!
 * \brief Returns the maximum number of nodes per cell side read from the mesh 
 *        file dimension data.
 * \return The maximum number of nodes per cell side.
 */
    int get_dims_nnodes_side_max() const { return dims.get_nnodes_side_max(); }

    // dimensions node data access
/*!
 * \brief Returns the number of spatial dimensions read from the mesh file 
 *        dimension data.
 * \return The number of spatial dimensions in the mesh.
 */
    int get_dims_ndim() const { return dims.get_ndim(); }
/*!
 * \brief Returns the number of topological dimensions read from the mesh file 
 *        dimension data.
 * \return The number of topological dimensions.
 */
    int get_dims_ndim_topo() const { return dims.get_ndim_topo(); }
/*!
 * \brief Returns the number of nodes read from the mesh file dimension data.
 * \return The number of nodes in the mesh.
 */
    int get_dims_nnodes() const { return dims.get_nnodes(); }
/*!
 * \brief Returns the number of node flag types read from the mesh file 
 *        dimension data.
 * \return The number of node flag types.
 */
    int get_dims_nnode_flag_types() const 
    { return dims.get_nnode_flag_types(); }
/*!
 * \brief Returns the number of node flags for the specified node flag type
 *        read from the mesh file dimension data.
 * \param i Node flag type number.
 * \return The number of node flags for the specified flag type number.
 */
    int get_dims_nnode_flags(int i) const { return dims.get_nnode_flags(i); }
/*!
 * \brief Returns the number of node data fields read from the mesh file 
 *        dimension data.
 * \return The number of node data fields.
 */
    int get_dims_nnode_data() const { return dims.get_nnode_data(); }

    // dimensions side data access
/*!
 * \brief Returns the number of sides read from the mesh file dimension data.
 * \return The number of sides in the mesh.
 */
    int get_dims_nsides() const { return dims.get_nsides(); }
/*!
 * \brief Returns the number of side types that are present in the "sides"
 *        block read from the mesh file dimension data.
 * \return The number of side types that are present in the "sides" block.
 */
    int get_dims_nside_types() const { return dims.get_nside_types(); }
/*!
 * \brief Returns the side type index for the specified side type read from 
 *        the mesh file dimension data.
 * \param i Side type number.
 * \return The side type index for the specified side type.
 */
    int get_dims_side_types(int i) const  { return  dims.get_side_types(i); }
/*!
 * \brief Returns the number of side flag types read from the mesh file 
 *        dimension data.
 * \return The number of side flag types.
 */
    int get_dims_nside_flag_types() const 
    { return dims.get_nside_flag_types(); }
/*!
 * \brief Returns the number of side flags for the specified side flag type
 *        read from the mesh file dimension data.
 * \param i Side flag type number.
 * \return The number of side flags for the specified flag type number.
 */
    int get_dims_nside_flags(int i) const { return dims.get_nside_flags(i); }
/*!
 * \brief Returns the number of side data fields read from the mesh file 
 *        dimension data.
 * \return The number of side data fields.
 */
    int get_dims_nside_data() const { return dims.get_nside_data(); }

    // dimensions cell data access
/*!
 * \brief Returns the number of cells read from the mesh file dimension data.
 * \return The number of cells in the mesh.
 */
    int get_dims_ncells() const { return dims.get_ncells(); }
/*!
 * \brief Returns the number of cell types that are present in the "cells"
 *        block read from the mesh file dimension data.
 * \return The number of cell types that are present in the "cells" block.
 */
    int get_dims_ncell_types() const { return dims.get_ncell_types(); }    
/*!
 * \brief Returns the cell type index for the specified cell type read from 
 *        the mesh file dimension data.
 * \param i Cell type number.
 * \return The cell type index for the specified cell type.
 */
    int get_dims_cell_types(int i)  const { return dims.get_cell_types(i); }
/*!
 * \brief Returns the number of cell flag types read from the mesh file 
 *        dimension data.
 * \return The number of cell flag types.
 */
    int get_dims_ncell_flag_types() const 
    { return dims.get_ncell_flag_types(); }
/*!
 * \brief Returns the number of cell flags for the specified cell flag type
 *        read from the mesh file dimension data.
 * \param i Cell flag type number.
 * \return The number of cell flags for the specified flag type number.
 */
    int get_dims_ncell_flags(int i) const { return dims.get_ncell_flags(i); }
/*!
 * \brief Returns the number of cell data fields read from the mesh file 
 *        dimension data.
 * \return The number of cell data fields.
 */
    int get_dims_ncell_data() const { return dims.get_ncell_data(); }

    // dimensions renumbering flag
/*!
 * \brief Returns the status of the flag indicating that the node, side, and 
 *        cell numbers are to be reassigned based upon coordinates in ascending
 *        order (x, y, and then z).
 * \return The status of the renumbering flag.
 */
    bool get_dims_renumber() const { return dims.get_renumber(); }

    // node flags access
/*!
 * \brief Returns the name of specified node flag type read from the mesh file
 *        node_flags data.
 * \param flagtype Node flag type number.
 * \return The name of the specified node flag type.
 */
    string get_node_flags_flag_type(int flagtype) const 
	{ return spNodeFlags->get_flag_type(flagtype); }
/*!
 * \brief Returns the node flag number associated with the specified node flag
 *        type and node flag index read from the mesh file node_flags data.
 * \param flagtype Node flag type number.
 * \param flag_index Node flag index.
 * \return The node flag number.
 */
    int get_node_flags_flag_number(int flagtype, int flag_index) const 
	{ return spNodeFlags->get_flag_number(flagtype, flag_index); }
/*!
 * \brief Returns the number of node flags for the specified node flag type 
 *        read from the mesh file node_flags data.
 * \param flagtype Node flag type number.
 * \return The number of node flags.
 */
    int get_node_flags_flag_size(int flagtype) const 
	{ return spNodeFlags->get_flag_size(flagtype); }
/*!
 * \brief Returns the node flag name associated with the specified node flag
 *        type and node flag type index from the mesh file node_flags data.
 * \param flagtype Node flag type number.
 * \param flag_index Node flag index.
 * \return The node flag name.
 */
    string get_node_flags_flag_name(int flagtype, int flag_index) const 
	{ return spNodeFlags->get_flag_name(flagtype, flag_index); }

    // side flags access
/*!
 * \brief Returns the name of specified side flag type read from the mesh file
 *        side_flags data.
 * \param flagtype Side flag type number.
 * \return The name of the specified side flag type.
 */
    string get_side_flags_flag_type(int flagtype) const 
	{ return spSideFlags->get_flag_type(flagtype); }
/*!
 * \brief Returns the side flag number associated with the specified side flag
 *        type and side flag index read from the mesh file side_flags data.
 * \param flagtype Side flag type number.
 * \param flag_index Side flag index.
 * \return The side flag number.
 */
    int get_side_flags_flag_number(int flagtype, int flag_index) const 
	{ return spSideFlags->get_flag_number(flagtype, flag_index); }
/*!
 * \brief Returns the number of side flags for the specified side flag type
 *        read from the mesh file side_flags data.
 * \param flagtype Side flag type number.
 * \return The number of side flags.
 */
    int get_side_flags_flag_size(int flagtype) const 
	{ return spSideFlags->get_flag_size(flagtype); }
/*!
 * \brief Returns the side flag name associated with the specified side flag
 *        index and side flag type read from the mesh file side_flags data.
 * \param flagtype Side flag type number.
 * \param flag_index Side flag index.
 * \return The side flag name.
 */
    string get_side_flags_flag_name(int flagtype, int flag_index) const 
	{ return spSideFlags->get_flag_name(flagtype, flag_index); }
/*!
 * \brief Returns the index to the required side flag type that contains the 
 *        problem boundary conditions read from the mesh file side_flags data.
 * \return The index to the boundary conditions side flag type.
 */
    int get_side_flags_boundary_flag_number() const 
	{ return spSideFlags->get_boundary_flag_number(); }
/*!
 * \brief Returns the index to the optional side flag type that contains the 
 *        problem external sources read from the mesh file side_flags data.
 * \return The index to the external source side flag type.
 */
    int get_side_flags_surface_src_flag_number() const 
	{ return spSideFlags->get_surface_src_flag_number(); }

    // cell flags access
/*!
 * \brief Returns the name of specified cell flag type read from the mesh file
 *        cell_flags data.
 * \param flagtype Cell flag type number.
 * \return The name of the specified cell flag type.
 */
    string get_cell_flags_flag_type(int flagtype) const 
	{ return spCellFlags->get_flag_type(flagtype); }
/*!
 * \brief Returns the cell flag number associated with the specified cell flag
 *        type and cell flag index read from the mesh file cell_flags data.
 * \param flagtype Cell flag type number.
 * \param flag_index Cell flag index.
 * \return The cell flag number.
 */
    int get_cell_flags_flag_number(int flagtype, int flag_index) const 
	{ return spCellFlags->get_flag_number(flagtype, flag_index); }
/*!
 * \brief Returns the number of cell flags for the specified cell flag type
 *        read from the mesh file cell_flags data.
 * \param flagtype Cell flag type number.
 * \return The number of cell flags.
 */
    int get_cell_flags_flag_size(int flagtype) const 
	{ return spCellFlags->get_flag_size(flagtype); }
/*!
 * \brief Returns the cell flag name associated with the specified cell flag
 *        type and cell flag index read from the mesh file cell_flags data.
 * \param flagtype Cell flag type number.
 * \param flag_index Cell flag index.
 * \return The cell flag name.
 */
    string get_cell_flags_flag_name(int flagtype, int flag_index) const 
	{ return spCellFlags->get_flag_name( flagtype, flag_index); }
/*!
 * \brief Returns the index to the required cell flag type that contains the 
 *        cell materials read from the mesh file cell_flags data.
 * \return The index to the material cell flag type.
 */
    int get_cell_flags_material_flag_number() const 
	{ return spCellFlags->get_material_flag_number(); }
/*!
 * \brief Returns the index to the optional cell flag type that contains the 
 *        cell volumetric sources read from the mesh file cell_flags data.
 * \return The index to the volumetric source cell flag type.
 */
    int get_cell_flags_volume_src_flag_number() const 
	{ return spCellFlags->get_volume_src_flag_number(); }
/*!
 * \brief Returns the index to the optional cell flag type that contains the 
 *        cell radiation sources read from the mesh file cell_flags data.
 * \return The index to the radiation source cell flag type.
 */
    int get_cell_flags_radiation_src_flag_number() const 
	{ return spCellFlags->get_radiation_src_flag_number(); }

    // node data ids access
/*!
 * \brief Returns the specified node_data_id name read from the mesh file 
 *        node_data_id data.
 * \param id_numb node_data_id index number.
 * \return The node_data_id name.
 */
    string get_node_data_id_name(int id_numb) const 
	{ return spNodeDataIds->get_data_id_name(id_numb) ; }
/*!
 * \brief Returns the units associated with the specified node_data_id read 
 *        from the mesh file node_data_id data.
 * \param id_numb node_data_id index number.
 * \return The node_data_id units.
 */
    string get_node_data_id_units(int id_numb) const 
	{ return spNodeDataIds->get_data_id_units(id_numb); }

    // side data ids access
/*!
 * \brief Returns the specified side_data_id name read from the mesh file 
 *        side_data_id data.
 * \param id_numb side_data_id index number.
 * \return The side_data_id name.
 */
    string get_side_data_id_name(int id_numb) const 
	{ return spSideDataIds->get_data_id_name(id_numb) ; }
/*!
 * \brief Returns the units associated with the specified side_data_id read 
 *        from the mesh file side_data_id data.
 * \param id_numb side_data_id index number.
 * \return The side_data_id units.
 */
    string get_side_data_id_units(int id_numb) const 
	{ return spSideDataIds->get_data_id_units(id_numb); }

    // cell data ids access
/*!
 * \brief Returns the specified cell_data_id name read from the mesh file 
 *        cell_data_id data.
 * \param id_numb cell_data_id index number.
 * \return The cell_data_id name.
 */
    string get_cell_data_id_name(int id_numb) const 
	{ return spCellDataIds->get_data_id_name(id_numb) ; }
/*!
 * \brief Returns the units associated with the specified cell_data_id read 
 *        from the mesh file cell_data_id data.
 * \param id_numb cell_data_id index number.
 * \return The cell_data_id units.
 */
    string get_cell_data_id_units(int id_numb) const 
	{ return spCellDataIds->get_data_id_units(id_numb); }

    // cell definitions access
/*!
 * \brief Returns the name of the specified cell definition read from the mesh
 *        file cell_defs data.
 * \param i Cell definition index number.
 * \return The cell definition name.
 */
    string get_cell_defs_name(int i) const { return spCellDefs->get_name(i); }
/*!
 * \brief Returns the number of nodes associated with the specified cell 
 *        definition read from the mesh file cell_defs data.
 * \param i Cell definition index number.
 * \return The number of nodes comprising the cell definition.
 */
    int get_cell_defs_nnodes(int i) const { return spCellDefs->get_nnodes(i); }
/*!
 * \brief Returns the number of sides associated with the specified cell 
 *        definition read from the mesh file cell_defs data.
 * \param i Cell definition index number.
 * \return The number of sides comprising the cell definition.
 */
    int get_cell_defs_nsides(int i) const { return spCellDefs->get_nsides(i); }
/*!
 * \brief Returns the side type number associated with the specified side 
 *        index and cell definition read from the mesh file cell_defs data.
 * \param i Cell definition index number.
 * \param s Side index number.
 * \return The side type number.
 */
    int get_cell_defs_side_types(int i, int s) const 
        { return spCellDefs->get_side_types(i,s); }
/*!
 * \brief Returns the side definition associated with the specified cell 
 *        definition and side index read from the mesh file cell_defs data 
 *        with the returned cell-node indexes in sorted order.
 * \param i Cell definition index number.
 * \param s Side index number.
 * \return The side definition (i.e., the cell-node indexes that comprise the 
 *         side in sorted order).
 */
    const set<int> & get_cell_defs_side(int i, int s) const 
        { return spCellDefs->get_side(i,s); }
/*!
 * \brief Returns the side definition associated with the specified cell  
 *        definition and side index read from the mesh file cell_defs data 
 *        with the returned cell-node indexes ordered to preserve the right 
 *        hand rule for the outward-directed normal).
 * \param i Cell definition index number.
 * \param s Side index number.
 * \return The side definition (i.e., the cell-node indexes that comprise the 
 *         side ordered to preserve the right hand rule for the 
 *         outward-directed normal).
 */
    const vector<int> & get_cell_defs_ordered_side(int i, int s) const 
	{ return spCellDefs->get_ordered_side(i,s); }

    // nodes access
/*!
 * \brief Returns the coordinate value for the specified node and direction 
 *        (i.e., x, y, and z) read from the mesh file node data.
 * \param node_numb Node number.
 * \param coord_index Coordinate index number (x = 0, y = 1, etc.).
 * \return The node coordinate value.
 */
    double get_nodes_coords(int node_numb, int coord_index) const
	{ return spNodes->get_coords(node_numb, coord_index); }
/*!
 * \brief Returns all of the coordinate values for the specified node read 
 *        from the mesh file node data.
 * \param node_numb Node number.
 * \return The node coordinate values.
 */
    vector<double> get_nodes_coords(int node_numb) const
	{ return spNodes->get_coords(node_numb); }

/*!
 * \brief Returns the node number that has the specified coordinate values.
 * \param node_coords Coordinate values.
 * \return The node number.
 */
    int get_nodes_node(vector<double> node_coords) const
        { return spNodes->get_node(node_coords); }
/*!
 * \brief Returns the node parent for the specified node read from the mesh 
 *        file node data.
 * \param node_numb Node number.
 * \return The node parent.
 */
    int get_nodes_parents(int node_numb) const
	{ return spNodes->get_parents(node_numb); }
/*!
 * \brief Returns the node flag for the specified node and flag read from the 
 *        mesh file node data.
 * \param node_numb Node number.
 * \param flag_numb Node flag number.
 * \return The node flag.
 */
    int get_nodes_flags(int node_numb, int flag_numb) const
	{ return spNodes->get_flags(node_numb, flag_numb); }

/*!
 * \brief Returns the new node number after sorting has been performed when
 *        the renumber flag is set true.
 * \param node_numb Original node number.
 * \return New node number.
 */
    int get_nodes_map(int node_numb) const
        { return spNodes->get_map(node_numb);}

    // sides access
/*!
 * \brief Returns the side type associated with the specified side read from 
 *        the mesh file side data.
 * \param side_numb Side number.
 * \return The side type.
 */
    int get_sides_type(int side_numb) const
	{ return spSides->get_type(side_numb); }
/*!
 * \brief Returns the node number associated with the specified side and 
 *        side-node index read from the mesh file side data.
 * \param side_numb Side number.
 * \param node_numb Side-node index number.
 * \return The node number.
 */
    int get_sides_nodes(int side_numb,int node_numb) const
	{ return spSides->get_nodes(side_numb, node_numb); }
/*!
 * \brief Returns the side flag for the specified side and flag read from the 
 *        mesh file side data.
 * \param side_numb Side number.
 * \param flag_numb Side flag number.
 * \return The side flag.
 */
    int get_sides_flags(int side_numb,int flag_numb) const
	{ return spSides->get_flags(side_numb, flag_numb); }

    // cells access
/*!
 * \brief Returns the cell type associated with the specified cell read from 
 *        the mesh file cell data.
 * \param cell_numb Cell number.
 * \return The cell type.
 */
    int get_cells_type(int cell_numb) const
	{ return spCells->get_type(cell_numb); }

/*!
 * \brief Returns the node number associated with the specified cell and 
 *        cell-node index read from the mesh file cell data.
 * \param cell_numb Cell number.
 * \param node_numb Cell-node index number.
 * \return The node number.
 */
    int get_cells_nodes(int cell_numb,int node_numb) const
	{ return spCells->get_nodes(cell_numb, node_numb); }
/*!
 * \brief Returns the cell flag for the specified cell and flag read from the 
 *        mesh file cell data.
 * \param cell_numb Cell number.
 * \param flag_numb Cell flag number.
 * \return The cell flag.
 */
    int get_cells_flags(int cell_numb,int flag_numb) const
	{ return spCells->get_flags(cell_numb, flag_numb); }

    // connectivity access
/*!
 * \brief Returns the number of the cell adjacent to the specified cell, face,
 *        and optional adjacent cell index (to allow multiple cells to connect
 *        to a single cell face for AMR type meshes).
 * \param cell Cell number.
 * \param face Face number.
 * \param adjcell Adjacent cell number (defaults to 0).
 * \return The adjacent cell number.
 */
    int get_adjCell(int cell, int face, int adjcell = 0) const
        { return spConnectivity->get_adjCell(cell, face, adjcell); }
/*!
 * \brief Returns the number of cells adjacent to the specified cell face 
 *        (allows multiple cells to connect to a single cell face for AMR type
 *         meshes).
 * \param cell Cell number.
 * \param face Face number.
 * \return The number of adjacent cells.
 */
    int get_adjCell_size(int cell, int face) const
        { return spConnectivity->get_adjCell_size(cell,face); }
/*!
 * \brief Returns the number of boundary faces (i.e., faces that are either
 *        on the outer boundary of the problem geometry or a connection between
 *        cells with different refinement levels in an AMR mesh) with the 
 *        specified face number
 * \param face Face number.
 * \return The number of boundary faces.
 */
    int get_bndryFaces_count(int face) const
	{ return spConnectivity->get_bndryFaces_count(face); }
/*!
 * \brief Returns the cells that have boundary faces (i.e., faces that are 
 *        either on the outer boundary of the problem geometry or a connection
 *        between cells with different refinement levels in an AMR mesh) with 
 *        the specified face number
 * \param face Face number.
 * \return The cells with these boundary faces.
 */
    set<int> get_bndryCells(int face) const
        { return spConnectivity->get_bndryCells(face); }
/*!
 * \brief Returns true if the specified cell face is a boundary face (i.e., 
 *        a faces that is either on the outer boundary of the problem geometry
 *        or a connection between cells with different refinement levels in an
 *        AMR mesh).
 * \param cell Cell number.
 * \param face Face number.
 * \return Boundary face status.
 */
    bool check_bndryFace(int cell, int face) const
        { return spConnectivity->check_bndryFace(cell,face); }
/*!
 * \brief Returns the cell number associated with the specified side number.
 * \param side Side number.
 * \return The cell number.
 */
    int get_Cell_from_Side(int side) const
        { return spConnectivity->get_Cell_from_Side(side); }
/*!
 * \brief Returns the cell face number associated with the specified side 
 *        number.
 * \param side Side number.
 * \return The face number.
 */
    int get_Cell_Face_from_Side(int side) const
        { return spConnectivity->get_Cell_Face_from_Side(side); }

    // IMPLEMENTATION

  private:
    
    void readMesh (const string & RTT_file, const bool & renumber = false);
    void readKeyword(ifstream & meshfile);
    void createMembers();
    void readFlagBlocks(ifstream & meshfile);
    void readDataIDs(ifstream & meshfile);
    void readEndKeyword(ifstream & meshfile);
    void calculateConnectivity();
};

} // end namespace rtt_format

#endif                          // __RTT_Format_hh__

/*!
 * \page rtt_format_overview Overview of the RTT_Format package
 *
 * \version 1_0_0
 *
 * <h3> Introduction </h3>
 * The RTT_Format package consists of: private member functions that are used 
 * to parse a mesh file in the \ref rtt_format_defined and determine the mesh
 * connectivity, and public member functions that are used to access the data.
 * The RTT_Format package was developed from the TychoMesh class that was 
 * originally developed by Shawn Pautz with the following additions:
 * <ul>
 *  <li> Extension of the tetrahedral-limited connectivity algorithm to 
 *       arbitrary cell definitions including adaptive mesh refinement (amr)
 *       capability,
 *  <li> Public accessor functions for all of the mesh file data,
 *  <li> Reassignment capability to allow the nodes, sides, and cells to be
 *       renumbered in ascending order based upon their coordinate values, 
 *       with x increasing first, then y, and finally z.
 * </ul> 
 *
 * <h3> Intended Usage </h3>
 * The RTT_Format class constructor automatically parses the specified input
 * file and determines the mesh connectivity via calls to the private member 
 * functions readMesh and calculateConnectivity, respectively. The mesh data 
 * can then be accessed using the public member accessor functions. A second, 
 * optional argument to the RTT_Format class can be used to specify 
 * renumbering. The class constructor defaults to no renumbering, but runtime 
 * testing indicates that there is a significant DECREASE in the total time 
 * required to read and connect the mesh if a continuous adaptive refinement
 * mesh is used with a few thousand cells or more and the renumbering option 
 * is implemented. This behaviour results from the fact that a bilinear search
 * routine can be used to connect adjacent cells with different refinement 
 * levels when renumbering is implemented, while a linear search of all the
 * nodes must be performed otherwise. The RTT_Format class contains several
 * nested classes that correspond primarily to the organization of the \ref
 * rtt_format_defined with the additions of some base classes and the 
 * Connectivity class:
 * <ul>
 *  <li> rtt_format::RTT_Format::Header
 *  <li> rtt_format::RTT_Format::Dims (dimensions)
 *  <li> rtt_format::RTT_Format::Flags (base class of NodeFlags, SideFlags,
 *                               and CellFlags)
 *  <li> rtt_format::RTT_Format::NodeFlags
 *  <li> rtt_format::RTT_Format::SideFlags
 *  <li> rtt_format::RTT_Format::CellFlags
 *  <li> rtt_format::RTT_Format::NodeDataIDs
 *  <li> rtt_format::RTT_Format::SideDataIDs
 *  <li> rtt_format::RTT_Format::CellDataIDs
 *  <li> rtt_format::RTT_Format::CellDef (base class of CellDefs)
 *  <li> rtt_format::RTT_Format::CellDefs (cell definitions)
 *  <li> rtt_format::RTT_Format::Nodes
 *  <li> rtt_format::RTT_Format::Sides
 *  <li> rtt_format::RTT_Format::Cells
 *  <li> rtt_format::RTT_Format::NodeData
 *  <li> rtt_format::RTT_Format::SideData
 *  <li> rtt_format::RTT_Format::CellData
 *  <li> rtt_format::RTT_Format::Connectivity
 * </ul> 
 * These contained classes provide a convenient grouping of the mesh data, and
 * the RTT_Format public member accessor functions reflect the name of the 
 * associated nested class.
 *
 */

/*!
 * \page rtt_format_defined RTT Format File Structure
 * The following example "mesh" documents the format of the RTT file and 
 * explains the associated nomenclature. Graphically depictions of the \ref 
 * rtt_stdcell and \ref rtt_sortcell are provided via the links.
 *
 * \include RTT_Format.defined
 */

/*!
 * \page rtt_stdcell RTT Format Standard Cell Definitions
 * The standard (i.e., default) RTT_Format side set numbering is depicted on 
 * this page. Note that the "right hand rule" is used to return the direction
 * of the outward-directed normal when the nodes are traversed in the order
 * that is specified in the side set node ordering. The RTT standard cell 
 * definitions do not assume any particular orientation relative to the system
 * coordinate system.
 *
 * <center>
 *   <table>
 *     <tr>
 *       <td align=center valign=center>
 *         <img src="../../draco/src/amr_mesh/stdcell.jpg"> 
 *       </td>
 *     </tr>
 *   </table>
 * </center>
 *
 */

/*!
 * \page rtt_sortcell RTT Format Renumbered Cell Definitions for AMR
 * The RTT_Format side set renumbered for AMR is depicted on this page. Note 
 * that the "right hand rule" is used to return the direction of the 
 * outward-directed normal when the nodes are traversed in the order that
 * is specified in the side set node ordering. The cell definitions assume 
 * the orientation relative to the coordinate system that is indicated in the 
 * figures.
 *
 * <center>
 *   <table>
 *     <tr>
 *       <td align=center valign=center>
 *         <img src="../../draco/src/amr_mesh/sortcell.jpg"> 
 *       </td>
 *     </tr>
 *   </table>
 * </center>
 *
 */

//---------------------------------------------------------------------------//
//                              end of RTT_Format.hh
//---------------------------------------------------------------------------//
