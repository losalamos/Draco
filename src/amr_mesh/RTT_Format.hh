//----------------------------------*-C++-*----------------------------------//
// RTT_Format.hh
// Shawn Pautz (TychoMesh.cc original) / B.T. Adams (Extended to RTT_Format.hh)
// 7 June 99
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
using dsxx::SP;

namespace rtt_format
{
 
//===========================================================================//
// class RTT_Format - 
//
// Purpose : A generalized input routine to read an RTT Format mesh file.
//
// revision history:
// -----------------
// 0) original (developed by extending the TychoMesh.hh file of Shawn Pautz)
// 
//===========================================================================//

class RTT_Format 
{

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
	dsxx::Mat1<string> comments;

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
	dsxx::Mat1<int> nnode_flags;
	int nnode_data;

	int nsides;
	int nside_types;
	dsxx::Mat1<int> side_types;
	int nside_flag_types;
	dsxx::Mat1<int> nside_flags;
	int nside_data;

	int ncells;
	int ncell_types;
	dsxx::Mat1<int> cell_types;
	int ncell_flag_types;
	dsxx::Mat1<int> ncell_flags;
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

    class Flags
    {
	int nflags;
	string name;
	dsxx::Mat1<int> flag_nums;
	dsxx::Mat1<string> flag_names;

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

    class NodeFlags
    {
	const Dims & dims;
	dsxx::Mat1<dsxx::SP<Flags> > flagTypes;

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

    class SideFlags
    {
	const Dims & dims;
	dsxx::Mat1<dsxx::SP<Flags> > flagTypes;

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
	    string bc("boundarycditBOUNDARYCDIT");
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
    };

    class CellFlags
    {
	const Dims & dims;
	dsxx::Mat1<dsxx::SP<Flags> > flagTypes;

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
	    string matl("materialMATERIAL");
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
    };

    class NodeDataIDs
    {
	const Dims & dims;
	dsxx::Mat1<string> names;
	dsxx::Mat1<string> units;

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

    class SideDataIDs
    {
	const Dims & dims;
	dsxx::Mat1<string> names;
	dsxx::Mat1<string> units;

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

    class CellDataIDs
    {
	const Dims & dims;
	dsxx::Mat1<string> names;
	dsxx::Mat1<string> units;

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

    class CellDef
    {
	const CellDefs & cellDefs;
	const string name;
	int nnodes;
	int nsides;
	dsxx::Mat1<int> side_types;
	dsxx::Mat1<set<int> > sides;
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

    class CellDefs
    {
	const Dims & dims;
	dsxx::Mat1<dsxx::SP<CellDef> > defs;

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

    class Nodes
    {
	const NodeFlags & nodeFlags;
	const Dims & dims;
	dsxx::Mat2<double> coords;
	dsxx::Mat1<int> parents;
	dsxx::Mat2<int> flags;
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

    class Sides
    {
	const SideFlags & sideFlags;
	const Dims & dims;
	const CellDefs & cellDefs;
        const Nodes & ptrNodes;
	dsxx::Mat1<int> sideType;
	dsxx::Mat2<int> nodes;
	dsxx::Mat2<int> flags;

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
    };

    class Cells
    {
	const CellFlags & cellFlags;
	const Dims & dims;
	const CellDefs & cellDefs;
        const Nodes & ptrNodes;
	dsxx::Mat1<int> cellType;
	dsxx::Mat2<int> nodes;
	dsxx::Mat2<int> flags;

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

        int get_material_flag_number() const 
	{ return cellFlags.get_material_flag_number(); }

    };

    class NodeData
    {
	const Dims & dims;
	dsxx::Mat2<double> data;
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

    class SideData
    {
	const Dims & dims;
	dsxx::Mat2<double> data; 
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

    class CellData
    {
	const Dims & dims;
	dsxx::Mat2<double> data; 
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

    class Connectivity
    {
	const Dims & dims;
	const CellDefs & cellDefs;
	const Cells & cells;
	const Sides & sides;
        const Nodes & nodes;
        vector<vector<vector<int> > > adjCell;
        multimap<int, int> bndryFaces;
      
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

    // Destructors
    ~RTT_Format() {}

    // ACCESSORS
    // header data access
    string get_header_version() const { return header.get_version(); }
    string get_header_title() const { return header.get_title(); }
    string get_header_date() const { return header.get_date(); }
    int get_header_cycle() const { return header.get_cycle(); }
    double get_header_time() const { return header.get_time(); }
    int get_header_ncomments() const { return header.get_ncomments(); }	
    string get_header_comments(int i) const { return header.get_comments(i); }

    // dimensions units and cell definition data access
    string get_dims_coor_units() const { return dims.get_coor_units(); }
    string get_dims_prob_time_units() const 
    { return dims.get_prob_time_units(); }
    int get_dims_ncell_defs() const { return dims.get_ncell_defs(); }
    int get_dims_nnodes_max() const { return dims.get_nnodes_max(); }
    int get_dims_nsides_max() const { return dims.get_nsides_max(); }
    int get_dims_nnodes_side_max() const { return dims.get_nnodes_side_max(); }

    // dimensions node data access
    int get_dims_ndim() const { return dims.get_ndim(); }
    int get_dims_ndim_topo() const { return dims.get_ndim_topo(); }
    int get_dims_nnodes() const { return dims.get_nnodes(); }
    int get_dims_nnode_flag_types() const 
    { return dims.get_nnode_flag_types(); }
    int get_dims_nnode_flags(int i) const { return dims.get_nnode_flags(i); }
    int get_dims_nnode_data() const { return dims.get_nnode_data(); }

    // dimensions side data access
    int get_dims_nsides() const { return dims.get_nsides(); }
    int get_dims_nside_types() const { return dims.get_nside_types(); }
    int get_dims_side_types(int i) const  { return  dims.get_side_types(i); }
    int get_dims_nside_flag_types() const 
    { return dims.get_nside_flag_types(); }
    int get_dims_nside_flags(int i) const { return dims.get_nside_flags(i); }
    int get_dims_nside_data() const { return dims.get_nside_data(); }

    // dimensions cell data access
    int get_dims_ncells() const { return dims.get_ncells(); }
    int get_dims_ncell_types() const { return dims.get_ncell_types(); }    
    int get_dims_cell_types(int i)  const { return dims.get_cell_types(i); }
    int get_dims_ncell_flag_types() const 
    { return dims.get_ncell_flag_types(); }
    int get_dims_ncell_flags(int i) const { return dims.get_ncell_flags(i); }
    int get_dims_ncell_data() const { return dims.get_ncell_data(); }

    // dimensions renumbering flag
    bool get_dims_renumber() const { return dims.get_renumber(); }

    // node flags access
    string get_node_flags_flag_type(int flagtype) const 
	{ return spNodeFlags->get_flag_type(flagtype); }
    int get_node_flags_flag_number(int flagtype, int flag_index) const 
	{ return spNodeFlags->get_flag_number(flagtype, flag_index); }
    int get_node_flags_flag_size(int flagtype) const 
	{ return spNodeFlags->get_flag_size(flagtype); }
    string get_node_flags_flag_name(int flagtype, int flag_index) const 
	{ return spNodeFlags->get_flag_name(flagtype, flag_index); }

    // side flags access
    string get_side_flags_flag_type(int flagtype) const 
	{ return spSideFlags->get_flag_type(flagtype); }
    int get_side_flags_flag_number(int flagtype, int flag_index) const 
	{ return spSideFlags->get_flag_number(flagtype, flag_index); }
    int get_side_flags_flag_size(int flagtype) const 
	{ return spSideFlags->get_flag_size(flagtype); }
    string get_side_flags_flag_name(int flagtype, int flag_index) const 
	{ return spSideFlags->get_flag_name(flagtype, flag_index); }
    int get_side_flags_boundary_flag_number() const 
	{ return spSideFlags->get_boundary_flag_number(); }

    // cell flags access
    string get_cell_flags_flag_type(int flagtype) const 
	{ return spCellFlags->get_flag_type(flagtype); }
    int get_cell_flags_flag_number(int flagtype, int flag_index) const 
	{ return spCellFlags->get_flag_number(flagtype, flag_index); }
    int get_cell_flags_flag_size(int flagtype) const 
	{ return spCellFlags->get_flag_size(flagtype); }
    string get_cell_flags_flag_name(int flagtype, int flag_index) const 
	{ return spCellFlags->get_flag_name( flagtype, flag_index); }
    int get_cell_flags_material_flag_number() const 
	{ return spCellFlags->get_material_flag_number(); }

    // node data ids access
    string get_node_data_id_name(int id_numb) const 
	{ return spNodeDataIds->get_data_id_name(id_numb) ; }
    string get_node_data_id_units(int id_numb) const 
	{ return spNodeDataIds->get_data_id_units(id_numb); }

    // side data ids access
    string get_side_data_id_name(int id_numb) const 
	{ return spSideDataIds->get_data_id_name(id_numb) ; }
    string get_side_data_id_units(int id_numb) const 
	{ return spSideDataIds->get_data_id_units(id_numb); }

    // cell data ids access
    string get_cell_data_id_name(int id_numb) const 
	{ return spCellDataIds->get_data_id_name(id_numb) ; }
    string get_cell_data_id_units(int id_numb) const 
	{ return spCellDataIds->get_data_id_units(id_numb); }

    // cell definitions access
    string get_cell_defs_name(int i) const { return spCellDefs->get_name(i); }
    int get_cell_defs_nnodes(int i) const { return spCellDefs->get_nnodes(i); }
    int get_cell_defs_nsides(int i) const { return spCellDefs->get_nsides(i); }
    int get_cell_defs_side_types(int i, int s) const 
        { return spCellDefs->get_side_types(i,s); }
    const set<int> & get_cell_defs_side(int i, int s) const 
        { return spCellDefs->get_side(i,s); }
    const vector<int> & get_cell_defs_ordered_side(int i, int s) const 
	{ return spCellDefs->get_ordered_side(i,s); }

    // nodes access
    double get_nodes_coords(int node_numb, int coord_index) const
	{ return spNodes->get_coords(node_numb, coord_index); }

    vector<double> get_nodes_coords(int node_numb) const
	{ return spNodes->get_coords(node_numb); }

    int get_nodes_node(vector<double> node_coords) const
        { return spNodes->get_node(node_coords); }

    int get_nodes_parents(int node_numb) const
	{ return spNodes->get_parents(node_numb); }

    int get_nodes_flags(int node_numb, int flag_numb) const
	{ return spNodes->get_flags(node_numb, flag_numb); }

    int get_nodes_map(int node_numb) const
        { return spNodes->get_map(node_numb);}

    // sides access
    int get_sides_type(int side_numb) const
	{ return spSides->get_type(side_numb); }

    int get_sides_nodes(int side_numb,int node_numb) const
	{ return spSides->get_nodes(side_numb, node_numb); }

    int get_sides_flags(int side_numb,int flag_numb) const
	{ return spSides->get_flags(side_numb, flag_numb); }

    // cells access
    int get_cells_type(int cell_numb) const
	{ return spCells->get_type(cell_numb); }

    int get_cells_nodes(int cell_numb,int node_numb) const
	{ return spCells->get_nodes(cell_numb, node_numb); }

    int get_cells_flags(int cell_numb,int flag_numb) const
	{ return spCells->get_flags(cell_numb, flag_numb); }

    // connectivity access
    int get_adjCell(int cell, int face, int adjcell = 0) const
        { return spConnectivity->get_adjCell(cell, face, adjcell); }
    int get_adjCell_size(int cell, int face) const
        { return spConnectivity->get_adjCell_size(cell,face); }
    int get_bndryFaces_count(int face) const
	{ return spConnectivity->get_bndryFaces_count(face); }
    set<int> get_bndryCells(int face) const
        { return spConnectivity->get_bndryCells(face); }
    bool check_bndryFace(int cell, int face) const
        { return spConnectivity->check_bndryFace(cell,face); }
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

//---------------------------------------------------------------------------//
//                              end of RTT_Format.hh
//---------------------------------------------------------------------------//
