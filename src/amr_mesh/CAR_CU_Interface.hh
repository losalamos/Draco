//----------------------------------*-C++-*----------------------------------//
// CAR_CU_Interface.hh
// B.T. Adams (bta@lanl.gov)
// 18 May 99
//---------------------------------------------------------------------------//
// @> CAR_CU_Interface class header file
//---------------------------------------------------------------------------//

#ifndef __milagro_CAR_CU_Interface_hh__
#define __milagro_CAR_CU_Interface_hh__

//===========================================================================//
// class CAR_CU_Interface - 
//
// Purpose : parses a continuous adaptive refinement cartesian unstructured 
//           mesh input file, interface between input file and IMCTEST package
//
// revision history:
// -----------------
//  0) original (developed from OS_Interface.hh)
// 
//===========================================================================//

#include "CAR_CU_Mesh.hh"
#include "imc/Global.hh"
#include "ds++/SP.hh"
#include "RTT_Format.hh"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

namespace rtt_imc 
{

// stl components
using std::string;
using std::vector;
using std::ifstream;

// draco components
using dsxx::SP;

using rtt_format::RTT_Format;

class CAR_CU_Interface 
{
  private:
    // typenames and typedefs
    typedef rtt_mc::CAR_CU_Mesh CAR_CU_Mesh;

  private:
    // input file particulars
    string input_file;
    string coord_system;

    // data required for CAR_CU_Mesh generation
    string surface_file;
    string rtt_mesh_file;

    // I need (want) the verbose switch
    bool verbose;

    // data required for Opacity<MT> generation

    // zone map
    vector<int> zone;

    // material zones
    vector<int> mat_zone;
    vector<double> density;
    vector<double> kappa;
    vector<double> kappa_thomson;
    vector<double> temperature;
    vector<double> specific_heat;
    double implicitness;
    string analytic_opacity;
    string analytic_sp_heat;

    // data required for Source_Init
    vector<double> evol_ext;
    vector<double> cell_evol;
    vector<double> rad_source;
    vector<double> cell_rsrc;
    double rad_s_tend;
    vector<string> ss_pos;
    vector<double> ss_temp;
    vector< vector<int> > defined_surcells;
    vector<double> rad_temp;
    double delta_t;
    int max_cycle;
    int npmax, npnom;
    double dnpdt;
    string ss_dist;
    int capacity;
    int print_f;
    int buffer;
    int seed;

    // Parser member functions

    // CAR_CU_Mesh parser functions
    SP<RTT_Format> parser_Mesh(ifstream &);

    // Source member functions
    void parser_Source(ifstream &);
    void zone_source_parser(ifstream &);

    // Opacity parser functions
    void parser_Opacity(ifstream &);

  public:
    // constructor
    explicit inline CAR_CU_Interface(const string & infile, 
				     const bool & verbose_);

    ~CAR_CU_Interface() {}

    // public Parser member functions
    SP<RTT_Format> parser();
    
    // public copy functions for mesh
    string get_coordinates() const { return coord_system; }

    // public copy functions for Opacity<MT>
    vector<double> get_density() const;
    vector<double> get_kappa() const;
    vector<double> get_kappa_thomson() const;
    vector<double> get_specific_heat() const;
    vector<double> get_temperature() const;
    string get_analytic_opacity() const { return analytic_opacity; }
    string get_analytic_sp_heat() const { return analytic_sp_heat; }

    // accessor function to get implicitness factor (Fleck's alpha)
    double get_implicit() const { return implicitness; }

    // public copy functions for Source_Init<MT>
    int get_zone_size()  { return zone.size();}
    vector<double> get_evol_ext() const { return cell_evol; }
    double get_evol_ext(int cell) const { return cell_evol[cell -1]; }
    vector<double> get_rad_source() const { return cell_rsrc; }
    double get_rad_source(int cell) const { return cell_rsrc[cell -1]; }
    double get_rad_s_tend() const { return rad_s_tend; }
    // return the number of grouped surface source cell sets
    int get_ss_pos_size() { return ss_pos.size(); }
    // return the position (lox, hix, etc.) of a set of grouped surface source 
    // cells
    string get_ss_pos(int surface) { return ss_pos[surface - 1]; }
    // return the positions (lox, hix, etc.) of the all of the grouped surface 
    // source cells
    const vector<string> & get_ss_pos() const { return ss_pos; }
    // return the temperature of a set of the grouped surface source cells
    const double & get_ss_temp(int surface) const 
    { return ss_temp[surface - 1]; }
    // return the temperature of the all of the grouped surface source cells
    const vector<double> & get_ss_temp() const { return ss_temp; }
    // return the number of grouped surface source cells in a given set
    int get_ss_cells_size(int surface) 
    { return defined_surcells[surface - 1].size(); }
    // return the defined surface source cells in a given set
    vector<int> get_defined_surcells(int surface) 
    {
        vector<int> source_set(defined_surcells[surface - 1].size());
	for (int cell = 0; cell < defined_surcells[surface - 1].size(); cell++)
	    source_set[cell] = defined_surcells[surface - 1][cell];

        return source_set; 
    }
    // return all of the defined surface cell sets.
    const vector< vector<int> >& get_defined_surcells() const {
	return defined_surcells; } 
    vector<double> get_rad_temp() const;
    double get_delta_t() const { return delta_t; }
    int get_npmax() const { return npmax; }
    int get_npnom() const { return npnom; }
    double get_dnpdt() const { return dnpdt; }
    int get_capacity() const { return capacity; }
    string get_ss_dist() const { return ss_dist; }
    int get_max_cycle() const { return max_cycle; }
    int get_printf() const { return print_f; }
    int get_buffer() const { return buffer; }
    int get_seed() const { return seed; }

    // public functions required for CAR_CU_Mesh generation
    string get_surface_file() const { return surface_file; }
    string get_mesh_file() const { return rtt_mesh_file; }

};


//---------------------------------------------------------------------------//
// inline functions for CAR_CU_Interface
//---------------------------------------------------------------------------//

inline CAR_CU_Interface::CAR_CU_Interface(const string &infile, 
					  const bool & verbose_)
      : input_file(infile), verbose(verbose_), coord_system(""), 
        zone(0), mat_zone(0), density(0), kappa(0), 
        kappa_thomson(0), temperature(0), implicitness(0), 
        analytic_opacity("straight"), analytic_sp_heat("straight"), 
        specific_heat(0), evol_ext(0), rad_source(0), rad_s_tend(0), 
        ss_pos(0), ss_temp(0), rad_temp(0), delta_t(0), max_cycle(0), 
        npmax(0), dnpdt(0), ss_dist("none"), capacity(0), print_f(1), 
        buffer(1000), seed(9836592), surface_file("undefined"),
        rtt_mesh_file("undefined")
{}

} // end namespace rtt_imc

#endif                          // __milagro_CAR_CU_Interface_hh__

//---------------------------------------------------------------------------//
//                              end of milagro/CAR_CU_Interface.hh
//---------------------------------------------------------------------------//
