//----------------------------------*-C++-*----------------------------------//
// OS_Interface.hh
// Thomas M. Evans
// Mon Feb 23 17:22:21 1998
//---------------------------------------------------------------------------//
// @> OS_Interface class header file
//---------------------------------------------------------------------------//

#ifndef __imc_OS_Interface_hh__
#define __imc_OS_Interface_hh__

//===========================================================================//
// class OS_Interface - 
//
// Purpose : parses a Orthogonal structured mesh input file, interface
//           between input file and IMCTEST package
//
// revision history:
// -----------------
//  0) original
//  1)  3-13-98 : changed name to OS_Interface because this class is the
//                interface between our input and the builders.
//  2)   5-6-98 : added parser for source
//  3)  6-16-98 : added material-necessary quantities to end-mat block
//  4)  6-17-98 : added card for input of analytic opacities
// 
//===========================================================================//

#include "imc/Names.hh"
#include "imc/OS_Mesh.hh"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

IMCSPACE

using std::string;
using std::vector;
using std::ifstream;

class OS_Interface 
{
private:
  // input file particulars
    string input_file;
    string coord_system;

  // data required for OS_Mesh generation

  // number of fine cells per coarse cell
    vector< vector<int> > fine_cells;
  // recursive total number of fine_cells per coarse cell
    vector< vector<int> > accum_cells;
  // coarse edges
    OS_Mesh::CCVF_a coarse_edge;
  // fine cell edges
    OS_Mesh::CCVF_a fine_edge;
  // boundary conditions
    vector<string> bnd_cond;

  // data required for Opacity<MT> generation

  // zone map
    vector<int> zone;

  // material zones
    vector<int> mat_zone;
    vector<double> density;
    vector<double> kappa;
    vector<double> temperature;
    vector<double> specific_heat;
    double implicitness;
    string analytic_opacity;

  // data required for Source_Init
    vector<double> evol_ext;
    vector<string> ss_pos;
    vector<double> ss_temp;
    vector<double> rad_temp;
    double delta_t;
    int max_cycle;
    int npmax, npnom;
    double dnpdt;
    string ss_dist;
    int capacity;

  // Parser member functions

  // OS_Mesh parser functions
    void parser_Mesh(ifstream &);
    void parser2D(ifstream &);
    void parser3D(ifstream &);

  // Source member functions
    void parser_Source(ifstream &);
    void zone_source_parser(ifstream &);

  // Opacity parser functions
    void parser_Opacity(ifstream &);
    void zone_mapper();
    void cell_zone(int, int);
    void cell_zone(int, int, int);
    void zone_opacity_parser(ifstream &);

public:
  // constructor
    explicit inline OS_Interface(const string &infile);

  // public Parser member functions
    void parser();
    
  // public copy functions for mesh
    string get_coordinates() const { return coord_system; }
    const vector<string>& get_boundaries() const { return bnd_cond; }
    const OS_Mesh::CCVF_a& get_fine_edge() const { return fine_edge; }
    
  // public copy functions for Opacity<MT>
    const vector<int>& get_zone() const { return zone; }
    const vector<int>& get_mat_zone() const { return mat_zone; }
    const vector<double>& get_density() const { return density; }
    const vector<double>& get_kappa() const { return kappa; }
    const vector<double>& get_temperature() const { return temperature; }
    double get_implicit() const { return implicitness; }
    const vector<double>& get_specific_heat() const { return specific_heat; }

  // public copy functions for Source_Init<MT>
    vector<double> get_evol_ext() const;
    const vector<string>& get_ss_pos() const { return ss_pos; }
    const vector<double>& get_ss_temp() const { return ss_temp; }
    vector<double> get_rad_temp() const;
    double get_delta_t() const { return delta_t; }
    int get_npmax() const { return npmax; }
    int get_npnom() const { return npnom; }
    double get_dnpdt() const { return dnpdt; }
    int get_capacity() const { return capacity; }
    string get_ss_dist() const { return ss_dist; }
    string get_analytic_opacity() const { return analytic_opacity; }
    int get_max_cycle() const { return max_cycle; }
};

//---------------------------------------------------------------------------//
// inline functions for OS_Interface
//---------------------------------------------------------------------------//

inline OS_Interface::OS_Interface(const string &infile)
    : input_file(infile), coord_system(""), fine_cells(0), accum_cells(0), 
      coarse_edge(0), fine_edge(0), bnd_cond(0), zone(0), mat_zone(0), 
      density(0), kappa(0), temperature(0), implicitness(0), 
      analytic_opacity("straight"), specific_heat(0), evol_ext(0), ss_pos(0), 
      ss_temp(0), rad_temp(0), delta_t(0), max_cycle(0), npmax(0), dnpdt(0), 
      ss_dist("none"), capacity(0)
{}

CSPACE

#endif                          // __imc_OS_Interface_hh__

//---------------------------------------------------------------------------//
//                              end of imc/OS_Interface.hh
//---------------------------------------------------------------------------//
