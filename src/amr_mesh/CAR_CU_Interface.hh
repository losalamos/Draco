//----------------------------------*-C++-*----------------------------------//
// CAR_CU_Interface.hh
// B.T. Adams (bta@lanl.gov)
// 18 May 99
/*! 
 * \file   amr_mesh/CAR_CU_Interface.hh
 * \author B.T. Adams
 * \date   Tue May 18 10:33:26 1999
 * \brief  Header file for CAR_CU_Interface class library.
 */
//---------------------------------------------------------------------------//
// @> CAR_CU_Interface class header file
//---------------------------------------------------------------------------//

#ifndef __amr_CAR_CU_Interface_hh__
#define __amr_CAR_CU_Interface_hh__

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

namespace rtt_amr 
{

// stl components
using std::string;
using std::vector;
using std::ifstream;

// draco components
using dsxx::SP;

using rtt_format::RTT_Format;

/*!
 * \brief  The continuous adaptive refinement (CAR) Cartesion unstructured 
 *         (CU) Interface class provides member functions that can be used to
 *         both read a user-input file and parse a mesh file in the \ref 
 *         rtt_format_defined by constructing an rtt_format::RTT_Format class 
 *         object. The CAR_CU_Interface serves as a "preprocessor" to be 
 *         implemented before the CAR_CU_Builder is instantiated and used to
 *         construct a CAR_CU_Mesh. A smart pointer to a CAR_CU_Interface 
 *         class object is also a required argument to the constructors for 
 *         the rtt_imc::Opacity_Builder and rtt_imc::Source_Init classes. 
 *         Finally, this class also contains the data needed to construct a
 *         rtt_imc::Mat_State class object.
 *
 *\sa The \ref amr_overview presents a summary of the capabilities and the
 *    intended usage of this mesh class. A \ref rtt_amr_input is also provided.
 */     
class CAR_CU_Interface 
{
  private:
    // typenames and typedefs
    typedef rtt_amr::CAR_CU_Mesh CAR_CU_Mesh;

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
/*!
 * \brief Private member function that parses the title-block of the 
 *        user-input file and the entire RTT_Format mesh file that is 
 *        designated therein. The latter is accomplished by constructing 
 *        an rtt_format::RTT_Format class object. The constructor for the
 *        RTT_Format class automatically invokes member functions to both 
 *        parse the mesh file and determine the mesh connectivity.
 * \param infile An ifstream class object that is attached to the user-input 
 *               file.    
 * \return Smart pointer to the new RTT_Format class object.
 */
SP<RTT_Format> parser_Mesh(ifstream & infile);

    // Opacity parser functions
/*!
 * \brief Private member function that parses the material information block
 *        of the user-input file and assigns the material and opacitity 
 *        variables.
 * \param infile An ifstream class object that is attached to the user-input 
 *               file.    
 */
    void parser_Opacity(ifstream & infile);

    // Source member functions
/*!
 * \brief Private member function that controls the parsing of the source block
 *        of the user-input file.
 * \param infile An ifstream class object that is attached to the user-input 
 *               file.    
 */
    void parser_Source(ifstream & infile);
/*!
 * \brief Private member function that parses the source block of the 
 *        user-input file and assigns the source variables.
 * \param infile An ifstream class object that is attached to the user-input
 *               file.    
 */
    void zone_source_parser(ifstream & infile);

  public:
    // constructor
/*!
 * \brief Constructs a CAR_CU_Mesh Interface class object.
 * \param infile User-input file name (must contain the RTT_format mesh file 
 *               name).
 * \param verbose Switch used to turn detailed run-time reporting on/off.
 */
    explicit inline CAR_CU_Interface(const string & infile, 
				     const bool & verbose);

/*!
 * \brief Destroys a CAR_CU_Mesh Interface class object.
 */
    ~CAR_CU_Interface() {}

    // public Parser member functions
/*!
 * \brief Controls the parsing of both the user-input and RTT_Format mesh 
 *        files. The latter is accomplished via a call to the parser_Mesh 
 *        private member function, which constructs an rtt_Format::RTT_Format
 *        class object. The constructor for the RTT_Format class automatically 
 *        invokes member functions to both parse the mesh file and determine 
 *        the mesh connectivity.
 * \return Smart pointer to the new RTT_Format class object.
 */
    SP<RTT_Format> parser();
    
    // public copy functions for mesh
/*!
 * \brief Returns the problem coordinate system (e.g., xyz).
 * \return Coordinate system.
 */
    string get_coordinates() const { return coord_system; }

    // public copy functions for Opacity<MT>
/*!
 * \brief Returns the cell-based density array.
 * \return Cell densities.
 */
    vector<double> get_density() const;
/*!
 * \brief Returns the cell-based kappa scattering cross section array.
 * \return Cell kappa scattering cross sections.
 */
    vector<double> get_kappa() const;
/*!
 * \brief Returns the cell-based kappa thomson scattering cross section array.
 * \return Cell kappa thomson scattering cross sections.
 */
    vector<double> get_kappa_thomson() const;
/*!
 * \brief Returns the cell-based specific heat array.
 * \return Cell specific heats.
 */
    vector<double> get_specific_heat() const;
/*!
 * \brief Returns the cell-based temperature array.
 * \return Cell temperatures.
 */
    vector<double> get_temperature() const;
/*!
 * \brief Returns the analytic opacity.
 * \return Analytic opacity.
 */
    string get_analytic_opacity() const { return analytic_opacity; }
/*!
 * \brief Returns the analytic specific heat.
 * \return Analytic specific heat.
 */
    string get_analytic_sp_heat() const { return analytic_sp_heat; }

    // accessor function to get implicitness factor (Fleck's alpha)
/*!
 * \brief Returns the implicitness factor (Fleck's alpha).
 * \return Implicitness factor.
 */
    double get_implicit() const { return implicitness; }

    // public copy functions for Source_Init<MT>
/*!
 * \brief Returns the size of the material zone vector (same as the number of
 *        cells in the mesh).
 * \return Number of cells in the mesh.
 */
    int get_zone_size()  { return zone.size();}
/*!
 * \brief Returns the cell-based initial radiation temperature array.
 * \return Cell initial radiation temperatures.
 */
    vector<double> get_rad_temp() const;
/*!
 * \brief Returns the cell-based volumetric source array.
 * \return Cell volumetric sources.
 */
    vector<double> get_evol_ext() const { return cell_evol; }
 /*!
 * \brief Returns the volumetric source for the specified cell.
 * \param cell Cell number.
 * \return Cell volumetric source.
 */
    double get_evol_ext(int cell) const { return cell_evol[cell -1]; }
/*!
 * \brief Returns the cell-based external radiation source array.
 * \return Cell external radiation sources.
 */
    vector<double> get_rad_source() const { return cell_rsrc; }
 /*!
 * \brief Returns the external radiation source for the specified cell.
 * \param cell Cell number.
 * \return Cell external radiation source.
 */
    double get_rad_source(int cell) const { return cell_rsrc[cell -1]; }
/*!
 * \brief Returns the cut-off time for the external radiation sources.
 * \return The external radiation source cut-off time.
 */
    double get_rad_s_tend() const { return rad_s_tend; }
    // return the number of grouped surface source cell sets
/*!
 * \brief Returns the number of a grouped surface source cell sets.
 * \return The number of grouped surface source cell sets.
 */
    int get_ss_pos_size() { return ss_pos.size(); }
    // return the position (lox, hix, etc.) of a set of grouped surface source 
    // cells
/*!
 * \brief Returns the position (either lox, loy, loz, hix, hiy, or hiz) of 
 *        the specified set of grouped surface source cells.
 * \param surface Surface source cells set number.
 * \return Surface source cells set position.
 */
    string get_ss_pos(int surface) { return ss_pos[surface - 1]; }
    // return the positions (lox, hix, etc.) of the all of the grouped surface 
    // source cells
/*!
 * \brief Returns the position (either lox, loy, loz, hix, hiy, or hiz) of 
 *        all of the specified sets of grouped surface source cells.
 * \return Surface source cells set positions.
 */
    const vector<string> & get_ss_pos() const { return ss_pos; }
    // return the temperature of a set of the grouped surface source cells
/*!
 * \brief Returns the temperature of the specified set of grouped surface 
 *        source cells.
 * \param surface Surface source cells set number.
 * \return Surface source cells set temperature.
 */
    const double & get_ss_temp(int surface) const 
    { return ss_temp[surface - 1]; }
    // return the temperature of the all of the grouped surface source cells
/*!
 * \brief Returns the temperature of all of the grouped surface source cells
 *        sets.
 * \return Surface source cells set temperatures.
 */
    const vector<double> & get_ss_temp() const { return ss_temp; }
    // return the number of grouped surface source cells in a given set
/*!
 * \brief Returns the number of grouped surface source cells in the specified
 *        set.
 * \param surface Surface source cells set number.
 * \return The number of grouped surface source cells in the set.
 */
    int get_ss_cells_size(int surface) 
    { return defined_surcells[surface - 1].size(); }
    // return the defined surface source cells in a given set
/*!
 * \brief Returns the surface source cells in the specified set.
 * \param surface Surface source cells set number.
 * \return Surface source cells in the set.
 */
    vector<int> get_defined_surcells(int surface) 
    {
        vector<int> source_set(defined_surcells[surface - 1].size());
	for (int cell = 0; cell < defined_surcells[surface - 1].size(); cell++)
	    source_set[cell] = defined_surcells[surface - 1][cell];

        return source_set; 
    }
    // return all of the defined surface cell sets.
/*!
 * \brief Returns all of the surface source cell sets.
 * \return All surface source cell sets.
 */
    const vector< vector<int> >& get_defined_surcells() const {
	return defined_surcells; } 
/*!
 * \brief Returns the initial time step size.
 * \return The initial time step size.
 */
    double get_delta_t() const { return delta_t; }
/*!
 * \brief Returns the maximum number of particles (for monte carlo).
 * \return The maximum number of particles.
 */
    int get_npmax() const { return npmax; }
/*!
 * \brief Returns the nominal number of particles (for monte carlo).
 * \return The nominal number of particles.
 */
    int get_npnom() const { return npnom; }
/*!
 * \brief Returns the rate of change in the number of particles (for monte 
 *        carlo).
 * \return The rate of change in the number of particles.
 */
    double get_dnpdt() const { return dnpdt; }
/*!
 * \brief Returns the number of cells each processor will hold for parallel
 *        execution.
 * \return The processor cell capacity.
 */
    int get_capacity() const { return capacity; }
/*!
 * \brief Returns the surface source angular distribution (e.g., cosine).
 * \return The surface source angular distribution.
 */
    string get_ss_dist() const { return ss_dist; }
/*!
 * \brief Returns the maximum number of cycles to run.
 * \return The maximum number of cycles.
 */
    int get_max_cycle() const { return max_cycle; }
/*!
 * \brief Returns the number of cycles between print outs.
 * \return The number of cycles between print outs.
 */
    int get_printf() const { return print_f; }
/*!
 * \brief Returns the buffer size in particles for in-cycle communication 
 *        between  processors for parallel execution.
 * \return The particle buffer size.
 */
    int get_buffer() const { return buffer; }
/*!
 * \brief Returns the input random number seed (for monte carlo).
 * \return The random number seed.
 */
    int get_seed() const { return seed; }

    // public functions required for CAR_CU_Mesh generation
/*!
 * \brief Returns the file name, read from the user-input file, that contains
 *        the analytic surface input data (optional input).
 * \return The analytic surface file name.
 */
    string get_surface_file() const { return surface_file; }
/*!
 * \brief Returns the file name, read from the user-input file, that contains
 *        the RTT_Format mesh file data (required input).
 * \return The RTT_Format mesh file name.
 */
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

} // end namespace rtt_amr

#endif                          // __amr_CAR_CU_Interface_hh__

/*!
 * \page rtt_amr_input Sample AMR_Mesh User-Input File
 * The following example input deck documents the format of the amr_mesh 
 * user-input file and explains the associated nomenclature.
 *
 * \include amr_mesh.inp
 */


//---------------------------------------------------------------------------//
//                              end of amr_mesh/CAR_CU_Interface.hh
//---------------------------------------------------------------------------//
