//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   viz/Ensight_Translator.hh
 * \author Thomas M. Evans
 * \date   Fri Jan 21 16:36:10 2000
 * \brief  Ensight_Translator header file.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef rtt_viz_Ensight_Translator_hh
#define rtt_viz_Ensight_Translator_hh

#include "traits/Viz_Traits.hh"
#include "ds++/Assert.hh"
#include "ds++/Check_Strings.hh"

#include <sstream>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <algorithm>
#include <cstring>
#include <errno.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>

#define ENSIGHT_DIR_MODE S_IRWXU | S_IRWXG | S_IRWXO

namespace rtt_viz
{

//===========================================================================//
/*!
 * \enum Ensight_Cell_Types
 *
 * \brief Supported Ensight cell types.
 *
 * Values from this enumeration are assigned (by the user) to each cell in
 * the mesh through the rtt_viz::Ensight_Translator::ensight_dump() function
 * input field argument iel_type.  See Section 2.5 (Figure 2-22) in the
 * Ensight 6.2 user manual for further information on Ensight cell types.  
 *
 * \sa Ensight example page for examples of how this enumeration is used to
 * generate the iel_type field.
 */
//===========================================================================//

enum Ensight_Cell_Types 
{
    point = 0,
    two_node_bar,
    three_node_bar,
    three_node_triangle,
    six_node_triangle,
    four_node_quadrangle,
    eight_node_quadrangle,
    four_node_tetrahedron,
    ten_node_tetrahedron,
    five_node_pyramid,
    thirteen_node_pyramid,
    eight_node_hexahedron,
    twenty_node_hexahedron,
    six_node_wedge,
    fifteen_node_wedge
};
 
//===========================================================================//
/*!
 * \class Ensight_Translator
 *
 * \brief A translator for dumping problem data in EnSight format.
 *
 * This class is used to create data dumps that are readable by EnSight.
 * Currently, only ASCII format dumps are supported.  The present incantation
 * of this class is a C++ version of the ensight_dump model from Dante.
 *
 * \anchor Ensight_Translator_strings
 *
 * A note about names: all strings containing names, be they problem names or
 * data names, should not contain the following characters:
 *
 * \arg ( ) [ ] + - \@ space ! \# * ^ \$ /
 *
 * Names in the ensight translation are often used to create directories, so
 * they must not have these characters in them.  Additionally, keep all names
 * under 19 characters.  Finally, all names of data should be unique.
 *
 * \anchor Ensight_Translator_description
 *
 * Ensight_Translator is not a templated class; however, the member functions
 * are templated on field types.  These types and there restrictions are as
 * follows:
 *
 * \arg \b SSF (string scalar field) can be any field type that supports
 * random access iterators and operator[] indexing
 *
 * \arg \b ISF (integer scalar field) can be any field type that supports
 * random access iterators and operator[] indexing
 *
 * \arg \b IVF (integer vector field) uses the rtt_traits::Viz_Traits class
 * to convert to i,j indexing to operator(int, int).  See the
 * rtt_traits::Viz_Traits class to determine which 2D vector containers are
 * specialized
 *
 * \arg \b FVF (floating point vector field) uses the rtt_traits::Viz_Traits
 * class to convert to i,j indexing to operator(int, int).  See the
 * rtt_traits::Viz_Traits class to determine which 2D vector containers are
 * specialized
 *
 * The vector-field (VF) types will most commonly by rtt_dsxx::Mat2 or
 * std::vector<std::vector> types.  Other types will require specialization
 * in the rtt_traits::Viz_Traits class.
 *
 * Data dumped to Ensight will be stored in the directory prefix argument
 * given to the constructor (Ensight_Translator::Ensight_Translator()).  The
 * directory created will have the name "prefix"_ensight. The directory is
 * placed in the path specified by the gd_wpath argument to the constructor.
 * This directory must exist.
 *
 * To launch ensight: select the "prefix".case file that resides in the
 * top-level ensight dump directory from the "file/Data (reader)" menu.  Set
 * the data file "Format" option in Ensight to "case" and hit the "(Set)
 * Geometry" button.  From there see the Ensight manual.  
 */
/*!
 * \example viz/test/tstEnsight_Translator.cc
 *
 * \anchor ensight_test
 *
 * Example usage of the Ensight_Translator class. A two part, 3-D, nine cell
 * mesh graphics dump is performed.  Two vertex data field and two cell data
 * fields are dumped.  After running the test look in testproblem_ensight to
 * see the output.
 */
// revision history:
// -----------------
// 0) original
// 1) 05-APR-03 : removed constructors where the dump_times are arguments;
//                these were confusing and they were causing trouble on the
//                IBM platforms because of a bug in xlC
// 
//===========================================================================//

class Ensight_Translator 
{
  public:
    // Ensight_Translator typedefs.
    typedef std::vector<int>           sf_int;
    typedef std::vector<double>        sf_double;
    typedef std::vector<std::string>   sf_string;
    typedef std::string                std_string;

  private:
    // >>> DATA

    // if true, geometry is static
    bool static_geom;

    // Number of Ensight cell types.
    int num_ensight_cell_types;

    // Ensight cell names.
    sf_string ensight_cell_names;

    // Number of vertices for a given Ensight cell type.
    sf_int vrtx_cnt;
    
    // Cell types.
    sf_int cell_type_index;

    // Vector of dump_times.
    sf_double dump_times;

    // Ensight file prefixes.
    std_string ens_prefix;

    // Names of vertex data.
    sf_string ens_vdata_names;
    
    // Names of cell data.
    sf_string ens_cdata_names;

    // Name of case file.
    std_string case_filename;

    // Name of geometry directory.
    std_string geo_dir;

    // Names of vdata directories.
    sf_string vdata_dirs;

    // Names of cdata directories.
    sf_string cdata_dirs;

  private:
    // >>> PRIVATE IMPLEMENTATION

    // Creates some of the file prefixes and names.
    void createFilenames(const std_string &prefix,
			 const std_string &gd_wpath);
    
    // Write out case file.
    void ensight_case(const double);

    // Write out geometry file.
    template<class IVF, class FVF>
    void ensight_geom(const std_string &, const int, 
		      const double, const double, 
		      const rtt_traits::Viz_Traits<IVF> &,
		      const rtt_traits::Viz_Traits<FVF> &,
		      const sf_string &, const sf_int &, const sf_int &); 

    // Write out vertex data.
    template<class FVF>
    void ensight_vrtx_data(const std_string &, 
			   const rtt_traits::Viz_Traits<FVF> &);

    // Write out cell data.
    template<class FVF>
    void ensight_cell_data(const std_string &,
			   const rtt_traits::Viz_Traits<FVF> &,
			   const sf_int &,
			   const sf_int &, const sf_string &);

    // Initializer used by constructors
    void initialize(const bool graphics_continue);


  public:
    // Constructor.
    template<class SSF>
    Ensight_Translator(const std_string &prefix, const std_string &gd_wpath,
		       const SSF &ens_vdata_names, 
		       const SSF &ens_cdata_names, 
		       const bool overwrite = false,
		       const bool static_geom_in = false); 

    // Do an Ensight_Dump.
    template<class ISF, class IVF, class SSF, class FVF>
    void ensight_dump(int icycle, double time, double dt,
		      const IVF &ipar_in, const ISF &iel_type, 
		      const ISF &cell_rgn_index, const FVF &pt_coor_in,
		      const FVF &ens_vrtx_data_in, 
		      const FVF &ens_cell_data_in, const ISF &rgn_numbers, 
		      const SSF &rgn_name); 

    // >>> ACCESSORS

    //! Get the list of dump times.
    const sf_double& get_dump_times() const { return dump_times; }
};

} // end namespace rtt_viz

//---------------------------------------------------------------------------//
// include template definitions so that template functions will be
// automatically instantiated in client code

#include "Ensight_Translator.t.hh"

#endif                          // rtt_viz_Ensight_Translator_hh

//---------------------------------------------------------------------------//
//                              end of viz/Ensight_Translator.hh
//---------------------------------------------------------------------------//
