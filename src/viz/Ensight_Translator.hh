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

#ifndef __viz_Ensight_Translator_hh__
#define __viz_Ensight_Translator_hh__

#include "traits/Viz_Traits.hh"

#include <sstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <algorithm>
#include <sys/types.h>
#include <sys/stat.h>

#define ENSIGHT_DIR_MODE S_IRWXU | S_IRWXG | S_IRWXO

namespace rtt_viz
{

//===========================================================================//
/*!
 * \enum FE_Cell_Types
 *
 * \brief Finite Element Cell Type indicators.
 *
 * This enumeration matches Finite Element methods to an index.
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
 */
// revision history:
// -----------------
// 0) original
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
    // Number of Ensight cell types.
    int num_ensight_cell_types;

    // Ensight cell names.
    sf_string ensight_cell_names;

    // Number of vertices for a given Ensight cell type.
    sf_int vrtx_cnt;
    
    // Cell types.
    sf_int cell_type_index;

    // Timestep index.
    int ntime_current;

    // Time values at each Ensight dump.
    sf_double dump_times;

    // Dump counter
    int igrdump_num;

  public:
    // Constructor.
    Ensight_Translator(int = 999);

    // Do an Ensight_Dump.
    template<class ISF, class IVF, class SSF, class FVF>
    void ensight_dump(const std_string &prefix, int icycle, double time,
		      double dt, const std_string &gd_wpath,
		      const IVF &ipar, const ISF &iel_type, 
		      const ISF &rgn_index, const FVF &pt_coor,
		      const FVF &ens_vrtx_data, const FVF &ens_cell_data,
		      const SSF &ens_vdata_names, const SSF &ens_cdata_names, 
		      const ISF &rgn_data, const SSF &rgn_name, 
		      const ISF &iproc); 

};

//---------------------------------------------------------------------------//
// include template definitions so that template functions will be
// automatically instantiated in client code

#include "Ensight_Translator.t.hh"

} // end namespace rtt_viz

#endif                          // __viz_Ensight_Translator_hh__

//---------------------------------------------------------------------------//
//                              end of viz/Ensight_Translator.hh
//---------------------------------------------------------------------------//
