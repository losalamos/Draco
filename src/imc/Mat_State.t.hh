//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/Mat_State.t.hh
 * \author Thomas M. Evans
 * \date   Mon Mar  9 16:06:28 1998
 * \brief  Mat_State class member template function definitions.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Mat_State.hh"

namespace rtt_imc 
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 *
 * \brief Mat_State constructor.
 *
 * \param density_ cell-centered field of material densities (g/cc)
 *
 * \param temp_ cell-centered field of material (electron) temperatures (keV)
 *
 * \param spec_heat_ cell-centered field of specific heat capacities
 * (jks/g/keV)
 *
 */
template<class MT>
Mat_State<MT>::Mat_State(const ccsf_double &density_, 
			 const ccsf_double &temp_,
			 const ccsf_double &spec_heat_)
    : density(density_),
      temperature(temp_), 
      spec_heat(spec_heat_)
{
    Ensure (density.get_Mesh() == temperature.get_Mesh());
    Ensure (density.get_Mesh() == spec_heat.get_Mesh());
}

//---------------------------------------------------------------------------//
// public member functions
//---------------------------------------------------------------------------//
// diagnostic for printing

template<class MT>
void Mat_State<MT>::print(std::ostream &output, int cell) const
{
  // print out material state of cell
    using std::endl;
    using std::setw;
    using std::ios;
    using std::setiosflags;
    
    output.precision(4);
    
    output << setw(8) << cell << setw(15) << setiosflags(ios::scientific)
	   << density(cell) << setw(15) << temperature(cell) << setw(15) 
	   << get_dedt(cell) << setw(15) << spec_heat(cell) << endl;
} 

} // end namespace rtt_imc

//---------------------------------------------------------------------------//
//                              end of Mat_State.t.hh
//---------------------------------------------------------------------------//
