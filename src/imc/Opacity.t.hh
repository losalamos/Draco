//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/Opacity.t.hh
 * \author Thomas M. Evans
 * \date   Fri Feb  6 13:52:29 1998
 * \brief  Opacity class member template function definitions.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Opacity.hh"
#include <iomanip>

namespace rtt_imc 
{

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
  
 * \brief Opacity constructor.

 * \param sigma_abs cell-centered field of absorption opacities
 * \param sigma_thomson cell-centered field of scattering cross sections
 * \param planck cell-centered field of Planck opacities
 * \param fleck cell-centered field of Fleck factors

 */
template<class MT>
inline Opacity<MT>::Opacity(const ccsf_double &sigma_abs_,
			    const ccsf_double &sigma_thomson_,
			    const ccsf_double &planck_,
			    const ccsf_double &fleck_)
    : sigma_abs(sigma_abs_),
      sigma_thomson(sigma_thomson_),
      planck(planck_), 
      fleck(fleck_)
{
    int num_cells = sigma_abs.get_Mesh().num_cells();
    
    Ensure (sigma_abs.size() == num_cells);
    Ensure (sigma_thomson.size() == num_cells);
    Ensure (planck.size() == num_cells);
    Ensure (fleck.size() == num_cells);
}

//---------------------------------------------------------------------------//
// PUBLIC INTERFACE
//---------------------------------------------------------------------------//
/*!
 * \brief Diagnostic print function.
 */
template<class MT>
void Opacity<MT>::print(std::ostream &output, int cell) const
{
  // print out opacities in cell
    using std::endl;
    using std::setw;
    using std::ios;
    using std::setiosflags;

    output.precision(4);

    output << setw(8) << cell << setw(15) << setiosflags(ios::scientific)
	   << sigma_abs(cell) << setw(15) << sigma_thomson(cell) << endl;
}

} // end namespace rtt_imc

//---------------------------------------------------------------------------//
//                              end of imc/Opacity.t.hh
//---------------------------------------------------------------------------//
