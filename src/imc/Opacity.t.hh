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

//===========================================================================//
// OPACITY<MT,GRAY_FREQUENCY> DEFINITIONS
//===========================================================================//

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 *
 * \brief Opacity constructor for Gray_Frequency.
 *
 * \param sigma_abs cell-centered field of absorption opacities
 *
 * \param sigma_thomson cell-centered field of scattering cross sections
 *
 * \param fleck cell-centered field of Fleck factors
 *
 */
template<class MT>
Opacity<MT,Gray_Frequency>::Opacity(SP_Frequency       freq,
				    const ccsf_double &sigma_abs_,
				    const ccsf_double &sigma_thomson_,
				    const ccsf_double &fleck_)
    : frequency(freq),
      sigma_abs(sigma_abs_),
      sigma_thomson(sigma_thomson_), 
      fleck(fleck_)
{
    int num_cells = sigma_abs.get_Mesh().num_cells();
    
    Ensure (sigma_abs.size() == num_cells);
    Ensure (sigma_thomson.size() == num_cells);
    Ensure (fleck.size() == num_cells);
    Ensure (frequency);
}

//---------------------------------------------------------------------------//
// PUBLIC INTERFACE
//---------------------------------------------------------------------------//
/*!
 * \brief Diagnostic print function for Opacity<MT,Gray_Frequency>
 * specialization.
 */
template<class MT>
void Opacity<MT,Gray_Frequency>::print(std::ostream &output) const
{
    // print out opacities in cell
    using std::endl;
    using std::setw;
    using std::ios;
    using std::setiosflags;

    output.precision(4);
    output.setf(ios::scientific, ios::floatfield);

    output << setw(5) << "Cell" << setw(20) << "Gray Absorption"
	   << setw(20) << "Gray Scattering" << endl;
    output << "---------------------------------------------" << endl;

    for (int cell = 1; cell <= num_cells(); cell++)
	output << setw(5) << cell << setw(20) << sigma_abs(cell) 
	       << setw(20) << sigma_thomson(cell); 
}

//===========================================================================//
// OPACITY<MT,MULTIGROUP_FREQUENCY> DEFINITIONS
//===========================================================================//

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 *
 * \brief Opacity constructor for Multigroup_Frequency.
 *
 * \param sigma_abs cell-centered field of vectors of absorption opacities
 *
 * \param sigma_thomson cell-centered field of vectors of scattering cross
 * sections
 *
 * \param fleck cell-centered field of Fleck factors
 *
 * \param int_planck cell-centered field of integrated normalized Planck
 * function
 *
 * \param group_cdf cell-centered-field of the group-wise Planck * sigma
 * absorption cumulative distribution function
 *
 */
template<class MT>
Opacity<MT,Multigroup_Frequency>::Opacity(SP_Frequency freq,
					  const ccsf_vector &sigma_abs_,
					  const ccsf_vector &sigma_thomson_,
					  const ccsf_double &fleck_,
					  const ccsf_double &int_planck,
					  const ccsf_vector &group_cdf)
    : frequency(freq),
      sigma_abs(sigma_abs_),
      sigma_thomson(sigma_thomson_), 
      fleck(fleck_),
      integrated_norm_planck(int_planck),
      emission_group_cdf(group_cdf)
{
    int num_cells = sigma_abs.get_Mesh().num_cells();
    
    Ensure (sigma_abs.size() == num_cells);
    Ensure (sigma_thomson.size() == num_cells);
    Ensure (fleck.size() == num_cells);
    Ensure (integrated_norm_planck.size() == num_cells);
    Ensure (emission_group_cdf.size() == num_cells);
    Ensure (frequency);
}

//---------------------------------------------------------------------------//
// PUBLIC INTERFACE
//---------------------------------------------------------------------------//
/*!
 * \brief Diagnostic print function for Opacity<MT, Multigroup_Frequency>
 * specialization.
 */
template<class MT>
void Opacity<MT,Multigroup_Frequency>::print(std::ostream &output) const
{
    // print out opacities in cell
    using std::endl;
    using std::setw;
    using std::ios;
    using std::setiosflags;

    output.precision(4);
    output.setf(ios::scientific, ios::floatfield);

    for (int cell = 1; cell <= num_cells(); cell++)
    {
	output << "Group opacities in cell " << setw(5) << cell << endl;
	output << "------------------------------" << endl;

	for (int j = 0; j < frequency->get_num_groups(); j++)
	    output << setw(5) << j+1 << setw(20) << sigma_abs(cell)[j] 
		   << setw(20) << sigma_thomson(cell)[j] << endl;

	output << endl;
    }
}

} // end namespace rtt_imc

//---------------------------------------------------------------------------//
//                              end of imc/Opacity.t.hh
//---------------------------------------------------------------------------//
