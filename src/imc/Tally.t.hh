//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/Tally.t.hh
 * \author Todd J. Urbatsch
 * \date   Mon Apr  6 14:38:03 1998
 * \brief  Tally class template member function definitions.
 * \note   Copyright © 2003 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Tally.hh"
#include "ds++/Assert.hh"
#include <iomanip>

namespace rtt_imc 
{

//---------------------------------------------------------------------------//
// CONSTRUCTORS
//---------------------------------------------------------------------------//
/*!
 * \brief Explicit constructor for Tally class.

 * The constructor builds the Tally with all of its appropriate fields
 * intialized to zero.

 * \param mesh on-processor mesh used for particle transport

 */
template<class MT>
Tally<MT>::Tally(SP_MT mesh)
    : energy_dep(mesh),
      energy_dep_tot(0),
      momentum_dep(mesh),
      eweighted_pathlen(mesh), 
      census_energy(mesh), 
      new_ecen_tot(0), 
      new_ncen(mesh),
      new_ncen_tot(0),
      n_effscat(0), 
      n_thomscat(0), 
      n_killed(0), 
      ew_killed(0),
      n_escaped(0),
      ew_escaped(0), 
      n_bndcross(0),
      n_reflections(0),
      ew_escaped_per_face(2*(mesh->get_Coord().get_dim()))
{
    Require (mesh);
}

//---------------------------------------------------------------------------//
// ACCUMULATE DATA FUNCTIONS
//---------------------------------------------------------------------------//
/*!
 * \brief Accumulate energy in a cell.

 * \param cell cell in which energy is accumulated
 * \param energy energy to accumulate in cell

 */
template<class MT>
void Tally<MT>::deposit_energy(const int cell, const double energy)
{
    Require (cell > 0);
    Require (cell <= energy_dep.get_Mesh().num_cells());

    Check (pos(energy));

    energy_dep(cell) += energy;
    energy_dep_tot   += energy;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Accumulate momentum in a cell.

 * \param cell cell in which momentum is accumulated 
 * \param energy_wt energy-weight of particle

 */
template<class MT>
void Tally<MT>::accumulate_momentum(const int        cell, 
				    const double     energy_wt, 
				    const sf_double &omega)
{
    Require (cell > 0);
    Require (cell <= energy_dep.get_Mesh().num_cells());

    Check (omega.size() >= momentum_dep.size());

    // do momentum deposition only in coordinate dimensions supplied by mesh, 
    // e.g. in an XY mesh we only do momentum deposition in the x and y
    // coordinates even though there is a z dimension in the direction
    // (Omega) vector
    for (int dim = 1; dim <= momentum_dep.size(); dim++)
	momentum_dep(dim, cell) += energy_wt * omega[dim-1];
}

//---------------------------------------------------------------------------//
/*!

 * \brief Accumulate the energy-weighted path length in a cell.

 * The energy-weight of a particle as a function of x is: ew = ew(0) *
 * exp(-sig * x).  Integrating this from 0 to d obtains the energy-weighted
 * path length which yields: ew(0) / sig * (1 - exp(-sig * d).  Since the
 * energy-weight of a particle that has travelled distance d is ew(0) *
 * exp(-sig * d), the energy-weighted path length is simply (delta ew) / sig.

 * \param cell cell in which energy-weighted path length is accumulated

 * \param ewpl energy-weighted path length of particle after streaming a
 * distance d

 */
template<class MT>
void Tally<MT>::accumulate_ewpl(const int cell, const double ewpl)
{
    Require (cell > 0);
    Require (cell <= energy_dep.get_Mesh().num_cells());
    
    Check (pos(ewpl));
    eweighted_pathlen(cell) += ewpl;
}

//---------------------------------------------------------------------------//
/*!

 * \brief Accumulate census statistics in a cell.

 * This function accumulates census data that is needed at the end of an IMC
 * timestep.

 * \param cell cell in which census statistics are accumulated

 * \param new_ecen energy of a particle that goes to census in the cell

 * \param num_new_ncen (default=1) number of census particles that new_ecen
 * describes

 */
template<class MT>
void Tally<MT>::accumulate_cen_info(const int    cell,
				    const double new_ecen,
				    const int    num_new_ncen)
{
    Require (cell > 0);
    Require (cell <= energy_dep.get_Mesh().num_cells());

    Check (pos(new_ecen));
    Check (pos(num_new_ncen));

    census_energy(cell) += new_ecen;
    new_ecen_tot        += new_ecen;
    new_ncen(cell)      += num_new_ncen;
    new_ncen_tot        += num_new_ncen;
}

//---------------------------------------------------------------------------//
// PRINT DIAGNOSTICS
//---------------------------------------------------------------------------//
/*!
 * \brief Print out all tally data.

 * \param out ostream to print data to

 */
template<class MT>
void Tally<MT>::print(std::ostream &out) const
{
    using std::endl;
    using std::setw;
    using std::ios;
    using std::setiosflags;

    out << setw(8) << setiosflags(ios::right) << "Cell"
	<< setw(15) << "energy-dep" << setw(15) << "ewpl" 
	<< setw(15) << "new_ecen" 
	<< setw(15) << "new_ncen" << endl;
    out << "------------------------------------" <<
           "------------------------------------" << endl;

    out.precision(4);
    for (int i = 1; i <= energy_dep.get_Mesh().num_cells(); i++)
	out << setw(8) << i 
	    << setw(15) << setiosflags(ios::scientific) << energy_dep(i) 
	    << setw(15) << setiosflags(ios::scientific) << eweighted_pathlen(i) 
	    << setw(15) << setiosflags(ios::scientific) << census_energy(i) 
	    << setw(15) << setiosflags(ios::scientific) << new_ncen(i) 
	    << endl;

    out << setw(8)  << setiosflags(ios::right) << "Cell"
	<< setw(20) << "momentum deposition"   << endl;
    out << "------------------------------------" 
	<< "------------------------------------" << endl;

    out.precision(4);
    for (int i = 1; i <= momentum_dep.get_Mesh().num_cells(); i++)
    {
	out << setw(8) << i;
	for (int dim = 1; dim <= momentum_dep.size(); dim++)
	    out << setw(15) << setiosflags(ios::scientific) 
		<< momentum_dep(dim, i);
	out << endl;
    }

    out << endl;
    out << "Total new census energy: " 
	<< setw(15) << setiosflags(ios::scientific) << new_ecen_tot
	<< ", new_ncen_tot: " << new_ncen_tot << endl;

    out << endl;
    out << setw(16) << " bnd crossings: " << setw(10) << n_bndcross 
	<< setw(16) << " reflections: "   << setw(10) << n_reflections
	<< setw(16) << " eff scatters: "  << setw(10) << n_effscat
	<< setw(16) << " thomson scats: " << setw(10) << n_thomscat
	<< endl;

    out << setw(16) << " killed: " << setw(10) << n_killed
	<< setw(16) << " ( " << ew_killed << " ) " << endl;

    out << setw(16) << " escaped: " << setw(10) << n_escaped
	<< setw(16) << " ( " << ew_escaped << " ) " << endl;

}

//---------------------------------------------------------------------------//
/*!
 * \brief Print particle transport statistics (over a cycle).

 * This function prints out particle transport statistics including number of
 * boundary crossings, number of reflections, number of effective scatters,
 * number killed, number escaped, and number to census.

 * \param out ostream to print data to

 */
template<class MT>
void Tally<MT>::cycle_print(std::ostream &out) const
{
    using std::endl;
    using std::setw;

    out.precision(4);

    out << endl;
    out << setw(16) << " bnd crossings: " << setw(10) << n_bndcross 
	<< setw(16) << " reflections: "   << setw(10) << n_reflections
	<< setw(16) << " eff scatters: "  << setw(10) << n_effscat
	<< endl;

    out << setw(16) << " killed: " << setw(10) << n_killed
	<< " (ew: " << setw(18) << ew_killed << " ) " << endl;

    out << setw(16) << " escaped: " << setw(10) << n_escaped
	<< " (ew: " << setw(18) << ew_escaped << " ) " << endl;

    out << setw(16) << " new_ncen_tot: " << setw(10) << new_ncen_tot
	<< " (ew: " << setw(18) << new_ecen_tot << " ) " << endl;
}

} // end namespace rtt_imc

//---------------------------------------------------------------------------//
//                              end of imc/Tally.t.hh
//---------------------------------------------------------------------------//
