//----------------------------------*-C++-*----------------------------------//
// Tally.t.hh
// Todd J. Urbatsch
// Mon Apr  6 14:38:03 1998
//---------------------------------------------------------------------------//
// @> Tally class implementation file.
//---------------------------------------------------------------------------//

#include "Tally.hh"
#include <iomanip>
#include <vector>

namespace rtt_imc 
{

// namespace objects
using std::setw;
using std::ios;
using std::setiosflags;
using std::endl;
using std::ostream;
using std::vector;

using rtt_dsxx::SP;

//---------------------------------------------------------------------------//
// constructors
//---------------------------------------------------------------------------//
// standard Tally constructor

template<class MT>
Tally<MT>::Tally(SP<MT> mesh)
    : energy_dep(mesh), energy_dep_tot(0), momentum_dep(mesh),
      eweighted_pathlen(mesh), census_energy(mesh), new_ecen_tot(0), 
      new_ncen(mesh), new_ncen_tot(0), evol_net(mesh), n_effscat(0), 
      n_thomscat(0), n_killed(0), ew_killed(0), n_escaped(0), ew_escaped(0), 
      n_bndcross(0), n_reflections(0)
{}

//---------------------------------------------------------------------------//
// Tally constructor that assigns evol_net

template<class MT>
Tally<MT>::Tally(SP<MT> mesh, typename MT::CCSF_double evol_net_)
    : energy_dep(mesh), energy_dep_tot(0), momentum_dep(mesh),
      eweighted_pathlen(mesh), census_energy(mesh), new_ecen_tot(0), 
      new_ncen(mesh), new_ncen_tot(0), evol_net(evol_net_), n_effscat(0), 
      n_thomscat(0), n_killed(0), ew_killed(0), n_escaped(0), ew_escaped(0), 
      n_bndcross(0), n_reflections(0)
{
    Ensure (*mesh == evol_net.get_Mesh());
}

//---------------------------------------------------------------------------//
// public member functions
//---------------------------------------------------------------------------//

template<class MT>
void Tally<MT>::deposit_energy(const int cell, const double energy)
{
    energy_dep(cell) += energy;
    energy_dep_tot   += energy;
}

template<class MT>
void Tally<MT>::accumulate_momentum(const int cell, const double energy_wt, 
				    const vector<double> omega)
{
    Check (omega.size() == momentum_dep.size());

    for (int dim = 1; dim <= omega.size(); dim++)
	momentum_dep(dim, cell) += energy_wt * omega[dim-1];
}


template<class MT>
void Tally<MT>::accumulate_ewpl(const int cell, const double ewpl)
{
    eweighted_pathlen(cell) += ewpl;
}


template<class MT>
void Tally<MT>::accumulate_cen_info(const int cell, const double new_ecen)
{
    census_energy(cell) += new_ecen;
    new_ecen_tot        += new_ecen;
    new_ncen(cell)      += 1;
    new_ncen_tot        += 1;
}

//---------------------------------------------------------------------------//
// print functions
//---------------------------------------------------------------------------//

template<class MT>
void Tally<MT>::print(ostream &out) const
{
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

template<class MT>
void Tally<MT>::cycle_print(ostream &out) const
{
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
//                              end of Tally.t.hh
//---------------------------------------------------------------------------//
