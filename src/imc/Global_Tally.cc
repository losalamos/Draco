//----------------------------------*-C++-*----------------------------------//
// Global_Tally.cc
// Thomas M. Evans
// Wed Jun 17 10:21:19 1998
//---------------------------------------------------------------------------//
// @> Global_Tally class implementation file
//---------------------------------------------------------------------------//

#include "imc/Global_Tally.hh"
#include <iomanip>

IMCSPACE

using std::setiosflags;
using std::ios;
using std::endl;
using std::setw;
using std::fabs;

//---------------------------------------------------------------------------//
// constructor
//---------------------------------------------------------------------------//

template<class MT, class PT>
Global_Tally<MT,PT>::Global_Tally(const MT &mesh, const Mat_State<MT>
				  &material, const Source_Init<MT> &source)
    : temperature(mesh.num_cells()), dedt(mesh.num_cells()), e_elec_tot(0),
      eint_initial(0), eint_begin(0), eint_end(0), e_in_probtot(0),
      e_esc_probtot(0), e_loss_probtot(0), evol_net(mesh.num_cells()),   
      evoltot(0), evolext(0), nvoltot(0), ncen(mesh.num_cells()), ncentot(0),  
      ncenrun(0), eradtot_b(0), eradtot_e(0), ecen(mesh.num_cells()),
      ecentot(0), ewpl(mesh.num_cells()), esstot(0), nsstot(0), eloss_vol(0), 
      eloss_ss(0), eloss_cen(0), e_escape(0)
{
    Require (mesh.num_cells() == material.num_cells());

  // assign the census and cell temperature data
    int ncentot = 0;
    for (int cell = 1; cell <= mesh.num_cells(); cell++)
    {
	temperature[cell-1] = material.get_T(cell);
	dedt[cell-1]        = material.get_dedt(cell);
	ncen[cell-1]        = source.get_ncen(cell);
	ncentot            += ncen[cell-1];
    }
    census    = source.get_census();
    eradtot_e = source.get_ecentot();

  // do the initial internal energy calculation
    set_energy_begin(source);
    eint_initial = eint_begin;

    Ensure (census);
    Ensure (census->size() == ncentot);
}

//---------------------------------------------------------------------------//
// functions to set the global-mesh data
//---------------------------------------------------------------------------//
// set the material temperatures at the end of a time-step

template<class MT, class PT>
void Global_Tally<MT,PT>::set_T(const vector<double> &tally)
{
    Require (tally.size() == temperature.size());

  // update cell temperatures
    e_elec_tot = 0.0;
    for (int i = 0; i < temperature.size(); i++)
    {
      // calculate net energy deposition from radiation to electrons
	double delta_E = tally[i] - evol_net[i];

      // calculate new electron temperature in cell
	temperature[i] += delta_E / dedt[i];
	e_elec_tot     += temperature[i] * dedt[i];
    }
}

//---------------------------------------------------------------------------//
// set energies (vol, ss) after source_initialization

template<class MT, class PT>
void Global_Tally<MT,PT>::set_energy_begin(const Source_Init<MT> &source)
{
    Require (num_cells() == source.num_cells());

  // get total energy from volume source
    evoltot   = source.get_evoltot();
    eloss_vol = source.get_eloss_vol();
    nvoltot   = source.get_nvoltot();

  // get new values of evol_net
    double sum = 0.0;
    eint_begin = 0.0;
    for(int i = 0; i < num_cells(); i++)
    {
	evol_net[i] = source.get_evol_net(i+1);
	sum += evol_net[i];
	eint_begin += temperature[i] * dedt[i];
    }
    evolext = evoltot - sum;

  // get total energy from surface source
    esstot    = source.get_esstot();
    eloss_ss  = source.get_eloss_ss();
    nsstot    = source.get_nsstot();

  // after updating the smart pointer to the census (it may still
  // be pointing to the pre-combed census), get the number of 
  // census particles that we are running this timestep
    census    = source.get_census(); 
    ncentot   = census->size();
    ncenrun   = ncentot;
    eradtot_b = eradtot_e;
    eloss_cen = source.get_eloss_cen();

  // calculate integrated energy at beginning of timestep
    eint_begin += eradtot_b;
}

//---------------------------------------------------------------------------//
// set the number of census particles per cell at the end of a timestep

template<class MT, class PT>
void Global_Tally<MT,PT>::set_energy_end(const vector<int> &ncen_,
					 const vector<double> &ecen_, 
					 const vector<double> &ewpl_,
					 double ecent, double eesc)
{
    Require (ncen_.size() == ncen.size());
    Require (ecen_.size() == ecen.size());
    Require (ecen.size()  == ncen.size());
    
  // get new numbers and energy of census particles per cell
    int ncentot = 0;
    ecentot = 0.0;
    for (int i = 0; i < ncen.size(); i++)
    {
	ncen[i]  = ncen_[i];
	ncentot += ncen[i];
      	ecen[i]  = ecen_[i];
	ecentot += ecen[i];
    }
    Check (ncentot == census->size());
    Check (fabs(ecentot - ecent) <= 1.0e-6 * ecent);

  // get the energy-weighted path length for cycle-average 
  // radiation temperature and energy
    for (int i = 0; i < ewpl.size(); i++)
	ewpl[i] = ewpl_[i];

  // update the ecentot energy
    eradtot_e = ecent;

  // update the energy loss in this timestep
    e_escape = eesc;

  // calculate the energy balance tables
    eint_end = e_elec_tot + eradtot_e;

  // do the cycle updates

  // energy in
    double e_in   = esstot + evolext;
    e_in_probtot += e_in; 

  // energy escaped
    e_esc_probtot += e_escape;

  // energy loss due to sampling
    double elosstot = eloss_vol + eloss_ss + eloss_cen;
    e_loss_probtot += elosstot; 
}

//---------------------------------------------------------------------------//
// functions to update the global-mesh objects
//---------------------------------------------------------------------------//
// update a mat_state

template<class MT, class PT>
void Global_Tally<MT,PT>::update_Mat(Mat_State<MT> &mat) const
{
    Require (num_cells() == mat.num_cells());

  // update the temperatures in the Mat_State
    for (int cell = 1; cell <= num_cells(); cell++)
	mat.get_T(cell) = temperature[cell-1];
}

//---------------------------------------------------------------------------//
// update the census part of Source_Init

template<class MT, class PT>
void Global_Tally<MT,PT>::update_Source_Init(Source_Init<MT> &source) const 
{
    Require (num_cells() == source.num_cells());

  // update the census part of Source_Init
    for (int cell = 1; cell <= num_cells(); cell++)
    {
	source.set_ncen(cell, ncen[cell-1]);
	source.set_ecen(cell, ecen[cell-1]);
    }
    source.set_ncentot(census->size());
    source.set_ecentot(ecentot);
    source.set_census(census);
}

//---------------------------------------------------------------------------//
// print diagnostics
//---------------------------------------------------------------------------//
// output the Global Tally

template<class MT, class PT>
void Global_Tally<MT,PT>::print(ostream &out) const
{
    out << endl;
    out << ">>> Material State <<<" << endl;
    out << "======================" << endl;
    out << setw(10) << setiosflags(ios::right) << "Cell"
	<< setw(20) << setiosflags(ios::right) << "Temperature"
	<< setw(20) << setiosflags(ios::right) << "Evol-net"
	<< setw(20) << setiosflags(ios::right) << "dE/dT" << endl;
    out << "=================================================="
	<< "====================" << endl;
    out.precision(6);
    out.setf(ios::scientific, ios::floatfield);
    for (int i = 0; i < num_cells(); i++)
    {	
	out << setw(10) << i+1 << setw(20) << temperature[i] << setw(20)
	    << evol_net[i] << setw(20) << dedt[i] << endl;
    }

  // do the energy conservation check

    out.precision(4);
    out.setf(ios::scientific, ios::floatfield);
    out << endl;
    out << ">>> Energy Conservation Check <<<" << endl;
    out << "=================================" << endl;
    out.setf(ios::right, ios::adjustfield);
    out << setw(30) << "Cycle energy check:" << setw(15) << calc_de_cyc()
	<< " (frac:" << setw(15) << calc_frac_cyc() << ")" << endl;
    out << setw(30) << "Accumulated energy check:" << setw(15) <<
	calc_de_tot() << " (frac:" << setw(15) << calc_frac_tot() << ")" 
	<< endl;
    out << setw(30) << "Cycle energy loss:" << setw(15) 
	<< eloss_vol + eloss_ss + eloss_cen << " (frac:" << setw(15) 
	<< (eloss_vol +	eloss_ss + eloss_cen) / eint_begin << ")" << endl;
    out << setw(30) << "Accumulated energy loss:" << setw(15) 
	<< e_loss_probtot << " (frac:" << setw(15) 
 	<< e_loss_probtot / eint_initial << ")" << endl;
    out << setw(30) << "Escaping energy this cycle:" << setw(15) << e_escape 
	<< endl;
    out << setw(30) << "Energy to census :" << setw(15) << eradtot_e << endl;
    out << setw(30) << "Number to census:" << setw(15) << census->size() 
	<< endl;
    out << setw(30) << "Initial Internal Energy:" << setw(15) << eint_begin 
	<< endl;
    out << setw(30) << "Ending Internal Energy:" << setw(15) << eint_end 
	<< endl;
    out << setw(30) << "Initial Material Energy:" << setw(15) 
	<< evoltot - evolext << endl;
    out << setw(30) << "Ending Material Energy:" << setw(15) << e_elec_tot
	<< endl;
    
    out << endl;   
    out << " ** Energy check by Source type per cycle **" << endl;
    out << " -------------------------------------------" << endl;
    out << setw(35) << "Census:" << setw(15) 
	<< "Volume" << setw(15) << "Surface" << endl;
    out << setw(20) << "Number:" << setw(15) << ncenrun << setw(15)
	<< nvoltot  << setw(15) << nsstot << endl;
    out << setw(20) << "Total energy:" << setw(15) << eradtot_b << setw(15)
	<< evoltot  << setw(15) << esstot << endl;
    out << setw(20) << "Energy loss:" << setw(15) << eloss_cen << setw(15)
	<< setw(15) << eloss_vol << setw(15) << eloss_ss << endl;
}

CSPACE

//---------------------------------------------------------------------------//
//                              end of Global_Tally.cc
//---------------------------------------------------------------------------//
