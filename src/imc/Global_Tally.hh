//----------------------------------*-C++-*----------------------------------//
// Global_Tally.hh
// Thomas M. Evans
// Wed Jun 17 10:21:19 1998
//---------------------------------------------------------------------------//
// @> Global_Tally class header file
//---------------------------------------------------------------------------//

#ifndef __imc_Global_Tally_hh__
#define __imc_Global_Tally_hh__

//===========================================================================//
// class Global_Tally - 
//
// Purpose : stores global-mesh data on the host processor, ie. temps,
//           specific-heat, energy balance etc
//
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

#include "imc/Names.hh"
#include "imc/Mat_State.hh"
#include "imc/Source_Init.hh"
#include "imc/Tally.hh"
#include "imc/Particle.hh"
#include "imc/Particle_Buffer.hh"
#include "ds++/SP.hh"
#include <vector>
#include <iostream>

IMCSPACE

using std::vector;
using std::ostream;

template<class MT, class PT = Particle<MT> >
class Global_Tally 
{
private:
  // tally info
    vector<double> temperature;
    vector<double> dedt;
    double e_elec_tot;

  // problem energies

  // internal energies
    double eint_initial;
    double eint_begin;
    double eint_end;

  // problem totals
    double e_in_probtot;
    double e_esc_probtot;
    double e_loss_probtot;

  // volume energies and source
    vector<double> evol_net;
    double evoltot;
    double evolext;
    int nvoltot;

  // census
    SP<typename Particle_Buffer<PT>::Census> census;
    vector<int> ncen;
    int ncentot;

  // radiation energies (census)
    double eradtot_b;
    double eradtot_e;

  // surface source energies and source
    double esstot;
    int nsstot;

  // total loss of energy from sampling the source
    double eloss_vol;
    double eloss_ss;
    double eloss_cen;
    double e_escape;

  // energy conservation checks
    inline double calc_de_cyc() const;
    inline double calc_de_tot() const;
    inline double calc_frac_cyc() const;
    inline double calc_frac_tot() const;
    
public:
  // constructor
    Global_Tally(const MT &, const Mat_State<MT> &,
		 const Source_Init<MT> &); 

  // update the global objects with timestep dependent values
    void update_Mat(Mat_State<MT> &) const;
    void update_Source_Init(Source_Init<MT> &) const;
    
  // set values in the global buffer
    void set_T(const vector<double> &);
    void set_cen(const vector<int> &);
    void set_energy_begin(const Source_Init<MT> &);
    void set_energy_end(const vector<int> &, double, double);

  // accessors
    double get_T(int cell) const { return temperature[cell-1]; }
    inline SP<typename Particle_Buffer<PT>::Census> get_census() const;
    int num_cells() const { return temperature.size(); }

  // print output
    void print(ostream &) const;
};

//---------------------------------------------------------------------------//
// overloaded operators
//---------------------------------------------------------------------------//
// stream insertion

template<class MT>
ostream& operator<<(ostream &out, const Global_Tally<MT> &object)
{
    object.print(out);
    return out;
}

//---------------------------------------------------------------------------//
// inline functions
//---------------------------------------------------------------------------//
// return the census

template<class MT, class PT> inline
SP<typename Particle_Buffer<PT>::Census> Global_Tally<MT,PT>::get_census()
    const 
{
    return census;
}

//---------------------------------------------------------------------------//
// calculate the energy conservation / cycle

template<class MT, class PT>
inline double Global_Tally<MT,PT>::calc_de_cyc() const
{   
    double elosstot = eloss_vol + eloss_ss + eloss_cen;
    return eint_end - (eint_begin + (esstot+evolext) - e_escape - elosstot);
}

//---------------------------------------------------------------------------//
// calculate the energy conservation fraction per cycle

template<class MT, class PT>
inline double Global_Tally<MT,PT>::calc_frac_cyc() const
{   
    return calc_de_cyc() / eint_end;
}

//---------------------------------------------------------------------------//
// calculate the energy conservation over the whole problem

template<class MT, class PT>
inline double Global_Tally<MT,PT>::calc_de_tot() const
{   
    return eint_end - (eint_initial + e_in_probtot - e_esc_probtot 
		       - e_loss_probtot);
}

//---------------------------------------------------------------------------//
// calculate the energy conservation fraction over the whole problem

template<class MT, class PT>
inline double Global_Tally<MT,PT>::calc_frac_tot() const
{   
    return calc_de_tot() / eint_end;
}

CSPACE

#endif                          // __imc_Global_Tally_hh__

//---------------------------------------------------------------------------//
//                              end of imc/Global_Tally.hh
//---------------------------------------------------------------------------//
