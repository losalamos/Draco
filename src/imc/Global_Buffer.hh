//----------------------------------*-C++-*----------------------------------//
// Global_Buffer.hh
// Thomas M. Evans
// Wed Jun 17 10:21:19 1998
//---------------------------------------------------------------------------//
// @> Global_Buffer class header file
//---------------------------------------------------------------------------//

#ifndef __imc_Global_Buffer_hh__
#define __imc_Global_Buffer_hh__

//===========================================================================//
// class Global_Buffer - 
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
class Global_Buffer 
{
private:
  // tally info
    vector<double> temperature;
    vector<double> Cv;
    vector<double> evol_net;

  // census info
    vector<int> ncen;
    vector<double> ecen;
    SP<typename Particle_Buffer<PT>::Census> census;
    
public:
  // constructor
    Global_Buffer(const MT &, const Mat_State<MT> &,
		  const Source_Init<MT> &); 

  // calculate the energy depositions and update the temperature
    void update_T(const vector<double> &);
    void update_cen(const vector<int> &);
    void update_Mat(Mat_State<MT> &) const;
    void update_Source_Init(Source_Init<MT> &) const;
    
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
ostream& operator<<(ostream &out, const Global_Buffer<MT> &object)
{
    object.print(out);
    return out;
}

//---------------------------------------------------------------------------//
// inline functions
//---------------------------------------------------------------------------//
// return the census

template<class MT, class PT> inline
SP<typename Particle_Buffer<PT>::Census> Global_Buffer<MT,PT>::get_census()
    const 
{
    return census;
}

CSPACE

#endif                          // __imc_Global_Buffer_hh__

//---------------------------------------------------------------------------//
//                              end of imc/Global_Buffer.hh
//---------------------------------------------------------------------------//
