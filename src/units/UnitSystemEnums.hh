//----------------------------------*-C++-*----------------------------------//
/*! \file   UnitSystemEnums.hh
 *  \author Kelly Thompson
 *  \brief  This file contains enums, conversion factors and labels that help
 *          define a UnitSystem. 
 *  \date   Fri Oct 24 15:57:09 2003
 *  \note   Copyright (C) 2016-2019 Triad National Security, LLC.
 *          All rights reserved.
 */
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//

#ifndef __units_UnitSystemEnums_hh__
#define __units_UnitSystemEnums_hh__

#include "ds++/config.h"
#include <string>

namespace rtt_units {

//========================================//
// ENUMERATED LENGTH TYPES
//========================================//

enum Ltype {
  L_null = 0, //!< no length type
  L_m = 1,    //!< meters (m)
  L_cm = 2    //!< centimeters (1 m = 100 cm)
};

int const num_Ltype(3);
double const L_cf[num_Ltype] = {0.0, 1.0, 100.0};
std::string const L_labels("NA,m,cm");
std::string const L_long_labels("no length unit specified,meter,centimeter");

//========================================//
// ENUMERATED MASS TYPES
//========================================//

enum Mtype {
  M_null = 0, //!< no mass type
  M_kg = 1,   //!< kilogram (kg)
  M_g = 2     //!< gram (1 g = 1e-3 kg)
};

int const num_Mtype(3);
double const M_cf[num_Mtype] = {0.0, 1.0, 1000.0};
std::string const M_labels("NA,kg,g");
std::string const M_long_labels("no mass unit specified,kilogram,gram");

//========================================//
// ENUMERATED TIME TYPES
//========================================//

enum ttype {
  t_null, //!< no time type
  t_s,    //!< seconds (s)
  t_ms,   //!< milliseconds (1 ms = 1e-3 s)
  t_us,   //!< microseconds (1 us = 1e-6 s)
  t_sh,   //!< shakes       (1 ns = 1e-8 s)
  t_ns    //!< nanoseconds  (1 ns = 1e-9 s)
};

int const num_ttype(6);
double const t_cf[num_ttype] = {0.0, 1.0, 1.0e3, 1.0e6, 1.0e8, 1.0e9};
std::string const t_labels("NA,s,ms,us,sh,ns");
std::string const t_long_labels(
    "no time unit specified,second,milisecond,microsecond,shake,nanosecond");

//========================================//
// ENUMERATED TEMPERATURE TYPES
//========================================//

enum Ttype {
  T_null, //!< no temperature type
  T_K,    //!< Kelvin
  T_keV, //!< keV     (1 K = 8.61733238496e-8 keV or 1 keV = 1.1604519302808940e7 K)
  // This conversion factor between K and keV must agree with the value
  // given in PhysicalConstants.hh.
  T_eV //!< eV      (1 K = 8.61733238496e-5 keV or 1 eV = 11604.519302808940 K)
};

int const num_Ttype(4);
double const T_cf[num_Ttype] = {0.0, 1.0, 1.0 / 1.1604519302808940e7,
                                1.0 / 1.1604519302808940e4};
std::string const T_labels("NA,K,keV,eV");

//========================================//
// ENUMERATED CURRENT TYPES
//========================================//

enum Itype {
  I_null, //!< no current type0.001
  I_amp   //!< Amp (SI)
};

int const num_Itype(2);
double const I_cf[num_Itype] = {0.0, 1.0};
std::string const I_labels("NA,Amp");

//========================================//
// ENUMERATED ANGLE TYPES
//========================================//

enum Atype {
  A_null, //!< no angle type
  A_rad,  //!< Radian (SI)
  A_deg   //!< degrees (PI rad = 180 deg)
};

int const num_Atype(3);
double const A_cf[num_Atype] = {0.0, 1.0, 57.295779512896171};
std::string const A_labels("NA,rad,deg");

//========================================//
// ENUMERATED QUANTITY TYPES
//========================================//

enum Qtype {
  Q_null, //!< no quantity type
  Q_mol   //!< mole (SI)
};

int const num_Qtype(2);
double const Q_cf[num_Qtype] = {0.0, 1.0};
std::string const Q_labels("NA,mol");

//---------------------------------------------------------------------------//
// HELPER FUNCTIONS
//---------------------------------------------------------------------------//

//! Extract unit labels from list in UnitSystemEnums.hh.
DLL_PUBLIC_units std::string setUnitLabel(size_t const pos,
                                          std::string const &labels);

} // end namespace rtt_units

#endif // __units_UnitSystemEnums_hh__

//---------------------------------------------------------------------------//
// end of UnitSystemEnums.hh
//---------------------------------------------------------------------------//
