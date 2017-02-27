//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ode/rkqs_pt.cc
 * \author Kent Budge
 * \date   Mon Sep 20 15:15:40 2004
 * \brief  Specializations of rkqs
 * \note   Copyright (C) 2016-2017 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

//#include <vector>
//#include "rkqs.hh"
//#include "quad.i.hh"

namespace rtt_ode {

//using std::vector;

//---------------------------------------------------------------------------//
// Field = double, Function = FPderivs
//---------------------------------------------------------------------------//
//
//typedef void (*FPderivs)(double x,
//                         vector<double> const &y,
//                         vector<double> &yderiv);

//template DLL_PUBLIC_ode
//void rkck( vector<double> const &y,
//		   vector<double> const &dydx,
//		   double x,
//		   double h,
//		   std::vector<double> &yout,
//		   std::vector<double> &yerr,
//		   FPderivs);
//
//template DLL_PUBLIC_ode
//void rkqs( vector<double> &y,
//		   vector<double> const &dydx,
//		   double &x,
//		   double htry,
//		   double eps,
//		   vector<double> const &yscal,
//		   double &hdid,
//		   double &hnext,
//		   FPderivs);

//---------------------------------------------------------------------------//
// Field = double, Function = Quad_To_ODE<double (*)(double)>
//---------------------------------------------------------------------------//

// template void rkck(vector<double> const &y,
// 		   vector<double> const &dydx,
// 		   double x,
// 		   double h,
// 		   std::vector<double> &yout,
// 		   std::vector<double> &yerr,
// 		   );
//
//template DLL_PUBLIC_ode
//void rkqs( vector<double> &y,
//		   vector<double> const &dydx,
//		   double &x,
//		   double htry,
//		   double eps,
//		   vector<double> const &yscal,
//		   double &hdid,
//		   double &hnext,
//           Quad_To_ODE<double (*)(double)>);

} // end namespace rtt_ode

//---------------------------------------------------------------------------//
// end of rkqs_pt.cc
//---------------------------------------------------------------------------//
