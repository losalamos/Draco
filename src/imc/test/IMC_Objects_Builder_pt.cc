//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   test/IMC_Objects_Builder_pt.cc
 * \author Thomas M. Evans
 * \date   Sat Aug 23 10:16:19 2003
 * \brief  IMC_Objects_Builder instantiation
 * \note   Copyright © 2003 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "mc/OS_Mesh.hh"
#include "../Frequency.hh"
#include "../Multigroup_Particle.hh"
#include "../IMC_Objects_Builder.t.hh"
#include "IMC_Test.hh"

namespace rtt_imc 
{

typedef rtt_mc::OS_Mesh                        MT;
typedef rtt_imc::Multigroup_Frequency          MG;
typedef rtt_imc::Multigroup_Particle<MT>       MGPT;
typedef rtt_imc_test::IMC_CDI_Interface<MGPT>  CIT;
typedef rtt_imc_test::IMC_Flat_Interface<MGPT> FIT;

template class IMC_Objects_Builder<CIT, MT, MG, MGPT>;
template class IMC_Objects_Builder<FIT, MT, MG, MGPT>;

} // end namespace rtt_imc

//---------------------------------------------------------------------------//
//                 end of IMC_Objects_Builder_pt.cc
//---------------------------------------------------------------------------//
