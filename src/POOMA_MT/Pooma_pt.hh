//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   POOMA_MT/Pooma_pt.hh
 * \author Randy M. Roberts
 * \date   Thu Nov 18 14:24:37 1999
 * \brief  
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __POOMA_MT_Pooma_pt_hh__
#define __POOMA_MT_Pooma_pt_hh__

// Configuration includes

#include <POOMA_MT/config.h>

// POOMA headers
#include "Index/Index.h"
#include "Index/NDIndex.h"
#include "AppTypes/Vektor.h"
#include "FieldLayout/CenteredFieldLayout.h"
#include "Meshes/UniformCartesian.h"
#include "Meshes/Cartesian.h"
#include "Meshes/Centering.h"
#include "Field/Field.h"
#include "Field/LField.h"
#include "Field/Assign.h"
#include "Field/GuardCellSizes.h"

// XTM headers
#include "VektorHelper.hh"
#include "ds++/config.hh"
#include "ds++/SP.hh"

#include "traits/MT_traits.hh"

// Standard C++ headers
#include <iterator.h>

const int Dimension_s = 3;

#define VEKTOR(TYPE, SIZE) Vektor<TYPE, SIZE>

#endif                          // __POOMA_MT_Pooma_pt_hh__

//---------------------------------------------------------------------------//
//                              end of POOMA_MT/Pooma_pt.hh
//---------------------------------------------------------------------------//
