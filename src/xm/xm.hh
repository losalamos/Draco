//----------------------------------*-C++-*----------------------------------//
// Copyright 1996 The Regents of the University of California. 
// All rights reserved.
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
// Created on: Mon Sep 16 15:21:48 1996
// Created by: Geoffrey Furnish
// Also maintained by:
//
//---------------------------------------------------------------------------//

#ifndef __xm_xm_hh__
#define __xm_xm_hh__

#define XM_NAMESPACE_BEG namespace xm {
#define XM_NAMESPACE_END }

// The purpose of this file is to implement an expression template math
// facility which is generally useful to Kull.

#include "applicative.hh"
#include "Xpr.hh"
#include "XprBin.hh"
#include "XprUnary.hh"

#undef XM_NAMESPACE_BEG
#undef XM_NAMESPACE_END

// This will not be necessary when the compiler adopt Koenig lookup
using namespace xm;

#endif // __xm_xm_hh__

//---------------------------------------------------------------------------//
// end of xm.hh
//---------------------------------------------------------------------------//
