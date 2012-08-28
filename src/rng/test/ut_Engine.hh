//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   rng/test/ut_Engine.hh
 * \author Gabriel Rockefeller
 * \date   Mon Aug 27 19:29:53 2012
 * \brief  Header file for ut_Engine test.
 * \note   Copyright (C) 2012 Los Alamos National Security, LLC
 */
//---------------------------------------------------------------------------//
#if defined __GNUC__
#pragma GCC system_header
#endif

#include <Random123/philox.h>
#include <Random123/aes.h>
#include <Random123/threefry.h>
#include <Random123/ars.h>
#include <Random123/conventional/Engine.hpp>
#include <Random123/ReinterpretCtr.hpp>
