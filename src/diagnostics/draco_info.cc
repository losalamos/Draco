//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   diagnostics/draco_info.cc
 * \author Kelly Thompson
 * \date   Wednesday, Nov 07, 2012, 18:49 pm
 * \brief  Small executable that prints the version and copyright strings.
 * \note   Copyright (C) 2012-2013 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id: tstScalarUnitTest.cc 6864 2012-11-08 01:34:45Z kellyt $
//---------------------------------------------------------------------------//

#include "diagnostics/config.h"
#include "c4/config.h"
#include "ds++/Release.hh"
#include "ds++/Assert.hh"
#include <iostream>
#include <string>
#include <stdexcept>

//---------------------------------------------------------------------------//
// Draco-6_5_20121113, build date 2012/11/13; build type: DEBUG; DBC: 7
//
// Draco Contributers: 
//     Kelly G. Thompson, Kent G. Budge, Tom M. Evans,
//     Rob Lowrie, B. Todd Adams, Mike W. Buksas,
//     James S. Warsa, John McGhee, Gabriel M. Rockefeller,
//     Paul J. Henning, Randy M. Roberts, Seth R. Johnson,
//     Allan B. Wollaber, Peter Ahrens, Jeff Furnish,
//     Paul W. Talbot, Jae H. Chang, and Benjamin K. Bergen.
//
// Copyright (C) 1995-2012 LANS, LLC
// Build information:
//     Library type   : shared
//     System type    : Linux
//     CUDA support   : disabled
//     MPI support    : enabled
//       mpirun cmd   : mpirun --mca mpi_paffinity_alone 0 -np 
//     OpenMPI support: enabled
//     Diagnostics    : disabled
//     Diagnostics Timing: disabled
//     C++11 Support  : enabled
//       Feature list : HAS_CXX11_AUTO 
//                      HAS_CXX11_NULLPTR 
//                      HAS_CXX11_LAMBDA 
//                      HAS_CXX11_STATIC_ASSERT 
//                      HAS_CXX11_SHARED_PTR 

int main( int /*argc*/, char *argv[] )
{
    using std::cout;
    using std::endl;
    try
    {
        // Print version and copyright information to the screen:
        cout << "\n"
             << rtt_dsxx::release() << "\n\n"
             << rtt_dsxx::copyright() << endl;

//---------------------------------------------------------------------------//
// Build Information
//---------------------------------------------------------------------------//

        cout << "Build information:"
             << "\n    Library type   : "
#ifdef DRACO_SHARED_LIBS
             << "shared"
#else
             << "static"
#endif
             << "\n    System type    : "
#if DRACO_UNAME == Linux
             << "Linux"
#else
             << "Unknown"
#endif
             << "\n    CUDA support   : "
#ifdef HAVE_CUDA
             << "enabled"
#else
             << "disabled"
#endif
           
             << "\n    MPI support    : "
#ifdef C4_MPI
             << "enabled"
             << "\n      mpirun cmd   : " << C4_MPICMD
#else
             << "disabled (c4 scalar mode)"
#endif
             << "\n    OpenMPI support: "
#if USE_OPENMP == ON
             << "enabled"
#else
             << "disabled (c4 scalar mode)"
#endif
             << "\n    Diagnostics    : "
#ifdef DRACO_DIAGNOSTICS
             << DRACO_DIAGNOSTICS
#else
             << "disabled"
#endif
             << "\n    Diagnostics Timing: " 
#ifdef DRACO_TIMING
             << "enabled"
#else
             << "disabled"
#endif
            
//---------------------------------------------------------------------------//
// C++11 Features
//---------------------------------------------------------------------------//
            
             << "\n    C++11 Support  : "
#ifdef DRACO_ENABLE_CXX11
             << "enabled"
             << "\n      Feature list : "
#ifdef HAS_CXX11_CSTDINT_H
             << "HAS_CXX11_CSTDINT_H "
#endif
#ifdef HAS_CXX11_AUTO
             << "\n                     HAS_CXX11_AUTO "
#endif
#ifdef HAS_CXX11_NULLPTR
             << "\n                     HAS_CXX11_NULLPTR "
#endif
#ifdef HAS_CXX11_LAMBDA
             << "\n                     HAS_CXX11_LAMBDA_TEMPLATES "
#endif
#ifdef HAS_CXX11_STATIC_ASSERT
             << "\n                     HAS_CXX11_STATIC_ASSERT "
#endif
#ifdef HAS_CXX11_RVALUE_REFERENCES
             << "\n                     HAS_CXX11_RVALUE_REFERENCES "
#endif
#ifdef HAS_CXX11_DECLTYPE
             << "\n                     HAS_CXX11_DECLTYPE "
#endif
#ifdef HAS_CXX11_CSTDINT_H
             << "\n                     HAS_CXX11_CSTDINT_H "
#endif
#ifdef HAS_CXX11_LONG_LONG
             << "\n                     HAS_CXX11_LONG_LONG "
#endif
#ifdef HAS_CXX11_VARIADIC_TEMPLATES
             << "\n                     HAS_CXX11_VARIADIC_TEMPLATES "
#endif
#ifdef HAS_CXX11_CONSTEXPR
             << "\n                     HAS_CXX11_CONSTEXPR "
#endif
#ifdef HAS_CXX11_SIZEOF_MEMBER
             << "\n                     HAS_CXX11_SIZEOF_MEMBER "
#endif
#ifdef HAS_CXX11_SHARED_PTR
             << "\n                     HAS_CXX11_SHARED_PTR "
#endif
#ifdef HAS_CXX11_ARRAY
             << "\n                     HAS_CXX11_ARRAY "
#endif
#else
             << "disabled"
#endif
             << "\n" << endl;
    }
    catch( rtt_dsxx::assertion &err )
    {
        std::string msg = err.what();
        std::cout << "ERROR: While running " << argv[0] << ", "
             << err.what() << std::endl;;
        return 1;
    }
    catch( std::exception &err )
    {
        std::cout << "ERROR: While running " << argv[0] << ", "
             << err.what() << std::endl;;
        return 1;
    }
    catch( ... )
    {
        std::cout << "ERROR: While running " << argv[0] << ", " 
             << "An unknown C++ exception was thrown" << std::endl;;
        return 1;
    }

    return 0;
}   

//---------------------------------------------------------------------------//
// end of draco_info.cc
//---------------------------------------------------------------------------//
