//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   diagnostics/draco_info_main.cc
 * \author Kelly Thompson
 * \date   Wednesday, Nov 07, 2012, 18:49 pm
 * \brief  Small executable that prints the version and copyright strings.
 * \note   Copyright (C) 2012-2015 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id: tstScalarUnitTest.cc 6864 2012-11-08 01:34:45Z kellyt $
//---------------------------------------------------------------------------//

#include "draco_info.hh"
#include "ds++/Assert.hh"
#include <iostream>
#include "ds++/XGetopt.hh"

int main( int argc, char *argv[] )
{
    using std::cout;
    using std::endl;
    try
    {
        bool version(false);
        bool brief(false);
        rtt_diagnostics::DracoInfo di;

        // Preparing to parse command line arguments.
	int c;
	// rtt_dsxx::optind=1; // resets global counter (see XGetopt.cc)
        std::map< std::string, char> long_options;
        long_options["version"] = 'v';
        long_options["brief"]   = 'b';

        // for( int iargc=1; iargc<argc; ++iargc )
        // {
            while ((c = rtt_dsxx::xgetopt (argc, argv, (char*)"vb", long_options)) != -1)
            {
                switch (c)
                {
                    case 'v': // --version
                        version = true;
                        break;

                    case 'b': // --brief
                        brief = true;
                        break;

                    default:
                        break;
                }
            }
            //}

        if( version )
            cout << di.versionReport();
        else if( brief )
            cout << di.briefReport();
        else
            cout << di.fullReport();
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
// end of draco_info_main.cc
//---------------------------------------------------------------------------//
