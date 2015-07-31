//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/XGetopt.cc
 * \author Katherine Wang
 * \date   Wed Nov 10 09:35:09 2010
 * \brief  Command line argument handling similar to getopt.
 * \note   Copyright (C) 2015 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <iostream>
#include <stdio.h>
#include <string.h>
#include <sstream>
#include <iterator>
#include <vector>
#include <map>
#include "XGetopt.hh"

namespace rtt_dsxx
{

//---------------------------------------------------------------------------//
//  getopt
//
//  Return code >0 if argv contains known shortopts.
//---------------------------------------------------------------------------//

char    *optarg;                // global argument pointer
int      optind = 1;            // global argv index (set to 1, to skip exe name)

int getopt( int argc, char **& argv, std::string const & shortopts,
            longopt_map map )
{
    // convert argv into a vector<string>...
    std::vector< std::string > options( argv, argv + argc );

    // convert shortopts into vector<char> + vector<int>...
    std::vector< char> vshortopts;
    std::vector< int > vshortopts_hasarg; // map<char,bool>
    for( size_t i = 0; i < shortopts.size(); ++i )
    {
        vshortopts.push_back( shortopts[ i ] );
        vshortopts_hasarg.push_back( 0 ); // assume no required argument.
        if( i + 1 < shortopts.size() && shortopts[ i + 1 ] == std::string( ":" )[ 0 ] )
        {
            vshortopts_hasarg[ i ] = 1;
            ++i;
        }
    }

    // Look for command line arguments that match provided list:
    for( ; optind < static_cast<int>( options.size() ); ++optind )
    {
        // stop processing command line arguments
        if( options[ optind ] == std::string( "--" ) )
        {
            return -1;
        }

        // consider single letter options here.
        for( size_t j = 0; j < vshortopts.size(); ++j )
        {
            if( options[ optind ] ==
                std::string( "-" ) + std::string( 1, vshortopts[ j ] ) )
            {

                if( vshortopts_hasarg[ j ] == 1 )
                {
                    ++optind;
                    optarg = argv[ optind ];
                }
                ++optind;
                std::cout << std::endl;
                return vshortopts[ j ];
            }
        }

        // consider string-based optons here.
        if( options[ optind ].substr( 0, 2 ) == std::string( "--" ) )
        {
            for( longopt_map::iterator it = map.begin(); it != map.end(); ++it )
            {
                if( options[ optind ] == std::string( "--" ) + ( it->first ) )
                {
                    // what is the short arg equivalent.
                    size_t j( 0 );
                    for( ; j < vshortopts.size(); ++j )
                        if( std::string( 1, vshortopts[ j ] ) == std::string( 1, it->second ) )
                            break;

                    if( vshortopts_hasarg[ j ] == 1 )
                    {
                        ++optind;
                        optarg = argv[ optind ];
                    }
                    ++optind;
                    std::cout << std::endl;
                    return vshortopts[ j ];
                }

            }

        }
    }

    return -1;
}

} // end namespace rtt_dsxx

//---------------------------------------------------------------------------//
// end of XGetopt.cc
//---------------------------------------------------------------------------//
