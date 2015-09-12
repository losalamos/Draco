//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/test/tstXGetopt.cc
 * \author Katherine Wang
 * \date   Wed Nov 10 09:35:09 2010
 * \brief  Test functions defined in ds++/XGetopt.cc
 * \note   Copyright (C) 2015 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <iostream>
#include <map>
#include <string.h>

#include "ds++/Release.hh"
#include "ds++/XGetopt.hh"

//---------------------------------------------------------------------------//
int main(int argc, char *argv[])
{
    int aflag = 0;
    int bflag = 0;
    char *cvalue = NULL;
    int index;

    printf ("aflag = %d, bflag = %d, cvalue = %s\n",
            aflag, bflag, cvalue);

    // Call with optional long_options

    rtt_dsxx::optind=1; // resets global counter (see XGetopt.cc)

    std::map< std::string, char> long_options;
    long_options["add"]    = 'a';
    long_options["append"] = 'b';
    long_options["create"] = 'c';

    int c(0);
    while ((c = rtt_dsxx::xgetopt (argc, argv, "abc:", long_options)) != -1)
    {
        switch (c)
        {
            case 'a':
                aflag = 1;
                break;

            case 'b':
                bflag = 1;
                break;

            case 'c':
                cvalue = rtt_dsxx::optarg;
                break;

            default:
                return 0; // nothin to do.
        }
    }

    printf ("aflag = %d, bflag = %d, cvalue = %s\n",
            aflag, bflag, cvalue);

    for (index = rtt_dsxx::optind; index < argc; index++)
        printf ("Non-option argument %s\n", argv[index]);

    return 0;
}

//---------------------------------------------------------------------------//
// end of tstXGetopt.cc
//---------------------------------------------------------------------------//
