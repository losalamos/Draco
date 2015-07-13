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
#include <stdio.h>
#include <string.h>

#include "ds++/Release.hh"
#include "ds++/XGetopt.h"

//---------------------------------------------------------------------------//
int main(int argc, char *argv[])
{
  int aflag = 0;
  int bflag = 0;
  char *cvalue = NULL;
  int index;
  int c;

  opterr = 0;

  while ((c = getopt (argc, argv, "abc:")) != -1)
    switch (c)
      {
      case 'a':
        aflag = 1;
        break;
      case 'b':
        bflag = 1;
        break;
      case 'c':
        cvalue = optarg;
        break;
     
      default:
        return 0; // nothin to do.
      }

  printf ("aflag = %d, bflag = %d, cvalue = %s\n",
          aflag, bflag, cvalue);

  for (index = optind; index < argc; index++)
    printf ("Non-option argument %s\n", argv[index]);
 
  return 0;
}

//---------------------------------------------------------------------------//
// end of tstXGetopt.cc
//---------------------------------------------------------------------------//



