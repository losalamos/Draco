//----------------------------------*-C++-*----------------------------------//
// read_input.cc
// Scott Turner
// 17 March 1998
//---------------------------------------------------------------------------//
// @> Read and test input.
//---------------------------------------------------------------------------//

#include "sn/read_input.hh"

#include <iostream.h>
#include <fstream.h>
#include <stdlib.h>

  void read_input::read_data()
  {

    ifstream input_file;
    input_file.open("snpp.in");

    if (input_file.bad())
    {
      cerr << "Error: Could not open snpp.in" << endl;
      exit (8);
    }

    input_file >> it   >> jt     >> kt   >> mm   >> isct;
    input_file >> dx   >> dy     >> dz   >> epsi;
    input_file >> ifxg >> iprint;
    input_file >> ibl  >> ibb    >> ibfr;

    input_file.close();

    // test for input errors and set parameters based on input

    if ( mm != 3 && mm != 6 )
    {
      cerr << "Error: mm must equal 3 or 6" << endl;
      exit (8);
    }

    if ( isct != 0 && isct != 1 )
    {
      cerr << "Error: isct must equal 0 or 1" << endl;
      exit (8);
    }

  }

  void read_input::get_basic_data( int &it_l, int &jt_l, int &kt_l, int &mm_l,
                                   int &isct_l )

  {
    it_l = it;
    jt_l = jt;
    kt_l = kt;
    mm_l = mm;
    isct_l = isct;
  }

  void read_input::get_all_data( int &it_l, int &jt_l, int &kt_l, int &mm_l,
                                 int &isct_l, int &ibl_l, int &ibb_l,
                                 int &ibfr_l, int &iprint_l, int &ifxg_l,
                                 REAL &dx_l, REAL &dy_l, REAL &dz_l,
                                 REAL &epsi_l )

  {
    it_l = it;
    jt_l = jt;
    kt_l = kt;
    mm_l = mm;
    isct_l = isct;
    ibl_l = ibl;
    ibb_l = ibb;
    ibfr_l = ibfr;
    iprint_l = iprint;
    ifxg_l = ifxg;
    dx_l = dx;
    dy_l = dy;
    dz_l = dz;
    epsi_l = epsi;
  }

//---------------------------------------------------------------------------//
//                              end of read_input.cc
//---------------------------------------------------------------------------//

