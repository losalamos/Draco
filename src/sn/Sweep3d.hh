//----------------------------------*-C++-*----------------------------------//
// Sweep3d.hh
// Scott Turner
// 17 April 1998
//---------------------------------------------------------------------------//
// @> Performs structured sweep.
//---------------------------------------------------------------------------//

#ifndef __sn_Sweep3d_hh__
#define __sn_Sweep3d_hh__

#include "sn/precision.hh"
#include "sn/Cross_section.hh"
#include "sn/Input_edit.hh"
#include "sn/Pre_calcs.hh"
#include "sn/Sn_constants.hh"

#include "ds++/Mat.hh"
using dsxx::Mat1;
using dsxx::Mat2;
using dsxx::Mat3;
using dsxx::Mat4;

class Sweep3d
{

    public:

        // use the default constructor
        // Sweep3d();

        // use the default destructor
        // ~Sweep3d();

        void do_sweep( Input_edit &data,   Cross_section &xsec,
                       Sn_constants &sn,   Pre_calcs &pre,
                       Mat1<REAL> &lkgs_l, Mat4<REAL> &src_mom,
                       Mat4<REAL> &flux                         );

        void octant_ordering( Input_edit &data );

        void build_angular_source( Input_edit &data,    Sn_constants &sn,
                                   Mat4<REAL> &src_mom, Mat3<REAL> &phi   );

        void edge_and_boundary_set( Input_edit &data,  Mat2<REAL> &phii,
                                    Mat1<REAL> &phij,  Mat1<REAL> &phik,
                                    Mat2<REAL> &bsavv, Mat3<REAL> &bsavz );

        void balance_eqn_no_fixup( Input_edit &data, Sn_constants &sn,
                                   Pre_calcs &pre,   Mat3<REAL> &phi,
                                   Mat2<REAL> &phii, Mat1<REAL> &phij,
                                   Mat1<REAL> &phik                    );

        void balance_eqn_with_fixup( Input_edit &data, Sn_constants &sn,
                                     Pre_calcs &pre,   Cross_section &xsec,
                                     Mat3<REAL> &phi,  Mat2<REAL> &phii,
                                     Mat1<REAL> &phij, Mat1<REAL> &phik     );

        void save_boundary( Input_edit &data,  Mat2<REAL> &bsavv,
                            Mat3<REAL> &bsavz, Mat1<REAL> &phij,
                            Mat1<REAL> &phik                      );
                           

        void cell_boundary_leakage( Input_edit &data, Sn_constants &sn,
                                    Mat3<REAL> &fh_i, Mat3<REAL> &fv_i,
                                    Mat3<REAL> &fz_i, Mat2<REAL> &phii,
                                    Mat1<REAL> &phij, Mat1<REAL> &phik  );

        void sweep_reversal( Input_edit &data );

        void flux_moments( Input_edit &data, Sn_constants &sn,
                           Mat4<REAL> &flux, Mat3<REAL> &phi   );

        void problem_boundary_leakage( Input_edit &data, Mat1<REAL> &lkgs_l,
                                       Mat3<REAL> &fh_i, Mat3<REAL> &fv_i,
                                       Mat3<REAL> &fz_i                      );

  private:

      int iq;      // loop variable for the number of quadrants, may be
                   //   calculated, rather than incremented in some cases
      int iqp;     // quadrant number(iq) + 1
      int i_cm;    // used in octant sweep order setup
      int mz;      // used in octant sweep order setup
      int ilim;    // used in octant sweep order setup
      int ishift;  // used in octant sweep order setup
      int jshift;  // used in octant sweep order setup
      int kshift;  // used in octant sweep order setup
      int jlow;    // used in octant sweep order setup
      int jhigh;   // used in octant sweep order setup
      int klow;    // used in octant sweep order setup
      int khigh;   // used in octant sweep order setup
      int ioff;    // used in octant sweep order setup
      int iop;     // loop index counter for sweep over pair of quadrants
      int k1;      // index for looping over tsi
      int ih;      // index for looping over eta
      int iz;      // first index of phi during balance equation
      int isw;     // calculated loop variable for cells in the x-direction

};

#endif                          // __sn_Sweep3d_hh__

//---------------------------------------------------------------------------//
//                              end of Sweep3d.hh
//---------------------------------------------------------------------------//

