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
#include "sn/Array.hh"
#include "sn/Cross_section.hh"
#include "sn/Input_edit.hh"
#include "sn/Pre_calcs.hh"
#include "sn/Sn_constants.hh"

class Sweep3d
{

    public:

        // use the default constructor
        // Sweep3d();

        // use the default destructor
        // ~Sweep3d();

        void do_sweep( Input_edit &data, Cross_section &xsec, Sn_constants &sn,
                       Pre_calcs &pre,   REAL *lkgs_l,        Array4D &src_mom,
                       Array4D &flux                                          );

        void edge_and_boundary_init( Input_edit &data, Array3D &fh_i,
                                     Array3D &fv_i,    Array3D &fz_i,
                                     Array2D &bsavv,   Array3D &bsavz );

        void octant_ordering( Input_edit &data );

        void build_angular_source( Input_edit &data, Sn_constants &sn,
                                   Array4D &src_mom, Array3D &phi      );

        void edge_and_boundary_set( Input_edit &data, Array2D &phii,
                                    REAL *phij,       REAL *phik,
                                    Array2D &bsavv,   Array3D &bsavz );

        void balance_eqn_no_fixup( Input_edit &data, Sn_constants &sn,
                                   Pre_calcs &pre,   Array3D &phi,
                                   Array2D &phii,    REAL *phij,
                                   REAL *phik                          );

        void balance_eqn_with_fixup( Input_edit &data, Sn_constants &sn,
                                     Pre_calcs &pre,   Cross_section &xsec,
                                     Array3D &phi,     Array2D &phii,
                                     REAL *phij,       REAL *phik           );

        void save_boundary( Input_edit &data, Array2D &bsavv, Array3D &bsavz,
                            REAL *phij,       REAL *phik                      );

        void cell_boundary_leakage( Input_edit &data, Sn_constants &sn,
                                    Array3D &fh_i,    Array3D &fv_i,
                                    Array3D &fz_i,    Array2D &phii,
                                    REAL *phij,       REAL *phik        );

        void sweep_reversal( Input_edit &data );

        void flux_moments( Input_edit &data, Sn_constants &sn, Array4D &flux,
                           Array3D &phi                                       );
                              
        void problem_boundary_leakage( Input_edit &data, REAL *lkgs_l,
                                       Array3D &fh_i,    Array3D &fv_i,
                                       Array3D &fz_i                    );

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

