//----------------------------------*-C++-*----------------------------------//
// Input_edit.hh
// Scott Turner
// 17 April 1998
//---------------------------------------------------------------------------//
// @> Read, test, and manipulate input.
//---------------------------------------------------------------------------//

#ifndef __sn_Input_edit_hh__
#define __sn_Input_edit_hh__

#include "sn/precision.hh"

#include "ds++/Mat.hh"
using dsxx::Mat1;

class Input_edit
{

    public:

        // use the default contructor
        // Input_edit();

        // use the default destructor
        // ~Input_edit();

        void read_edit_data();

        int it()     const { return it_p;     }
        int jt()     const { return jt_p;     }
        int kt()     const { return kt_p;     }
        int mm()     const { return mm_p;     }
        int isct()   const { return isct_p;   }
        int isctp()  const { return isctp_p;  }
        int nm()     const { return nm_p;     }
        int ibl()    const { return ibl_p;    }
        int ibb()    const { return ibb_p;    }
        int ibfr()   const { return ibfr_p;   }
        int iprint() const { return iprint_p; }
        int ifxg()   const { return ifxg_p;   }
        int jbdim()  const { return jbdim_p;  }
        int kbdim()  const { return kbdim_p;  }
        int itmm()   const { return itmm_p;   }
        int maxop()  const { return maxop_p;  }

        REAL dx()        const { return dx_p;    }
        REAL dy()        const { return dy_p;    }
        REAL dz()        const { return dz_p;    }
        REAL epsi()      const { return epsi_p;  }
        REAL hi( int i ) const { return hi_p(i); }
        REAL hj( int j ) const { return hj_p(j); }
        REAL hk( int k ) const { return hk_p(k); }

    private:

        int it_p;      // total number of mesh cells in the x-direction
        int jt_p;      // total number of mesh cells in the y-direction
        int kt_p;      // total number of mesh cells in the z-direction
        int mm_p;      // number of   quadrature points (angles) per quadrant
        int isct_p;    // legendre order of scattering
        int isctp_p;   // isct_p + 1
        int nm_p;      // isctp_p^2, number of flux and source angular moments
        int ibl_p;     // left boundary condition 0/1 = vacuum/reflective
        int ibb_p;     // bottom boundary condition 0/1 = vacuum/reflective
        int ibfr_p;    // front boundary condition 0/1 = vacuum/reflective
        int iprint_p;  // 0/1 = no/yes print fluxes
        int ifxg_p;    // 0/1 = no/yes set to zero flux fixup
        int jbdim_p;   // dimension parameter for boundary arrays
        int kbdim_p;   // dimension parameter for boundary arrays
        int itmm_p;    // number of x-direction cells * number of angles
        int maxop_p;   // loop index for quadrant sweep

        REAL dx_p;     // mesh size in the x-direction
        REAL dy_p;     // mesh size in the y-direction
        REAL dz_p;     // mesh size in the z-direction
        REAL epsi_p;   // convergence precision or,
                       //   if negative, then the number of iterations to do

        Mat1<REAL> hi_p;  // commonly used differencing factor
        Mat1<REAL> hj_p;  // commonly used differencing factor
        Mat1<REAL> hk_p;  // commonly used differencing factor

};

#endif                          // __sn_Input_edit_hh__

//---------------------------------------------------------------------------//
//                              end of Input_edit.hh
//---------------------------------------------------------------------------//

