//----------------------------------*-C++-*----------------------------------//
// Source.cc
// Thomas M. Evans
// Thu May 14 08:45:49 1998
//---------------------------------------------------------------------------//
// @> Source class implementation file
//---------------------------------------------------------------------------//

#include "imctest/Source.hh"
#include "ds++/Assert.hh"

IMCSPACE

using std::ios;

//---------------------------------------------------------------------------//
// constructor
//---------------------------------------------------------------------------//
// standard constructor for the Source class

template<class MT>
Source<MT>::Source(typename MT::CCSF_int &vol_rnnum_, 
		   typename MT::CCSF_int &nvol_,
		   typename MT::CCSF_int &ss_rnnum_, 
		   typename MT::CCSF_int &nss_, 
		   string title, int nvoltot_, int nsstot_, int ncentot_)
    : vol_rnnum(vol_rnnum_), nvol(nvol_), ss_rnnum(ss_rnnum_), nss(nss_),
      census(title.c_str(), ios::in), nvoltot(nvoltot_), nsstot(nsstot_), 
      ncentot(ncentot_), buffer(vol_rnnum.get_Mesh())
{
  // some assertions
    Check (vol_rnnum.get_Mesh() == nvol.get_Mesh());
    Check (nvol.get_Mesh() == ss_rnnum.get_Mesh());
    Check (nvol.get_Mesh() == nss.get_Mesh());
    Check (census);
}
	
CSPACE

//---------------------------------------------------------------------------//
//                              end of Source.cc
//---------------------------------------------------------------------------//
