//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   LAMG/CompressedRowStorage.hh
 * \author Randy M. Roberts
 * \date   Fri Jan 21 10:57:12 2000
 * \brief  
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __LAMG_CompressedRowStorage_hh__
#define __LAMG_CompressedRowStorage_hh__

#include <LAMG/config.h>
#include "ds++/Mat.hh"

#include <algorithm>

namespace rtt_LAMG
{
 
//===========================================================================//
/*!
 * \class CompressedRowStorage
 *
 */
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

class CompressedRowStorage 
{

    // NESTED CLASSES AND TYPEDEFS

  public:
    
    typedef dsxx::Mat1<RTT_F90_INTEGER> IMat1;
    typedef dsxx::Mat1<RTT_F90_DOUBLE> DMat1;

    // DATA

  private:
    
    IMat1 rowPointer_m;
    IMat1 colIndex_m;
    DMat1 val_m;
    
  public:

    // CREATORS

    CompressedRowStorage(const IMat1 &rowPointer_in,
			 const IMat1 &colIndex_in,
			 const DMat1 &val_in)
	: rowPointer_m(rowPointer_in),
	  colIndex_m(colIndex_in),
	  val_m(val_in)
    {
	// empty
    }

    template<class IFT, class DFT>
    CompressedRowStorage(const IFT &rowPointer_in,
			 const IFT &colIndex_in,
			 const DFT &val_in)
	: rowPointer_m(rowPointer_in.size()),
	  colIndex_m(colIndex_in.size()),
	  val_m(val_in.size())
    {
	std::copy(rowPointer_in.begin(), rowPointer_in.end(),
		  rowPointer_m.begin());
	std::copy(colIndex_in.begin(), colIndex_in.end(),
		  colIndex_m.begin());
	std::copy(val_in.begin(), val_in.end(),
		  val_m.begin());
    }
    
    // MANIPULATORS
    
    // ACCESSORS

    const IMat1 &rowPointer() const { return rowPointer_m; }
    const IMat1 &colIndex() const { return colIndex_m; }
    const DMat1 &val() const { return val_m; }
    
  private:
    
    // IMPLEMENTATION
};

} // end namespace rtt_LAMG

#endif                          // __LAMG_CompressedRowStorage_hh__

//---------------------------------------------------------------------------//
//                              end of LAMG/CompressedRowStorage.hh
//---------------------------------------------------------------------------//
