//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   meshTest/DoubleContainer.hh
 * \author Shawn Pautz, Randy M. Roberts
 * \date   Mon Aug 23 16:30:39 1999
 * \brief  Header and implementation file for DoubleContainer class.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __meshTest_DoubleContainer_hh__
#define __meshTest_DoubleContainer_hh__

namespace rtt_meshTest
{
 
//===========================================================================//
/*! \class DoubleContainer
 *
 * \brief The following class exists to test the MT field types with a
 *            non-trivial type.
 */
//
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

class DoubleContainer
{
  public:
    double data;

    DoubleContainer() : data(0.) {}

    DoubleContainer(double _data) : data(_data) {}

    DoubleContainer(const DoubleContainer& dc) : data(dc.data) {}

    DoubleContainer& operator=(const DoubleContainer& dc)
    {
        data = dc.data;
        return *this;
    }
    
    bool operator==(const DoubleContainer& dc) const
    {
	return data == dc.data;
    }

    bool operator!=(const DoubleContainer& dc) const
    {
	return data != dc.data;
    }

    bool operator<(const DoubleContainer& dc) const
    {
	return data < dc.data;
    }
    
    bool operator<=(const DoubleContainer& dc) const
    {
	return data <= dc.data;
    }

    bool operator>(const DoubleContainer& dc) const
    {
	return data > dc.data;
    }
    
    bool operator>=(const DoubleContainer& dc) const
    {
	return data >= dc.data;
    }
};

} // end namespace rtt_meshTest

#endif                          // __meshTest_DoubleContainer_hh__

//---------------------------------------------------------------------------//
//                              end of meshTest/DoubleContainer.hh
//---------------------------------------------------------------------------//
