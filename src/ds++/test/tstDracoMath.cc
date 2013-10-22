//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   ds++/test/tstDracoMath.cc
 * \author Kent G. Budge
 * \date   Wed Nov 10 09:35:09 2010
 * \brief  Test functions defined in ds++/draco_math.hh.
 * \note   Copyright (C) 2010-2013 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "../ScalarUnitTest.hh"
#include "../Release.hh"
#include "../DracoMath.hh"
#include "../Soft_Equivalence.hh"

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void tstabs( rtt_dsxx::UnitTest &ut )
{
    using rtt_dsxx::abs;
    
    if (abs(-2.2)==abs(2.2))
        ut.passes("Correctly calculated abs(double)");
    else
        ut.failure("Did NOT correctly calculate abs(double)");
    if (abs(-2.2f)==abs(2.2f))
        ut.passes("Correctly calculated abs(float)");
    else
        ut.failure("Did NOT correctly calculate abs(float)");
    if (abs(-2)==abs(2))
        ut.passes("Correctly calculated abs(int)");
    else
        ut.failure("Did NOT correctly calculate abs(int)");
    if (abs(-2L)==abs(2L))
        ut.passes("Correctly calculated abs(long)");
    else
        ut.failure("Did NOT correctly calculate abs(long)");
    return;
}

//---------------------------------------------------------------------------//
void tstconj(rtt_dsxx::UnitTest &ut)
{
    if (rtt_dsxx::soft_equiv(rtt_dsxx::conj(3.5), 3.5))
	ut.passes("conj(double) is correct");
    else
	ut.failure("conj(double) is NOT correct");

    std::complex<double> c(2.7, -1.4);
    if (rtt_dsxx::soft_equiv((rtt_dsxx::conj(c)*c).real(), rtt_dsxx::square(abs(c))))
	ut.passes("conj(std::complex) is correct");
    else
	ut.failure("conj(std::complex) is NOT correct");
    return;
}

//---------------------------------------------------------------------------//
void tstcube(rtt_dsxx::UnitTest &ut)
{
    if (rtt_dsxx::soft_equiv(rtt_dsxx::cube(2.0), 8.0))
        ut.passes("rtt_dsxx::square function returned correct double");
    else
        ut.failure("rtt_dsxx::square function did NOT return correct double.");
    return;
}

//---------------------------------------------------------------------------//
void tstpythag( rtt_dsxx::UnitTest & ut )
{
    using rtt_dsxx::pythag;
    if (rtt_dsxx::soft_equiv(pythag(3.0e307, 4.0e307), 5.0e307))
	ut.passes("pythag correct");
    else
	ut.failure("pythag NOT correct");
    if (rtt_dsxx::soft_equiv(pythag(4.0e307, 3.0e307), 5.0e307))
	ut.passes("pythag correct");
    else
	ut.failure("pythag NOT correct");
    if (rtt_dsxx::soft_equiv(pythag(0.0, 0.0), 0.0))
	ut.passes("pythag correct");
    else
	ut.failure("pythag NOT correct");
    return;
}

//---------------------------------------------------------------------------//
void tstsign( rtt_dsxx::UnitTest & ut )
{
    using rtt_dsxx::sign;
    if (sign(3.2, 5.6)!=3.2)
	ut.failure("sign: FAILED");
    else
	ut.passes("sign: passed");
    if (sign(4.1, -0.3)!=-4.1)
	ut.failure("sign: FAILED");
    else
	ut.passes("sign: passed");
    return;
}

//---------------------------------------------------------------------------//
void tstsquare( rtt_dsxx::UnitTest & ut )
{
    if (rtt_dsxx::soft_equiv(rtt_dsxx::square(3.0), 9.0))
        ut.passes("square function returned correct double");
    else
        ut.failure("square function did NOT return correct double.");
    return;
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    rtt_dsxx::ScalarUnitTest ut(argc, argv, rtt_dsxx::release);
    try
    {
        tstabs(ut);
	tstconj(ut);
        tstcube(ut);
        tstpythag(ut);
	tstsign(ut);
        tstsquare(ut);
    }
    UT_EPILOG(ut);
}   

//---------------------------------------------------------------------------//
// end of tstDracoMath.cc
//---------------------------------------------------------------------------//
