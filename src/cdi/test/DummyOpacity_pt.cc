//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi/test/DummyOpacity_pt.cc
 * \author Kelly Thompson
 * \date   Wed Nov 14 9:50 2000
 * \brief  DummyOpacity< EnergyPolicy > class instantiation file.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "DummyOpacity.hh"
#include "DummyOpacity.t.hh"

#include "../Opacity.hh"

#include <vector>

namespace rtt_dummy_opacity
{
    // Instantiate both energy policies.

    template class DummyOpacity<Gray>;
    template class DummyOpacity<Multigroup>;

    // We must also instantiate the getOpacity() accessor for all
    // common STL containers. This is done below.  

    // Note that Mat1<double> and vector<double> have the same
    // iterator definition so we only instatiate one of these
    // containers. 

    // --------------------------- //
    // Const Vector Input Iterator //
    // --------------------------- //

//     namespace
// 	{
// 	    namespace constVectorDoubleNS
// 		{
// 		    typedef std::vector<double>::const_iterator InputIterator;
// 		    typedef std::vector<double>::iterator OutputIterator;
// 		}
//  	}
//     template constVectorDoubleNS::OutputIterator DummyOpacity<Gray>::getOpacity(
// 	constVectorDoubleNS::InputIterator temperatureIterator, 
// 	constVectorDoubleNS::InputIterator temperatureIteratorEnd,
// 	constVectorDoubleNS::InputIterator densityIterator, 
// 	constVectorDoubleNS::InputIterator densityIteratorEnd,
// 	constVectorDoubleNS::OutputIterator opacityIterator ) const;
//     template constVectorDoubleNS::OutputIterator DummyOpacity<Multigroup>::getOpacity(
// 	constVectorDoubleNS::InputIterator temperatureIterator, 
// 	constVectorDoubleNS::InputIterator temperatureIteratorEnd,
// 	constVectorDoubleNS::InputIterator densityIterator, 
// 	constVectorDoubleNS::InputIterator densityIteratorEnd,
// 	constVectorDoubleNS::OutputIterator opacityIterator ) const;
//     template constVectorDoubleNS::OutputIterator DummyOpacity<Gray>::getOpacity(
// 	const double targetTemperature,
// 	constVectorDoubleNS::InputIterator densityIterator, 
// 	constVectorDoubleNS::InputIterator densityIteratorEnd,
// 	constVectorDoubleNS::OutputIterator opacityIterator ) const;
//     template constVectorDoubleNS::OutputIterator DummyOpacity<Multigroup>::getOpacity(
// 	const double targetTemperature,
// 	constVectorDoubleNS::InputIterator densityIterator, 
// 	constVectorDoubleNS::InputIterator densityIteratorEnd,
// 	constVectorDoubleNS::OutputIterator opacityIterator ) const;
//     template constVectorDoubleNS::OutputIterator DummyOpacity<Gray>::getOpacity(
// 	constVectorDoubleNS::InputIterator temperatureIterator, 
// 	constVectorDoubleNS::InputIterator temperatureIteratorEnd,
// 	const double targetDensity,
// 	constVectorDoubleNS::OutputIterator opacityIterator ) const;
//     template constVectorDoubleNS::OutputIterator DummyOpacity<Multigroup>::getOpacity(
// 	constVectorDoubleNS::InputIterator temperatureIterator, 
// 	constVectorDoubleNS::InputIterator temperatureIteratorEnd,
// 	const double targetDensity,
// 	constVectorDoubleNS::OutputIterator opacityIterator ) const;

    // ------------------------------- //
    // Non-Const Vector Input Iterator //
    // ------------------------------- //

//      namespace
//  	{
// 	    namespace vectorDoubleNS
// 		{
// 		    typedef std::vector<double>::iterator InputIterator;
// 		    typedef std::vector<double>::iterator OutputIterator;
// 		}
//  	}
//     template vectorDoubleNS::OutputIterator DummyOpacity<Gray>::getOpacity(
// 	vectorDoubleNS::InputIterator temperatureIterator, 
// 	vectorDoubleNS::InputIterator temperatureIteratorEnd,
// 	vectorDoubleNS::InputIterator densityIterator, 
// 	vectorDoubleNS::InputIterator densityIteratorEnd,
// 	vectorDoubleNS::OutputIterator opacityIterator ) const;
//     template vectorDoubleNS::OutputIterator DummyOpacity<Multigroup>::getOpacity(
// 	vectorDoubleNS::InputIterator temperatureIterator, 
// 	vectorDoubleNS::InputIterator temperatureIteratorEnd,
// 	vectorDoubleNS::InputIterator densityIterator, 
// 	vectorDoubleNS::InputIterator densityIteratorEnd,
// 	vectorDoubleNS::OutputIterator opacityIterator ) const;
//     template vectorDoubleNS::OutputIterator DummyOpacity<Gray>::getOpacity(
// 	const double targetTemperature,
// 	vectorDoubleNS::InputIterator densityIterator, 
// 	vectorDoubleNS::InputIterator densityIteratorEnd,
// 	vectorDoubleNS::OutputIterator opacityIterator ) const;
//     template vectorDoubleNS::OutputIterator DummyOpacity<Multigroup>::getOpacity(
// 	const double targetTemperature,
// 	vectorDoubleNS::InputIterator densityIterator, 
// 	vectorDoubleNS::InputIterator densityIteratorEnd,
// 	vectorDoubleNS::OutputIterator opacityIterator ) const;
//     template vectorDoubleNS::OutputIterator DummyOpacity<Gray>::getOpacity(
// 	vectorDoubleNS::InputIterator temperatureIterator, 
// 	vectorDoubleNS::InputIterator temperatureIteratorEnd,
// 	const double targetDensity,
// 	vectorDoubleNS::OutputIterator opacityIterator ) const;
//     template vectorDoubleNS::OutputIterator DummyOpacity<Multigroup>::getOpacity(
// 	vectorDoubleNS::InputIterator temperatureIterator, 
// 	vectorDoubleNS::InputIterator temperatureIteratorEnd,
// 	const double targetDensity,
// 	vectorDoubleNS::OutputIterator opacityIterator ) const;

    // ----------------------------- //
    // Non-Const Mat1 Input Iterator //
    // ----------------------------- //

    // Mat1<double>::iterator == vector<double>::iterator == double *
    // Declaring a Mat1 iterator is redundant (???).

//     namespace
//  	{
// 	    namespace Mat1DoubleNS
// 		{
// 		    typedef rtt_dsxx::Mat1<double>::iterator InputIterator;
// 		    typedef rtt_dsxx::Mat1<double>::iterator OutputIterator;
// 		}
// 	}
//     template Mat1DoubleNS::OutputIterator DummyOpacity<Gray>::getOpacity(
// 	Mat1DoubleNS::InputIterator temperatureIterator, 
// 	Mat1DoubleNS::InputIterator temperatureIteratorEnd,
// 	Mat1DoubleNS::InputIterator densityIterator, 
// 	Mat1DoubleNS::InputIterator densityIteratorEnd,
// 	Mat1DoubleNS::OutputIterator opacityIterator ) const;
//     template Mat1DoubleNS::OutputIterator DummyOpacity<Multigroup>::getOpacity(
// 	Mat1DoubleNS::InputIterator temperatureIterator, 
// 	Mat1DoubleNS::InputIterator temperatureIteratorEnd,
// 	Mat1DoubleNS::InputIterator densityIterator, 
// 	Mat1DoubleNS::InputIterator densityIteratorEnd,
// 	Mat1DoubleNS::OutputIterator opacityIterator ) const;

} // end namespace rtt_cdi_gandolf

//---------------------------------------------------------------------------//
//                          end of DummyOpacity_pt.cc
//---------------------------------------------------------------------------//
