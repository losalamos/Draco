//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   test/tstSampler.cc
 * \author Todd J. Urbatsch
 * \date   Wed Apr  5 09:13:05 2000
 * \brief  Monte Carlo Sampler tests.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "MC_Test.hh"
#include "../Sampler.hh"
#include "../Release.hh"
#include "rng/Rnd_Control.hh"
#include "c4/global.hh"

using namespace std;
using rtt_rng::Rnd_Control;
using rtt_rng::Sprng;
using rtt_dsxx::SP;
#include <vector>


bool passed = true;
#define ITFAILS passed = rtt_mc_test::fail(__LINE__);

//---------------------------------------------------------------------------//
// test the function sample_general_linear_test
//---------------------------------------------------------------------------//
/*!  
 * \brief Tests the sample_general_linear function on a uniform pdf, f(x)
 * = 1, and a linear pdf, f(x) = x - 0.5, for x in [1,2].
 */
void sample_general_linear_test()
{
    // set parameters
    int num_samples        = 10000;
    int num_bins           = 100;
    double avg_samples_per_bin = static_cast<double>(num_samples) /
	static_cast<double>(num_bins);

    // a rough estimate of 4 sigma relative error
    double rough_eps = 4.0 / sqrt(avg_samples_per_bin);

    // make sampling bins 
    vector<double> bin(num_bins, 0.0);
    int bindex = 0;
 
    // make a random number generator controller, rng object
    int seed = 1234567;
    Rnd_Control rand_control(seed);
    Sprng ran_object = rand_control.get_rn(10);

    //-----------------------------------------------------------------------//    
    // uniform pdf
    //-----------------------------------------------------------------------//    
    double x_lo = 0.0;
    double x_hi = 1.0;
    double y_lo = 1.0;
    double y_hi = 1.0;

    double x;
    for (int i = 0; i <= num_samples; i++)
    {
	x = rtt_mc::sampler::sample_general_linear(ran_object, x_lo, x_hi, 
						   y_lo, y_hi);
	// calculate bin index
	bindex = static_cast<int>(x*num_bins);
	Check (bindex >= 0 && bindex < num_bins);
	bin[bindex] += 1.0;
    }

    // normalize the bins
    for (int b = 0; b < num_bins; b++)
	bin[b] = bin[b] / avg_samples_per_bin;

    // check that bins replicate original pdf (to w/in roughly 4 sigma)
    for (int b = 0; b < num_bins; b++)
	if (abs(bin[b]-1.0) > rough_eps) ITFAILS;

    //-----------------------------------------------------------------------//    
    // linear pdf
    //-----------------------------------------------------------------------//    
    x_lo = 1.0;
    x_hi = 2.0;
    y_lo = 0.5;
    y_hi = 1.5;

    double delta_x = x_hi - x_lo;

    // reset bins
    for (int b = 0; b < num_bins; b++)
	bin[b] = 0.0;


    // call sample_general_linear num_samples times
    for (int i = 0; i <= num_samples; i++)
    {
	x = rtt_mc::sampler::sample_general_linear(ran_object, x_lo, x_hi, 
						   y_lo, y_hi);
	// calculate bin index
	bindex = static_cast<int>((x-1.0)*num_bins);
	Check (bindex >= 0 && bindex < num_bins);
	bin[bindex] += 1.0;
    }

    // normalize the bins
    for (int b = 0; b < num_bins; b++)
	bin[b] = bin[b] / avg_samples_per_bin;

    // check that bins replicate original pdf (to w/in roughly 4 sigma)
    for (int b = 0; b < num_bins; b++)
    {
	double func_value = y_lo + (b+0.5)*delta_x/num_bins;
	if (abs(bin[b]-func_value) > rough_eps*func_value) ITFAILS;
    }
}

//---------------------------------------------------------------------------//
// main function for tstSampler
//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    C4::Init(argc, argv);

    // this is a serial test
    if (C4::node())
    {
        C4::Finalize();
        return 0;
    }

    // version tag
    for (int arg = 1; arg < argc; arg++)
        if (string(argv[arg]) == "--version")
        {
            cout << argv[0] << ": version " << rtt_mc::release() << endl;
            C4::Finalize();
            return 0;
        }

    try
    {
        // test the sampling of a general linear pdf
	sample_general_linear_test();
    }
    catch (rtt_dsxx::assertion &ass)
    {
        cout << "Assertion during sample_general_linear_test: " 
	     << ass.what() << endl;
        C4::Finalize();
        return 1;
    }

    // status of test
    cout << endl;
    cout <<     "***********************************" << endl;
    if (passed)
    {
        cout << "**** Sampler Self Test: PASSED ****" << endl;
    }
    cout <<     "***********************************" << endl;

    cout << "Done testing Sampler (always serial)" << endl;

    C4::Finalize();
}

//---------------------------------------------------------------------------//
//                              end of tstSampler.cc
//---------------------------------------------------------------------------//
