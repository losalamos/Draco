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
using rtt_mc::global::soft_equiv;
#include <vector>


bool passed = true;
#define ITFAILS passed = rtt_mc_test::fail(__LINE__);

//---------------------------------------------------------------------------//
// make a fake random number generator that always returns "1"
//---------------------------------------------------------------------------//
class Ran_1
{
  public:
    inline double ran() const {return 1;}
};

//---------------------------------------------------------------------------//
// make a fake random number generator that always returns 1/2
//---------------------------------------------------------------------------//
class Ran_half
{
  public:
    inline double ran() const {return 0.5;}
};

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
// Test the sampling a frequency from a planckian distribution at kT.
//---------------------------------------------------------------------------//

void test_sampling_planckian_frequency()
{
    using rtt_mc::global::pi;

    // make a random number generator
    int seed = 1234567;
    Rnd_Control rand_control(seed);
    Sprng ran_object = rand_control.get_rn(10);

    // make a simple check that kT=0 produces hnu=0
    double kT  = 0.0;
    double hnu = rtt_mc::sampler::sample_planckian_frequency(ran_object, kT); 
    if (!soft_equiv(hnu, 0.0))  ITFAILS;

    // test the planck freq sampler when all random numbers are 1
    {
	Ran_1 ran1;
	kT  = 1.0;
	hnu = rtt_mc::sampler::sample_planckian_frequency(ran1, kT);
	if (!soft_equiv(hnu, 0.0))  ITFAILS;
    }

    // test the planck freq sampler when all random numbers are 1/2
    {
	Ran_half ranh;
	kT  = 1.0;
	hnu = rtt_mc::sampler::sample_planckian_frequency(ranh, kT);
	if (!soft_equiv(hnu, -kT*log(0.0625)))  ITFAILS;
    }

    // do some binning and compare to analytic Planck at the 3-sigma level.
    int num_bins     = 100;
    double bin_width = 0.1;
    int num_samples  = 10000;
    kT               = 1.0;
    vector<double> freq_bins(num_bins, 0.0);

    for (int i = 0; i < num_samples; i++)
    {
	hnu = rtt_mc::sampler::sample_planckian_frequency(ran_object, kT);
	int index = hnu/kT/bin_width;
	if (index < num_bins)
	    freq_bins[index] += 1.0;
    }

    for (int b = 0; b < num_bins; b++)
	freq_bins[b] /= num_samples;
    
    for (int b = 0; b < num_bins; b++)
    {
	double freq     = (b+0.5)*bin_width;
	double x        = freq/kT; 

	double analytic = 15.0/std::pow(pi,4)*x*x*x/(exp(x)-1.0) * bin_width;

	double expected_samples = analytic * num_samples;
	double three_sigma_eps  = 3.0 / sqrt(expected_samples);

	if (!soft_equiv(freq_bins[b], analytic, three_sigma_eps))     ITFAILS;
    }
}

//---------------------------------------------------------------------------//
// Test the boolean check for validity of a cdf
//---------------------------------------------------------------------------//

void test_cdf_validity()
{
    // STL typedef for scalar field of doubles
    typedef std::vector<double> sf_double;

    // make a cdf of zero size (that's bad)
    sf_double my_cdf;

    if (rtt_mc::sampler::is_this_cdf_valid(my_cdf))  ITFAILS;

    // now size it and make some negative values (that's bad)
    my_cdf.resize(3);
    my_cdf[0] =  0.0;
    my_cdf[1] = -3.0;
    my_cdf[2] =  6.0;

    if (rtt_mc::sampler::is_this_cdf_valid(my_cdf))  ITFAILS;

    // give it a negative slope such that it's not monotonically increasing
    // (that's bad)
    my_cdf[1] =  9.0;

    if (rtt_mc::sampler::is_this_cdf_valid(my_cdf))  ITFAILS;

    // now make it a nice cdf (that's good)
    my_cdf[0] =  0.0;
    my_cdf[1] =  0.5;
    my_cdf[2] =  1.0;

    if (!rtt_mc::sampler::is_this_cdf_valid(my_cdf))  ITFAILS;
}

//---------------------------------------------------------------------------//
// Test the sampling of a cdf
//---------------------------------------------------------------------------//

void test_cdf_sampling()
{
    // STL typedef for scalar field of doubles
    typedef std::vector<double> sf_double;

    int sampled_bin = -1;
	
    // make a random number generator that always returns 0.5 
    Ran_half ranhalf_gen;
    
    // make a random number generator that always returns 1.0 
    Ran_1    ranone_gen;
    
    // cdf = 5/8, 5/8, 1.0  =>  pdf = 5/8, 0, 3/8
    {
	// make a normalized cdf and check it once
	sf_double my_cdf(3, 0.0);
	my_cdf[0] = 0.625;
	my_cdf[1] = 0.625;
	my_cdf[2] = 1.000;
	
	if (!rtt_mc::sampler::is_this_cdf_valid(my_cdf))  ITFAILS;
	
	// "sample" the cdf with ran=1/2; it had better be the first
	// bin--second bin has zero probability 
	sampled_bin = 
	    rtt_mc::sampler::sample_bin_from_discrete_cdf(ranhalf_gen, my_cdf);

	if (sampled_bin != 0)  ITFAILS;

	// "sample" the cdf with ran=1; it had better be the last bin
	sampled_bin = 
	    rtt_mc::sampler::sample_bin_from_discrete_cdf(ranone_gen, my_cdf);

	if (sampled_bin != 2)  ITFAILS;

    }

    // cdf = 1/2, 1/2, 1/2  =>  pdf = 1, 0, 0
    {
	// make an unnormalized, flat cdf and check it once
	sf_double my_cdf(3, 0.0);
	my_cdf[0] = 0.5;
	my_cdf[1] = 0.5;
	my_cdf[2] = 0.5;
	
	if (!rtt_mc::sampler::is_this_cdf_valid(my_cdf))  ITFAILS;
	
	// "sample" the cdf
	sampled_bin = 
	    rtt_mc::sampler::sample_bin_from_discrete_cdf(ranone_gen, my_cdf);

	// it had better be the first bin since the other bins have no prob
	if (sampled_bin != 0)  ITFAILS;

	// "sample" the cdf
	sampled_bin = 
	    rtt_mc::sampler::sample_bin_from_discrete_cdf(ranhalf_gen, my_cdf);

	// it had better be the first bin since the other bins have no prob
	if (sampled_bin != 0)  ITFAILS;

    }

    // cdf = 1, 2, 3, 4  =>  pdf = 1/4, 1/4, 1/4, 1/4
    {
	// make an unnormalized cdf and check it once
	sf_double my_cdf(4, 0.0);
	my_cdf[0] = 1.0;
	my_cdf[1] = 2.0;
	my_cdf[2] = 3.0;
	my_cdf[3] = 4.0;
	
	if (!rtt_mc::sampler::is_this_cdf_valid(my_cdf))  ITFAILS;
	
	// "sample" the cdf
	sampled_bin = 
	    rtt_mc::sampler::sample_bin_from_discrete_cdf(ranone_gen, my_cdf);

	// it had better be the last bin
	if (sampled_bin != 3)  ITFAILS;
	
	// "sample" the cdf
	sampled_bin = 
	    rtt_mc::sampler::sample_bin_from_discrete_cdf(ranhalf_gen, my_cdf);

	// it had better be the middle bin
	if (sampled_bin != 1)  ITFAILS;
    }	

    // cdf = 1, 2, 3, ..., 1000  =>  pdf = 1/1000, ... 
    {
	// make an unnormalized cdf and check it once
	sf_double my_cdf(1000);
	for (int i = 0; i < 1000; i++)
	    my_cdf[i] = i + 1;

	if (!rtt_mc::sampler::is_this_cdf_valid(my_cdf))  ITFAILS;
	
	// "sample" the cdf
	sampled_bin = 
	    rtt_mc::sampler::sample_bin_from_discrete_cdf(ranone_gen, my_cdf);

	// it had better be the last bin
	if (sampled_bin != 999)  ITFAILS;
	
	// "sample" the cdf
	sampled_bin = 
	    rtt_mc::sampler::sample_bin_from_discrete_cdf(ranhalf_gen, my_cdf);

	// it had better be the bin just below the midpoint
	if (sampled_bin != 499)  ITFAILS;
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

	// test the sampling of a Planckian frequency
	test_sampling_planckian_frequency();

	// test validity of a cumulative distribution function (cdf)
	test_cdf_validity();

	// test sampling a discrete cdf
	test_cdf_sampling();
    }
    catch (rtt_dsxx::assertion &ass)
    {
        cout << "Assertion during tstSampler: " 
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
