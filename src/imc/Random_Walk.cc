//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   imc/Random_Walk.cc
 * \author Thomas M. Evans
 * \date   Fri Jan 24 15:13:23 2003
 * \brief  Random_Walk_Sampling_Tables definitions.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Random_Walk.hh"
#include "mc/Math.hh"
#include <algorithm>

namespace rtt_imc
{

//===========================================================================//
// RANDOM_WALK_SAMPLING_TABLES
//===========================================================================//

//---------------------------------------------------------------------------//
// CONSTRUCTOR
//---------------------------------------------------------------------------//
/*!
 * \brief Constructor.
 */
Random_Walk_Sampling_Tables::Random_Walk_Sampling_Tables()
{
    // build abcissas
    set_abcissas();
    
    // build probabilities for exiting sphere
    set_prob_exit();

    // build probabilities as a function of (a,R1) of staying in sphere
    set_prob_R_cen();
}

//---------------------------------------------------------------------------//
// PUBLIC INTERFACE
//---------------------------------------------------------------------------//
/*!
 * \brief Return the exiting probability for a particle.
 *
 * Here we calculate the probability that a particle will exit a sphere in a
 * given time, \e t.  The exiting probability is given by:
 * \f[ 
 * P_{T}(t) = 1 - 4\pi\int_{0}^{R_{0}}\psi(r,t)r^{2}dr.
 * \f]
 * Here,
 * \f[
 * \psi(r, t) = \frac{1}{2R_{0}^{2}}\sum_{n}
 * \left(\frac{n}{r}\right)e^{-(n\pi)^2 a}
 * \sin\frac{n\pi r}{R_{o}}.
 * \f]
 * We define \e a as follows
 * \f[
 * a = \frac{Dt}{R_{0}^{2}},
 * \f]
 * where 
 * \f[
 * D = \frac{c}{3(1-f)\sigma_{R}},
 * \f]
 * with \e f begin the Fleck factor.  Taking the integral, we arrive at the
 * following equation for the probability,
 * \f[
 * P_{T}(t) = 1 - \sum_{n}2(-1)^{n-1}e^{-n^{2}\pi^{2}a}.
 * \f]
 * Thus, for a given \e a we can calculate the exit probability by log-log
 * interpolation on precomputed tables of probability versus \e a.
 * 
 * \param t time particle is in flight [shakes]
 * \param D random walk diffusion coefficient [cm^2/shake]
 * \param Ro sphere radius [cm]
 * \return probability of a random walk particle exiting the sphere
 */
double Random_Walk_Sampling_Tables::get_prob_exit(double t, 
						  double D, 
						  double Ro) const
{
    using rtt_mc::global::linear_interpolate;
    using std::lower_bound;

    Require (t >= 0.0);
    Require (Ro >= 0.0);
    Require (D >= 0.0);

    // calculate A
    double A = D * t / (Ro * Ro);
    Check (A >= 0.0);

    // if a > 20 then probability is  1 that the particle will exit sphere
    if (A >= 20.0)
	return 1.0;

    // do a binary search on A
    const double *ptr = lower_bound(a, a + 41, A);
    
    // calculate the index of a (index points to first value >= A)
    int index = ptr - a;
    Check (index >= 0 && index < 41);

    // calculate probability from table; linear interpolation appears to work
    // better than log-log interpolation
    double value = linear_interpolate(a[index-1], a[index], 
				      prob_exit[index-1], prob_exit[index],
				      A);

    // return value
    Ensure(value >= prob_exit[index-1] && value <= prob_exit[index]);
    return value;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Return the elapsed time for a particle that exits the sphere.
 *
 * If a particle exits the surface of a Random Walk sphere, it is necessary
 * to determine how much time was spent in the random walk.  We do this be
 * sampling the exit-probability table for \e a where:
 * \f[
 * a = \frac{Dt}{R_{0}^{2}}.
 * \f]
 * We use the same random number that was used to succesfully determine the
 * probability of exiting the sphere:
 * \f[
 * \xi=P_{T}(t).
 * \f]
 * We invert the table to get \e a; then,
 * \f[
 * t = \frac{a R_{0}^{2}}{D}.
 * \f]
 *
 * \param D random walk diffusion coefficient [cm^2/shake]
 * \param Ro sphere radius [cm]
 * \param ran random number that was used to sample the probability of
 * exiting the sphere
 * \return elapsed time spent in random walk for a particle that exits the
 * sphere 
 */
double Random_Walk_Sampling_Tables::get_elapsed_time(double D,
						     double Ro,
						     double ran) const
{
    using rtt_mc::global::linear_interpolate;
    using std::lower_bound;

    Require (ran >= 0.0 && ran <= 1.0);
    Require (Ro >= 0.0);
    Require (D > 0.0);

    // look up A from the table
    const double *ptr = lower_bound(prob_exit, prob_exit + 41, ran);

    // calculate the index of prob_exit (index points to first value >=
    // P_exit 
    int index = ptr - prob_exit;
    Check (index >= 0 && index < 41);

    // calculate a from table using linear interpolation
    double value = linear_interpolate(prob_exit[index-1], prob_exit[index],
				      a[index-1], a[index], 
				      ran);
    Check (value >= a[index-1] && value <= a[index]);

    // calculate the elapsed time
    value = value * Ro * Ro / D;

    // return value
    return value;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Return the elapsed time for a particle that exits the sphere.
 *
 * When a particle goes to census while still doing random walk we must
 * determine the radius of the sphere that the particle reaches in \f$
 * t_{cen} \f$ (the time to census).  We can do this by sampling,
 * \f[
 * \xi = \frac{4\pi\int_{0}^{R_{1}}\psi(r,t_{cen})r^{2}dr}{P_{R}(t_{cen})}, 
 * \f]
 * where 
 * \f[
 * P_{R}(t_{cen}) = 4\pi\int_{0}^{R_{0}}\psi(r,t_{cen})r^{2}dr.
 * \f]
 * This sampling is done by a two-dimensional table lookup,
 * \f[
 * \xi = \frac{\sum_{n}2e^{-n^{2}\pi^{2}a}\left[
 * \frac{\sin(n\pi b)}{n\pi} - b\cos(n\pi b)\right]}
 * {\sum_{n}2(-1)^{n-1}e^{-n^{2}\pi^{2}a}},
 * \f]
 * where,
 * \f[
 * a = \frac{Dt_{cen}}{R_{0}^{2}},
 * \f]
 * and 
 * \f[
 * b = \frac{R_{1}}{R_{0}}.
 * \f]
 * Thus,
 * \f[
 * \xi = \mbox{Prob}(a,b) = \frac{P_{R}(a,b)}{P_{R}(a)}.
 * \f]
 * To summarize, we do the following:
 * -# calculate \e a
 * -# invert the tables to get Prob(\e a, \e b )
 * -# interpolate to calculate R1
 * -# MORE DETAILS!!!!!
 * .
 *
 * \param t time to census in [shakes]
 * \param D random walk diffusion coefficient [cm^2/shake]
 * \param Ro sphere radius [cm]
 * \param ran random number used to sample tables
 * \return elapsed time spent in random walk for a particle that exits the
 * sphere 
 */
double Random_Walk_Sampling_Tables::get_radius(double t,
					       double D,
					       double Ro,
					       double ran) const
{
    using rtt_mc::global::linear_interpolate;
    using std::lower_bound;

    Require (t >= 0.0);
    Require (Ro >= 0.0);
    Require (D >= 0.0);

    // calculate A
    double A = D * t / (Ro * Ro);
    Check (A >= 0.0);

    // indices
    int a_index = 0;
    int b_index = 0;

    // find (high) index on abcissa table

    // if A > 20 then hard-cap A (so it doesn't go off the table)
    if (A > 20.0)
    {
	A       = 20.0;
	a_index = 33;
    }

    // if A > 0.4 then use the last index
    else if (A > 0.4)
    {
	a_index = 33;
    }

    // else we need to do a binary search to find high a index
    else
    {
	// do a binary search on A
	const double *ptr = lower_bound(a_R, a_R + 34, A);
	
	// calculate the index of a (index points to first value >= A)
	a_index = ptr - a_R;
	Check (a_index >= 0 && a_index < 33);
    }

    // now we do our two-dimensional lookup, this may take several attempts
    // because we don't a priori know the shape of Prob for each b
    bool not_done = true;
    double low_R  = 0.0;
    double high_R = 0.0;
    double B      = 0.0;

    // get an initial b_index to search on
    const double *b_ptr = lower_bound(prob_R_cen[a_index-1],
				      prob_R_cen[a_index-1] + 13,
				      ran);
    
    // calculate b_index
    b_index = b_ptr - prob_R_cen[a_index-1];
    Check (b_index >= 0 && b_index < 13);

    // now invert the tables
    while (not_done)
    {

	// now interpolate with A to get low and high R probs
	low_R  = linear_interpolate(a_R[a_index-1], a_R[a_index], 
				    prob_R_cen[a_index-1][b_index-1],
				    prob_R_cen[a_index][b_index-1],
				    A);
	high_R = linear_interpolate(a_R[a_index-1], a_R[a_index], 
				    prob_R_cen[a_index-1][b_index],
				    prob_R_cen[a_index][b_index],
				    A);
	
	// if low_R <= ran <= high_R then we can interpolate to get b and we
	// are done
	if (ran >= low_R && ran <= high_R)
	{
	    // interpolate to get B
	    B = linear_interpolate(low_R, high_R, 
				   b[b_index-1], b[b_index], ran);
	    Check (B >= 0 && B <= 1.0);

	    // we are done
	    not_done = false;
	}
	
	// if ran < low_R then we need to lower b_index
	else if (ran < low_R)
	{
	    b_index--;
	    Check (b_index >= 0);
	}

	// if ran > high_R then we need to raise b_index
	else if (ran > high_R)
	{
	    b_index++;
	    Check (b_index < 13); 
	}
    }
    
    // now calculate the radius
    double radius = B * Ro;
    return radius;
}


//---------------------------------------------------------------------------//
// PRIVATE IMPLEMENTATION
//---------------------------------------------------------------------------//
/*!
 * \brief Set abcissas.
 */
void Random_Walk_Sampling_Tables::set_abcissas()
{
    // abcissas for prob_exit sampling table
    a[0]  = 0.0; 
    a[1]  = 0.01; 
    a[2]  = 0.02; 
    a[3]  = 0.03; 
    a[4]  = 0.04;
    a[5]  = 0.05; 
    a[6]  = 0.06;
    a[7]  = 0.07; 
    a[8]  = 0.08;
    a[9]  = 0.09; 
    a[10] = 0.10;
    a[11] = 0.11;
    a[12] = 0.12;
    a[13] = 0.13;
    a[14] = 0.14;
    a[15] = 0.15; 
    a[16] = 0.16; 
    a[17] = 0.17;
    a[18] = 0.18;
    a[19] = 0.19; 
    a[20] = 0.20;
    a[21] = 0.21;
    a[22] = 0.22;
    a[23] = 0.23; 
    a[24] = 0.24;
    a[25] = 0.25;
    a[26] = 0.26; 
    a[27] = 0.27;
    a[28] = 0.28;
    a[29] = 0.29;
    a[30] = 0.30;
    a[31] = 0.35;
    a[32] = 0.40;
    a[33] = 0.45; 
    a[34] = 0.50;
    a[35] = 0.60;
    a[36] = 0.70;
    a[37] = 0.80;
    a[38] = 1.00;
    a[39] = 1.20; 
    a[40] = 20.0;

    // abcissas for prob_R_cen sampling table, we only use the first 33
    // points because the curves are so flat afterward
    for (int i = 0; i < 33; i++)
	a_R[i] = a[i];
    a_R[33] = 20.0;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Set probability of exiting sphere.
 */
void Random_Walk_Sampling_Tables::set_prob_exit()
{
    // set table for probability of escaping the sphere
    prob_exit[0]  = 0.00000E+00;
    prob_exit[1]  = 0.15671E-09;
    prob_exit[2]  = 0.29734E-04;
    prob_exit[3]  = 0.15659E-02;
    prob_exit[4]  = 0.10891E-01;
    prob_exit[5]  = 0.34001E-01;
    prob_exit[6]  = 0.71420E-01;
    prob_exit[7]  = 0.11991E+00;
    prob_exit[8]  = 0.17528E+00;
    prob_exit[9]  = 0.23386E+00;
    prob_exit[10] = 0.29290E+00;
    prob_exit[11] = 0.35053E+00; 
    prob_exit[12] = 0.40559E+00;
    prob_exit[13] = 0.45741E+00;
    prob_exit[14] = 0.50567E+00;
    prob_exit[15] = 0.55028E+00;
    prob_exit[16] = 0.59131E+00;
    prob_exit[17] = 0.62888E+00;
    prob_exit[18] = 0.66319E+00;
    prob_exit[19] = 0.69446E+00;
    prob_exit[20] = 0.72292E+00;
    prob_exit[21] = 0.74879E+00;
    prob_exit[22] = 0.77228E+00;
    prob_exit[23] = 0.79361E+00;
    prob_exit[24] = 0.81295E+00;
    prob_exit[25] = 0.83049E+00;
    prob_exit[26] = 0.84640E+00;
    prob_exit[27] = 0.86082E+00;
    prob_exit[28] = 0.87389E+00;
    prob_exit[29] = 0.88573E+00;
    prob_exit[30] = 0.89647E+00;
    prob_exit[31] = 0.93679E+00;
    prob_exit[32] = 0.96141E+00;
    prob_exit[33] = 0.97644E+00;
    prob_exit[34] = 0.98562E+00;
    prob_exit[35] = 0.99464E+00;
    prob_exit[36] = 0.99800E+00;
    prob_exit[37] = 0.99926E+00;
    prob_exit[38] = 0.99990E+00;
    prob_exit[39] = 0.99999E+00;
    prob_exit[40] = 1.000000000;
}

//---------------------------------------------------------------------------//
/*!
 * \brief Set probability of staying in sphere as a function of R1.
 */
void Random_Walk_Sampling_Tables::set_prob_R_cen()
{
    // set table for probability of remaining in sphere as a function of
    // a = Dtcen/Ro^2 and b = R1/Ro 

    // tabulated ratios, b
    b[0]  = 0.0;
    b[1]  = 0.05;
    b[2]  = 0.1;
    b[3]  = 0.15;
    b[4]  = 0.2;
    b[5]  = 0.3;
    b[6]  = 0.4;
    b[7]  = 0.5;
    b[8]  = 0.6;
    b[9]  = 0.7;
    b[10] = 0.8;
    b[11] = 0.9;
    b[12] = 1.0;

    // b = 0.0
    for (int i = 0; i < 34; i++)
	prob_R_cen[i][0] = 0.0;

    // b = 0.05
    prob_R_cen[0][1]   = 1.0;
    prob_R_cen[1][1]   = 0.0113228587;
    prob_R_cen[2][1]   = 0.00407871456;
    prob_R_cen[3][1]   = 0.00223748624;
    prob_R_cen[4][1]   = 0.00147157375;
    prob_R_cen[5][1]   = 0.00108018318;
    prob_R_cen[6][1]   = 0.000855900244;
    prob_R_cen[7][1]   = 0.000717246578;
    prob_R_cen[8][1]   = 0.000626799424;
    prob_R_cen[9][1]   = 0.000565484961;
    prob_R_cen[10][1]  = 0.000522766786;
    prob_R_cen[11][1]  = 0.000492424932;
    prob_R_cen[12][1]  = 0.000470580219;
    prob_R_cen[13][1]  = 0.000454703486;
    prob_R_cen[14][1]  = 0.000443087503;
    prob_R_cen[15][1]  = 0.000434549133;
    prob_R_cen[16][1]  = 0.00042825228;
    prob_R_cen[17][1]  = 0.000423597605;
    prob_R_cen[18][1]  = 0.00042015109;
    prob_R_cen[19][1]  = 0.000417596084;
    prob_R_cen[20][1]  = 0.000415700335;
    prob_R_cen[21][1]  = 0.000414292856;
    prob_R_cen[22][1]  = 0.000413247408;
    prob_R_cen[23][1]  = 0.000412470611;
    prob_R_cen[24][1]  = 0.000411893289;
    prob_R_cen[25][1]  = 0.000411464139;
    prob_R_cen[26][1]  = 0.000411145091;
    prob_R_cen[27][1]  = 0.000410907874;
    prob_R_cen[28][1]  = 0.000410731487;
    prob_R_cen[29][1]  = 0.000410600324;
    prob_R_cen[30][1]  = 0.000410502787;
    prob_R_cen[31][1]  = 0.000410284171;
    prob_R_cen[32][1]  = 0.000410234435;
    prob_R_cen[33][1]  = 0.000410219785;

    // b = 0.1
    prob_R_cen[0][2]   = 1.0;
    prob_R_cen[1][2]   = 0.0811085941;
    prob_R_cen[2][2]   = 0.0308605157;
    prob_R_cen[3][2]   = 0.0172449327;
    prob_R_cen[4][2]   = 0.0114475386;
    prob_R_cen[5][2]   = 0.00844988307;
    prob_R_cen[6][2]   = 0.00672038305;
    prob_R_cen[7][2]   = 0.00564669417;
    prob_R_cen[8][2]   = 0.00494442854;
    prob_R_cen[9][2]   = 0.00446754522;
    prob_R_cen[10][2]  = 0.00413493507;
    prob_R_cen[11][2]  = 0.00389852585;
    prob_R_cen[12][2]  = 0.0037282487;
    prob_R_cen[13][2]  = 0.003604458;
    prob_R_cen[14][2]  = 0.00351387324;
    prob_R_cen[15][2]  = 0.00344728175;
    prob_R_cen[16][2]  = 0.00339816895;
    prob_R_cen[17][2]  = 0.00336186303;
    prob_R_cen[18][2]  = 0.00333497999;
    prob_R_cen[19][2]  = 0.00331505049;
    prob_R_cen[20][2]  = 0.00330026317;
    prob_R_cen[21][2]  = 0.00328928443;
    prob_R_cen[22][2]  = 0.0032811296;
    prob_R_cen[23][2]  = 0.00327507033;
    prob_R_cen[24][2]  = 0.00327056701;
    prob_R_cen[25][2]  = 0.0032672195;
    prob_R_cen[26][2]  = 0.00326473081;
    prob_R_cen[27][2]  = 0.00326288043;
    prob_R_cen[28][2]  = 0.00326150455;
    prob_R_cen[29][2]  = 0.00326048144;
    prob_R_cen[30][2]  = 0.00325972061;
    prob_R_cen[31][2]  = 0.00325801533;
    prob_R_cen[32][2]  = 0.00325762737;
    prob_R_cen[33][2]  = 0.00325751309;
    
    // b = 0.15
    prob_R_cen[0][3]   = 1.0;
    prob_R_cen[1][3]   = 0.228957359;
    prob_R_cen[2][3]   = 0.095041986;
    prob_R_cen[3][3]   = 0.054728625;
    prob_R_cen[4][3]   = 0.0368863559;
    prob_R_cen[5][3]   = 0.0274783216;
    prob_R_cen[6][3]   = 0.0219887621;
    prob_R_cen[7][3]   = 0.0185570195;
    prob_R_cen[8][3]   = 0.0163024989;
    prob_R_cen[9][3]   = 0.0147672101;
    prob_R_cen[10][3]  = 0.0136944728;
    prob_R_cen[11][3]  = 0.0129311386;
    prob_R_cen[12][3]  = 0.0123809455;
    prob_R_cen[13][3]  = 0.011980781;
    prob_R_cen[14][3]  = 0.0116878776;
    prob_R_cen[15][3]  = 0.0114725195;
    prob_R_cen[16][3]  = 0.0113136714;
    prob_R_cen[17][3]  = 0.0111962378;
    prob_R_cen[18][3]  = 0.0111092796;
    prob_R_cen[19][3]  = 0.0110448125;
    prob_R_cen[20][3]  = 0.0109969784;
    prob_R_cen[21][3]  = 0.0109614639;
    prob_R_cen[22][3]  = 0.0109350843;
    prob_R_cen[23][3]  = 0.0109154833;
    prob_R_cen[24][3]  = 0.0109009157;
    prob_R_cen[25][3]  = 0.0108900869;
    prob_R_cen[26][3]  = 0.0108820363;
    prob_R_cen[27][3]  = 0.0108760505;
    prob_R_cen[28][3]  = 0.0108715997;
    prob_R_cen[29][3]  = 0.0108682901;
    prob_R_cen[30][3]  = 0.0108658289;
    prob_R_cen[31][3]  = 0.0108603125;
    prob_R_cen[32][3]  = 0.0108590575;
    prob_R_cen[33][3]  = 0.0108586878;

    // b = 0.2
    prob_R_cen[0][4]   = 1.0;
    prob_R_cen[1][4]   = 0.427593317;
    prob_R_cen[2][4]   = 0.198753966;
    prob_R_cen[3][4]   = 0.119171781;
    prob_R_cen[4][4]   = 0.0820017087;
    prob_R_cen[5][4]   = 0.0618608321;
    prob_R_cen[6][4]   = 0.0499228881;
    prob_R_cen[7][4]   = 0.0423877049;
    prob_R_cen[8][4]   = 0.0374070921;
    prob_R_cen[9][4]   = 0.0340021706;
    prob_R_cen[10][4]  = 0.0316172037;
    prob_R_cen[11][4]  = 0.0299174762;
    prob_R_cen[12][4]  = 0.028691161;
    prob_R_cen[13][4]  = 0.0277987023;
    prob_R_cen[14][4]  = 0.0271452164;
    prob_R_cen[15][4]  = 0.0266646284;
    prob_R_cen[16][4]  = 0.0263100966;
    prob_R_cen[17][4]  = 0.0260479748;
    prob_R_cen[18][4]  = 0.0258538664;
    prob_R_cen[19][4]  = 0.0257099579;
    prob_R_cen[20][4]  = 0.0256031768;
    prob_R_cen[21][4]  = 0.0255238962;
    prob_R_cen[22][4]  = 0.0254650071;
    prob_R_cen[23][4]  = 0.0254212505;
    prob_R_cen[24][4]  = 0.02538873;
    prob_R_cen[25][4]  = 0.0253645559;
    prob_R_cen[26][4]  = 0.0253465839;
    prob_R_cen[27][4]  = 0.0253332214;
    prob_R_cen[28][4]  = 0.0253232854;
    prob_R_cen[29][4]  = 0.025315897;
    prob_R_cen[30][4]  = 0.0253104026;
    prob_R_cen[31][4]  = 0.025298088;
    prob_R_cen[32][4]  = 0.0252952863;
    prob_R_cen[33][4]  = 0.025294461;

    // b = 0.3
    prob_R_cen[0][5]   = 1.0;
    prob_R_cen[1][5]   = 0.787709754;
    prob_R_cen[2][5]   = 0.477847063;
    prob_R_cen[3][5]   = 0.318228029;
    prob_R_cen[4][5]   = 0.231478479;
    prob_R_cen[5][5]   = 0.180716508;
    prob_R_cen[6][5]   = 0.149272837;
    prob_R_cen[7][5]   = 0.128884069;
    prob_R_cen[8][5]   = 0.115177586;
    prob_R_cen[9][5]   = 0.105706503;
    prob_R_cen[10][5]  = 0.0990275313;
    prob_R_cen[11][5]  = 0.094247312;
    prob_R_cen[12][5]  = 0.0907893636;
    prob_R_cen[13][5]  = 0.0882686924;
    prob_R_cen[14][5]  = 0.0864211104;
    prob_R_cen[15][5]  = 0.0850615116;
    prob_R_cen[16][5]  = 0.0840581461;
    prob_R_cen[17][5]  = 0.0833161376;
    prob_R_cen[18][5]  = 0.0827665812;
    prob_R_cen[19][5]  = 0.082359114;
    prob_R_cen[20][5]  = 0.0820567542;
    prob_R_cen[21][5]  = 0.0818322572;
    prob_R_cen[22][5]  = 0.0816654991;
    prob_R_cen[23][5]  = 0.0815415905;
    prob_R_cen[24][5]  = 0.0814494992;
    prob_R_cen[25][5]  = 0.081381043;
    prob_R_cen[26][5]  = 0.0813301497;
    prob_R_cen[27][5]  = 0.0812923095;
    prob_R_cen[28][5]  = 0.0812641727;
    prob_R_cen[29][5]  = 0.08124325;
    prob_R_cen[30][5]  = 0.081227691;
    prob_R_cen[31][5]  = 0.081192818;
    prob_R_cen[32][5]  = 0.0811848841;
    prob_R_cen[33][5]  = 0.0811825472;

    // b = 0.4
    prob_R_cen[0][6]   = 1.0;
    prob_R_cen[1][6]   = 0.953988303;
    prob_R_cen[2][6]   = 0.738557853;
    prob_R_cen[3][6]   = 0.554947337;
    prob_R_cen[4][6]   = 0.432301585;
    prob_R_cen[5][6]   = 0.352596626;
    prob_R_cen[6][6]   = 0.300188451;
    prob_R_cen[7][6]   = 0.264959427;
    prob_R_cen[8][6]   = 0.240741646;
    prob_R_cen[9][6]   = 0.223772057;
    prob_R_cen[10][6]  = 0.211700291;
    prob_R_cen[11][6]  = 0.20301331;
    prob_R_cen[12][6]  = 0.19670805;
    prob_R_cen[13][6]  = 0.192102253;
    prob_R_cen[14][6]  = 0.188722;
    prob_R_cen[15][6]  = 0.186232579;
    prob_R_cen[16][6]  = 0.184394531;
    prob_R_cen[17][6]  = 0.183034856;
    prob_R_cen[18][6]  = 0.182027653;
    prob_R_cen[19][6]  = 0.181280782;
    prob_R_cen[20][6]  = 0.180726531;
    prob_R_cen[21][6]  = 0.180314992;
    prob_R_cen[22][6]  = 0.18000929;
    prob_R_cen[23][6]  = 0.179782136;
    prob_R_cen[24][6]  = 0.17961331;
    prob_R_cen[25][6]  = 0.179487812;
    prob_R_cen[26][6]  = 0.17939451;
    prob_R_cen[27][6]  = 0.179325139;
    prob_R_cen[28][6]  = 0.179273556;
    prob_R_cen[29][6]  = 0.179235199;
    prob_R_cen[30][6]  = 0.179206675;
    prob_R_cen[31][6]  = 0.179142743;
    prob_R_cen[32][6]  = 0.179128198;
    prob_R_cen[33][6]  = 0.179123914;

    // b = 0.5
    prob_R_cen[0][7]   = 1.0;
    prob_R_cen[1][7]   = 0.994147338;
    prob_R_cen[2][7]   = 0.899965933;
    prob_R_cen[3][7]   = 0.757165225;
    prob_R_cen[4][7]   = 0.634154594;
    prob_R_cen[5][7]   = 0.543162944;
    prob_R_cen[6][7]   = 0.478822805;
    prob_R_cen[7][7]   = 0.433668326;
    prob_R_cen[8][7]   = 0.40180219;
    prob_R_cen[9][7]   = 0.379109921;
    prob_R_cen[10][7]  = 0.362805512;
    prob_R_cen[11][7]  = 0.351000304;
    prob_R_cen[12][7]  = 0.342399216;
    prob_R_cen[13][7]  = 0.336101724;
    prob_R_cen[14][7]  = 0.331473294;
    prob_R_cen[15][7]  = 0.328061644;
    prob_R_cen[16][7]  = 0.32554132;
    prob_R_cen[17][7]  = 0.323676325;
    prob_R_cen[18][7]  = 0.322294518;
    prob_R_cen[19][7]  = 0.321269741;
    prob_R_cen[20][7]  = 0.320509199;
    prob_R_cen[21][7]  = 0.319944461;
    prob_R_cen[22][7]  = 0.319524946;
    prob_R_cen[23][7]  = 0.319213218;
    prob_R_cen[24][7]  = 0.318981531;
    prob_R_cen[25][7]  = 0.318809305;
    prob_R_cen[26][7]  = 0.318681263;
    prob_R_cen[27][7]  = 0.31858606;
    prob_R_cen[28][7]  = 0.318515271;
    prob_R_cen[29][7]  = 0.318462631;
    prob_R_cen[30][7]  = 0.318423486;
    prob_R_cen[31][7]  = 0.318335749;
    prob_R_cen[32][7]  = 0.318315788;
    prob_R_cen[33][7]  = 0.318309908;

    // b = 0.6
    prob_R_cen[0][8]   = 1.0;
    prob_R_cen[1][8]   = 0.999560151;
    prob_R_cen[2][8]   = 0.97073799;
    prob_R_cen[3][8]   = 0.889783001;
    prob_R_cen[4][8]   = 0.796376066;
    prob_R_cen[5][8]   = 0.716257262;
    prob_R_cen[6][8]   = 0.654801402;
    prob_R_cen[7][8]   = 0.609591563;
    prob_R_cen[8][8]   = 0.576779761;
    prob_R_cen[9][8]   = 0.553015334;
    prob_R_cen[10][8]  = 0.535763835;
    prob_R_cen[11][8]  = 0.523194052;
    prob_R_cen[12][8]  = 0.514000577;
    prob_R_cen[13][8]  = 0.507253465;
    prob_R_cen[14][8]  = 0.502287415;
    prob_R_cen[15][8]  = 0.498623666;
    prob_R_cen[16][8]  = 0.495915641;
    prob_R_cen[17][8]  = 0.493911087;
    prob_R_cen[18][8]  = 0.492425578;
    prob_R_cen[19][8]  = 0.491323757;
    prob_R_cen[20][8]  = 0.490505976;
    prob_R_cen[21][8]  = 0.489898706;
    prob_R_cen[22][8]  = 0.489447584;
    prob_R_cen[23][8]  = 0.489112364;
    prob_R_cen[24][8]  = 0.488863214;
    prob_R_cen[25][8]  = 0.488678005;
    prob_R_cen[26][8]  = 0.488540311;
    prob_R_cen[27][8]  = 0.488437932;
    prob_R_cen[28][8]  = 0.488361806;
    prob_R_cen[29][8]  = 0.488305198;
    prob_R_cen[30][8]  = 0.488263101;
    prob_R_cen[31][8]  = 0.488168749;
    prob_R_cen[32][8]  = 0.488147284;
    prob_R_cen[33][8]  = 0.488140961;

    // b = 0.7
    prob_R_cen[0][9]   = 1.0;
    prob_R_cen[1][9]   = 0.999980359;
    prob_R_cen[2][9]   = 0.993455501;
    prob_R_cen[3][9]   = 0.958809891;
    prob_R_cen[4][9]   = 0.904115754;
    prob_R_cen[5][9]   = 0.849270232;
    prob_R_cen[6][9]   = 0.803548923;
    prob_R_cen[7][9]   = 0.768307125;
    prob_R_cen[8][9]   = 0.74202843;
    prob_R_cen[9][9]   = 0.722688465;
    prob_R_cen[10][9]  = 0.708513334;
    prob_R_cen[11][9]  = 0.698124868;
    prob_R_cen[12][9]  = 0.69049991;
    prob_R_cen[13][9]  = 0.684891872;
    prob_R_cen[14][9]  = 0.680758779;
    prob_R_cen[15][9]  = 0.6777071;
    prob_R_cen[16][9]  = 0.675450373;
    prob_R_cen[17][9]  = 0.673779381;
    prob_R_cen[18][9]  = 0.672540836;
    prob_R_cen[19][9]  = 0.671622088;
    prob_R_cen[20][9]  = 0.670940139;
    prob_R_cen[21][9]  = 0.670433714;
    prob_R_cen[22][9]  = 0.670057498;
    prob_R_cen[23][9]  = 0.669777934;
    prob_R_cen[24][9]  = 0.669570148;
    prob_R_cen[25][9]  = 0.669415687;
    prob_R_cen[26][9]  = 0.669300852;
    prob_R_cen[27][9]  = 0.669215469;
    prob_R_cen[28][9]  = 0.66915198;
    prob_R_cen[29][9]  = 0.66910477;
    prob_R_cen[30][9]  = 0.669069662;
    prob_R_cen[31][9]  = 0.668990973;
    prob_R_cen[32][9]  = 0.668973071;
    prob_R_cen[33][9]  = 0.668967797;

    // b = 0.8
    prob_R_cen[0][10]  = 1.0;
    prob_R_cen[1][10]  = 0.999999477;
    prob_R_cen[2][10]  = 0.998895671;
    prob_R_cen[3][10]  = 0.987859987;
    prob_R_cen[4][10]  = 0.964233765;
    prob_R_cen[5][10]  = 0.936802702;
    prob_R_cen[6][10]  = 0.912156704;
    prob_R_cen[7][10]  = 0.892373691;
    prob_R_cen[8][10]  = 0.877280699;
    prob_R_cen[9][10]  = 0.866024321;
    prob_R_cen[10][10] = 0.857708846;
    prob_R_cen[11][10] = 0.851585906;
    prob_R_cen[12][10] = 0.847078942;
    prob_R_cen[13][10] = 0.843758397;
    prob_R_cen[14][10] = 0.84130859;
    prob_R_cen[15][10] = 0.839498606;
    prob_R_cen[16][10] = 0.838159591;
    prob_R_cen[17][10] = 0.837167879;
    prob_R_cen[18][10] = 0.836432712;
    prob_R_cen[19][10] = 0.835887319;
    prob_R_cen[20][10] = 0.835482473;
    prob_R_cen[21][10] = 0.835181819;
    prob_R_cen[22][10] = 0.834958463;
    prob_R_cen[23][10] = 0.834792486;
    prob_R_cen[24][10] = 0.834669123;
    prob_R_cen[25][10] = 0.834577419;
    prob_R_cen[26][10] = 0.83450924;
    prob_R_cen[27][10] = 0.834458548;
    prob_R_cen[28][10] = 0.834420854;
    prob_R_cen[29][10] = 0.834392824;
    prob_R_cen[30][10] = 0.834371981;
    prob_R_cen[31][10] = 0.834325262;
    prob_R_cen[32][10] = 0.834314633;
    prob_R_cen[33][10] = 0.834311503;

    // b = 0.9
    prob_R_cen[0][11]  = 1.0;
    prob_R_cen[1][11]  = 0.999999992;
    prob_R_cen[2][11]  = 0.99987815;
    prob_R_cen[3][11]  = 0.997775965;
    prob_R_cen[4][11]  = 0.992053435;
    prob_R_cen[5][11]  = 0.984639536;
    prob_R_cen[6][11]  = 0.977611849;
    prob_R_cen[7][11]  = 0.971810265;
    prob_R_cen[8][11]  = 0.967314996;
    prob_R_cen[9][11]  = 0.963932581;
    prob_R_cen[10][11] = 0.961420866;
    prob_R_cen[11][11] = 0.959565685;
    prob_R_cen[12][11] = 0.958197584;
    prob_R_cen[13][11] = 0.957188486;
    prob_R_cen[14][11] = 0.956443492;
    prob_R_cen[15][11] = 0.95589284;
    prob_R_cen[16][11] = 0.955485367;
    prob_R_cen[17][11] = 0.955183534;
    prob_R_cen[18][11] = 0.954959761;
    prob_R_cen[19][11] = 0.954793742;
    prob_R_cen[20][11] = 0.954670502;
    prob_R_cen[21][11] = 0.954578977;
    prob_R_cen[22][11] = 0.954510983;
    prob_R_cen[23][11] = 0.954460455;
    prob_R_cen[24][11] = 0.9544229;
    prob_R_cen[25][11] = 0.954394983;
    prob_R_cen[26][11] = 0.954374227;
    prob_R_cen[27][11] = 0.954358795;
    prob_R_cen[28][11] = 0.95434732;
    prob_R_cen[29][11] = 0.954338787;
    prob_R_cen[30][11] = 0.954332441;
    prob_R_cen[31][11] = 0.954318219;
    prob_R_cen[32][11] = 0.954314983;
    prob_R_cen[33][11] = 0.95431403;
	
    // b = 1.0
    for (int i = 0; i < 34; i++)
	prob_R_cen[i][12] = 1.0;
}

} // end namespace rtt_imc

//---------------------------------------------------------------------------//
//                              end of Random_Walk.cc
//---------------------------------------------------------------------------//
