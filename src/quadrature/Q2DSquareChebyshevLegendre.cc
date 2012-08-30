//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   quadrature/Q2DSquareChebyshevLegendre.cc
 * \author James Warsa
 * \date   Fri Jun  9 13:52:25 2006
 * \brief  
 * \note   Copyright 2006 Los Alamos National Security, LLC.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>

#include "ds++/Soft_Equivalence.hh"
#include "units/PhysicalConstants.hh"

#include "Q1DGaussLeg.hh"
#include "Q2DSquareChebyshevLegendre.hh"
#include "Ordinate.hh"

namespace rtt_quadrature
{
//---------------------------------------------------------------------------//
/*!
 * \brief Constructs a 2D Square Chebyshev Legendre quadrature object.
 *
 * \param snOrder_ Integer specifying the order of the SN set to be
 *                 constructed.  Number of ordinates = (snOrder*snOrder)
 * \param norm_    A normalization constant.  The sum of the quadrature
 *                 weights will be equal to this value (default = 2*PI).
 */
//---------------------------------------------------------------------------//
Q2DSquareChebyshevLegendre::Q2DSquareChebyshevLegendre( size_t sn_order_, double norm_, Quadrature::QIM qm_ ) 
    : Quadrature( sn_order_, norm_, qm_ ), numOrdinates (sn_order_ * sn_order_)
{
    using std::fabs;
    using std::sqrt;
    using std::cos;
    using rtt_dsxx::soft_equiv;
    
    Require ( snOrder > 0 );
    Require ( norm > 0.0 );
    Insist( snOrder%2 == 0, "SN order must be even");

    // Force the direction vectors to be the correct length.
    mu.resize(numOrdinates);
    xi.resize(numOrdinates);
    wt.resize(numOrdinates);

//     if ( snOrder == 8)
//     {
//         eta[ 0] =   -9.602898564975362317E-01;
//         eta[ 1] =   -9.602898564975362317E-01;
//         eta[ 2] =   -9.602898564975362317E-01;
//         eta[ 3] =   -9.602898564975362317E-01;
//         eta[ 4] =   -7.966664774136267396E-01;
//         eta[ 5] =   -7.966664774136267396E-01;
//         eta[ 6] =   -7.966664774136267396E-01;
//         eta[ 7] =   -7.966664774136267396E-01;
//         eta[ 8] =   -5.255324099163289858E-01;
//         eta[ 9] =   -5.255324099163289858E-01;
//         eta[10] =   -5.255324099163289858E-01;
//         eta[11] =   -5.255324099163289858E-01;
//         eta[12] =   -1.834346424956498049E-01;
//         eta[13] =   -1.834346424956498049E-01;
//         eta[14] =   -1.834346424956498049E-01;
//         eta[15] =   -1.834346424956498049E-01;
//         eta[16] =   -9.602898564975362317E-01;
//         eta[17] =   -9.602898564975362317E-01;
//         eta[18] =   -9.602898564975362317E-01;
//         eta[19] =   -9.602898564975362317E-01;
//         eta[20] =   -7.966664774136267396E-01;
//         eta[21] =   -7.966664774136267396E-01;
//         eta[22] =   -7.966664774136267396E-01;
//         eta[23] =   -7.966664774136267396E-01;
//         eta[24] =   -5.255324099163289858E-01;
//         eta[25] =   -5.255324099163289858E-01;
//         eta[26] =   -5.255324099163289858E-01;
//         eta[27] =   -5.255324099163289858E-01;
//         eta[28] =   -1.834346424956498049E-01;
//         eta[29] =   -1.834346424956498049E-01;
//         eta[30] =   -1.834346424956498049E-01;
//         eta[31] =   -1.834346424956498049E-01;
//         eta[32] =    9.602898564975362317E-01;
//         eta[33] =    9.602898564975362317E-01;
//         eta[34] =    9.602898564975362317E-01;
//         eta[35] =    9.602898564975362317E-01;
//         eta[36] =    7.966664774136267396E-01;
//         eta[37] =    7.966664774136267396E-01;
//         eta[38] =    7.966664774136267396E-01;
//         eta[39] =    7.966664774136267396E-01;
//         eta[40] =    5.255324099163289858E-01;
//         eta[41] =    5.255324099163289858E-01;
//         eta[42] =    5.255324099163289858E-01;
//         eta[43] =    5.255324099163289858E-01;
//         eta[44] =    1.834346424956498049E-01;
//         eta[45] =    1.834346424956498049E-01;
//         eta[46] =    1.834346424956498049E-01;
//         eta[47] =    1.834346424956498049E-01;
//         eta[48] =    9.602898564975362317E-01;
//         eta[49] =    9.602898564975362317E-01;
//         eta[50] =    9.602898564975362317E-01;
//         eta[51] =    9.602898564975362317E-01;
//         eta[52] =    7.966664774136267396E-01;
//         eta[53] =    7.966664774136267396E-01;
//         eta[54] =    7.966664774136267396E-01;
//         eta[55] =    7.966664774136267396E-01;
//         eta[56] =    5.255324099163289858E-01;
//         eta[57] =    5.255324099163289858E-01;
//         eta[58] =    5.255324099163289858E-01;
//         eta[59] =    5.255324099163289858E-01;
//         eta[60] =    1.834346424956498049E-01;
//         eta[61] =    1.834346424956498049E-01;
//         eta[62] =    1.834346424956498049E-01;
//         eta[63] =    1.834346424956498049E-01;
//         mu[ 0] =   -2.736432967052099398E-01;
//         mu[ 1] =   -2.319835853645034495E-01;
//         mu[ 2] =   -1.550064760884893226E-01;
//         mu[ 3] =   -5.443103596520747900E-02;
//         mu[ 4] =   -5.928054175855475004E-01;
//         mu[ 5] =   -5.025561665526400068E-01;
//         mu[ 6] =   -3.357972948450873296E-01;
//         mu[ 7] =   -1.179163290074280016E-01;
//         mu[ 8] =   -8.344262052007293529E-01;
//         mu[ 9] =   -7.073923795513043837E-01;
//         mu[10] =   -4.726644766430822852E-01;
//         mu[11] =   -1.659776918801013757E-01;
//         mu[12] =   -9.641432254265330918E-01;
//         mu[13] =   -8.173611593354460190E-01;
//         mu[14] =   -5.461432661328973778E-01;
//         mu[15] =   -1.917800114626509772E-01;
//         mu[16] =    2.736432967052099398E-01;
//         mu[17] =    2.319835853645034495E-01;
//         mu[18] =    1.550064760884893226E-01;
//         mu[19] =    5.443103596520747900E-02;
//         mu[20] =    5.928054175855475004E-01;
//         mu[21] =    5.025561665526400068E-01;
//         mu[22] =    3.357972948450873296E-01;
//         mu[23] =    1.179163290074280016E-01;
//         mu[24] =    8.344262052007293529E-01;
//         mu[25] =    7.073923795513043837E-01;
//         mu[26] =    4.726644766430822852E-01;
//         mu[27] =    1.659776918801013757E-01;
//         mu[28] =    9.641432254265330918E-01;
//         mu[29] =    8.173611593354460190E-01;
//         mu[30] =    5.461432661328973778E-01;
//         mu[31] =    1.917800114626509772E-01;
//         mu[32] =   -2.736432967052099398E-01;
//         mu[33] =   -2.319835853645034495E-01;
//         mu[34] =   -1.550064760884893226E-01;
//         mu[35] =   -5.443103596520747900E-02;
//         mu[36] =   -5.928054175855475004E-01;
//         mu[37] =   -5.025561665526400068E-01;
//         mu[38] =   -3.357972948450873296E-01;
//         mu[39] =   -1.179163290074280016E-01;
//         mu[40] =   -8.344262052007293529E-01;
//         mu[41] =   -7.073923795513043837E-01;
//         mu[42] =   -4.726644766430822852E-01;
//         mu[43] =   -1.659776918801013757E-01;
//         mu[44] =   -9.641432254265330918E-01;
//         mu[45] =   -8.173611593354460190E-01;
//         mu[46] =   -5.461432661328973778E-01;
//         mu[47] =   -1.917800114626509772E-01;
//         mu[48] =    2.736432967052099398E-01;
//         mu[49] =    2.319835853645034495E-01;
//         mu[50] =    1.550064760884893226E-01;
//         mu[51] =    5.443103596520747900E-02;
//         mu[52] =    5.928054175855475004E-01;
//         mu[53] =    5.025561665526400068E-01;
//         mu[54] =    3.357972948450873296E-01;
//         mu[55] =    1.179163290074280016E-01;
//         mu[56] =    8.344262052007293529E-01;
//         mu[57] =    7.073923795513043837E-01;
//         mu[58] =    4.726644766430822852E-01;
//         mu[59] =    1.659776918801013757E-01;
//         mu[60] =    9.641432254265330918E-01;
//         mu[61] =    8.173611593354460190E-01;
//         mu[62] =    5.461432661328973778E-01;
//         mu[63] =    1.917800114626509772E-01;
//         wt[ 0] =    6.326783518148516197E-03;
//         wt[ 1] =    6.326783518148516197E-03;
//         wt[ 2] =    6.326783518148516197E-03;
//         wt[ 3] =    6.326783518148516197E-03;
//         wt[ 4] =    1.389881465333590441E-02;
//         wt[ 5] =    1.389881465333590441E-02;
//         wt[ 6] =    1.389881465333590441E-02;
//         wt[ 7] =    1.389881465333590441E-02;
//         wt[ 8] =    1.960666536736795546E-02;
//         wt[ 9] =    1.960666536736795546E-02;
//         wt[10] =    1.960666536736795546E-02;
//         wt[11] =    1.960666536736795546E-02;
//         wt[12] =    2.266773646114762394E-02;
//         wt[13] =    2.266773646114762394E-02;
//         wt[14] =    2.266773646114762394E-02;
//         wt[15] =    2.266773646114762394E-02;
//         wt[16] =    6.326783518148516197E-03;
//         wt[17] =    6.326783518148516197E-03;
//         wt[18] =    6.326783518148516197E-03;
//         wt[19] =    6.326783518148516197E-03;
//         wt[20] =    1.389881465333590441E-02;
//         wt[21] =    1.389881465333590441E-02;
//         wt[22] =    1.389881465333590441E-02;
//         wt[23] =    1.389881465333590441E-02;
//         wt[24] =    1.960666536736795546E-02;
//         wt[25] =    1.960666536736795546E-02;
//         wt[26] =    1.960666536736795546E-02;
//         wt[27] =    1.960666536736795546E-02;
//         wt[28] =    2.266773646114762394E-02;
//         wt[29] =    2.266773646114762394E-02;
//         wt[30] =    2.266773646114762394E-02;
//         wt[31] =    2.266773646114762394E-02;
//         wt[32] =    6.326783518148516197E-03;
//         wt[33] =    6.326783518148516197E-03;
//         wt[34] =    6.326783518148516197E-03;
//         wt[35] =    6.326783518148516197E-03;
//         wt[36] =    1.389881465333590441E-02;
//         wt[37] =    1.389881465333590441E-02;
//         wt[38] =    1.389881465333590441E-02;
//         wt[39] =    1.389881465333590441E-02;
//         wt[40] =    1.960666536736795546E-02;
//         wt[41] =    1.960666536736795546E-02;
//         wt[42] =    1.960666536736795546E-02;
//         wt[43] =    1.960666536736795546E-02;
//         wt[44] =    2.266773646114762394E-02;
//         wt[45] =    2.266773646114762394E-02;
//         wt[46] =    2.266773646114762394E-02;
//         wt[47] =    2.266773646114762394E-02;
//         wt[48] =    6.326783518148516197E-03;
//         wt[49] =    6.326783518148516197E-03;
//         wt[50] =    6.326783518148516197E-03;
//         wt[51] =    6.326783518148516197E-03;
//         wt[52] =    1.389881465333590441E-02;
//         wt[53] =    1.389881465333590441E-02;
//         wt[54] =    1.389881465333590441E-02;
//         wt[55] =    1.389881465333590441E-02;
//         wt[56] =    1.960666536736795546E-02;
//         wt[57] =    1.960666536736795546E-02;
//         wt[58] =    1.960666536736795546E-02;
//         wt[59] =    1.960666536736795546E-02;
//         wt[60] =    2.266773646114762394E-02;
//         wt[61] =    2.266773646114762394E-02;
//         wt[62] =    2.266773646114762394E-02;
//         wt[63] =    2.266773646114762394E-02;
//     }
//     else
    {
        Q1DGaussLeg gauss(snOrder, 2.0, interpModel);    
        
        // NOTE: this aligns the gauss points with the x-axis (r-axis in cylindrical coords)
        
        for (unsigned i=0; i<snOrder; ++i)
        {
            double const xmu=gauss.getMu(i);

            double const xwt=gauss.getWt(i);
            double const xsr=sqrt(1.0-xmu*xmu);
            
            for (unsigned j=0; j<snOrder; ++j)
            {
                size_t ordinate=j+i*snOrder;
                
                xi[ordinate] = xmu;
                mu[ordinate]  = xsr*cos(rtt_units::PI*(2.0*j+1.0)/snOrder/2.0);
                wt[ordinate]  = xwt/snOrder;
            }
        }
    }

    // Normalize the quadrature set
    double wsum = 0.0;
    for(size_t ordinate = 0; ordinate < numOrdinates; ++ordinate)
	wsum = wsum + wt[ordinate];
    
    for(size_t ordinate = 0; ordinate < numOrdinates; ++ordinate)
	wt[ordinate] = wt[ordinate]*(norm/wsum);

    // Sort the directions by xi and then by mu
    sortOrdinates();
    
    // Verify that the quadrature meets our integration requirements.
    Ensure( soft_equiv(iDomega(),norm) );

    // check each component of the vector result
    vector<double> iod = iOmegaDomega();
    Ensure( soft_equiv(iod[0],0.0) );
    Ensure( soft_equiv(iod[1],0.0) );
		    
    // check each component of the tensor result
    vector<double> iood = iOmegaOmegaDomega();
    Ensure( soft_equiv(iood[0],norm/3.0) );  // mu*mu
    Ensure( soft_equiv(iood[1],0.0) ); // mu*eta
    Ensure( soft_equiv(iood[2],0.0) ); // eta*mu
    Ensure( soft_equiv(iood[3],norm/3.0) ); // eta*eta

    // Copy quadrature data { mu, eta } into the vector omega.
    omega.resize( numOrdinates );
    size_t ndims = dimensionality();
    for ( size_t ordinate = 0; ordinate < numOrdinates; ++ordinate )
    {
	omega[ordinate].resize(ndims);
	omega[ordinate][0] = mu[ordinate];
	omega[ordinate][1] = xi[ordinate];
    }

    //display();

} // end of Q2DLevelSym() constructor.

//---------------------------------------------------------------------------//
/*!
 * \brief Resort all of the ordinates by xi and then by mu.
 *
 * The ctor for OrdinateSet sorts automatically.
 */
void Q2DSquareChebyshevLegendre::sortOrdinates(void)
{
    size_t len( mu.size() );

    // temporary storage
    vector<Ordinate> omega;
    for( size_t m=0; m<len; ++m )
    {
        double eta=std::sqrt(1.0-mu[m]*mu[m]-xi[m]*xi[m]);
        omega.push_back( Ordinate(mu[m],eta,xi[m],wt[m] ) );    
    }
    
    std::sort(omega.begin(),omega.end(),Ordinate::SnCompare);
    
    // Save sorted data
    for( size_t m=0; m<len; ++m )
    {
        mu[m]=omega[m].mu();
        xi[m]=omega[m].xi();
        wt[m]=omega[m].wt();        
    }
    
    return;
}

//---------------------------------------------------------------------------//

void Q2DSquareChebyshevLegendre::display() const 
{
    using std::cout;
    using std::endl;
    using std::setprecision;

    cout << endl << "The Quadrature directions and weights are:" 
	 << endl << endl;
    cout << "   m  \t    mu        \t    xi        \t     wt      " << endl;
    cout << "  --- \t------------- \t------------- \t-------------" << endl;
    double sum_wt = 0.0;
    for ( size_t ordinate = 0; ordinate < mu.size(); ++ordinate )
    {
	cout << "   "
	     << ordinate << "\t"
	     << setprecision(10) << mu[ordinate]  << "\t"
	     << setprecision(10) << xi[ordinate] << "\t"
	     << setprecision(10) << wt[ordinate]  << endl;
	sum_wt += wt[ordinate];
    }
    cout << endl << "  The sum of the weights is " << sum_wt << endl;
    cout << endl;
}

} // end namespace rtt_quadrature

//---------------------------------------------------------------------------//
//                 end of Q2DSquareChebyshevLegendre.cc
//---------------------------------------------------------------------------//
