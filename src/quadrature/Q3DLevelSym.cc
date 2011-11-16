//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   quadrature/Q3DLevelSym.cc
 * \author Kelly Thompson
 * \date   Wed Sep  1 10:19:52 2004
 * \brief  
 * \note   Copyright 2004 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <iostream>
#include <iomanip>
#include <cmath>

#include "ds++/Soft_Equivalence.hh"
#include "units/PhysicalConstants.hh"
#include "Q3DLevelSym.hh"

namespace rtt_quadrature
{

/*!
 * \brief Constructs a 3D Level Symmetric quadrature object.
 *
 * \param snOrder_ Integer specifying the order of the SN set to be
 *                 constructed.  Number of ordinates = (snOrder+2)*snOrder.
 * \param norm_    A normalization constant.  The sum of the quadrature
 *                 weights will be equal to this value (default = 4*PI).
 */
Q3DLevelSym::Q3DLevelSym( size_t sn_order_, double norm_ ) 
    : Quadrature( sn_order_, norm_ ), 
      numOrdinates ( sn_order_ * (sn_order_+2) )

{ 
    using rtt_dsxx::soft_equiv;

    Require ( snOrder > 0 );
    Require ( norm > 0.0 );
    Require ( snOrder%2 == 0 );

    if ( snOrder > 24 )
    {
        throw std::range_error("Sn order out of range");
    }

    // The number of quadrature levels is equal to the requested SN order.
    size_t levels( snOrder );
    size_t octantOrdinates( numOrdinates/8 );  // 8 octants in 3D.

    // Force the direction vectors to be the correct length.
    mu.resize(numOrdinates);
    eta.resize(numOrdinates);
    xi.resize(numOrdinates);
    wt.resize(numOrdinates);

    // Set up some temporaries.
    vector<double> att;
    vector<double> wtt;

    // Set LS values for one quadrant based on the SN order specified.
    switch ( levels ) {
    case 2: // LQn S-2 Quadrature Set.
	att.resize(1);
	att[0] = 0.577350269189625764509149;
	wt[0]  = 1.000000000000000000000000;
	break;

    case 4: { // LQn S-4 Quadrature Set
	att.resize(2);
	wtt.resize(1);
	att[0] = 0.350021174581540677777041;
	att[1] = 0.868890300722201205229788;
	wtt[0] = 0.333333333333333333333333;
 	size_t wp4[3] = {0,0,0};
	for(size_t m=0;m<=octantOrdinates-1;++m)
	    wt[m] = wtt[wp4[m]];
	break;
    }	

    case 6: { // LQn S-6 Quadrature Set
	att.resize(3);
	wtt.resize(2);
	att[0] = 0.266635401516704720331535;
	att[1] = 0.681507726536546927403750;
	att[2] = 0.926180935517489107558380;
	wtt[0] = 0.176126130863383433783565;
	wtt[1] = 0.157207202469949899549768;
 	size_t wp6[6] = {0,1,0,1,1,0};
	for( size_t m=0;m<=octantOrdinates-1;++m)
	    wt[m] = wtt[wp6[m]];
	break;
    }

    case 8: { // LQn S-8 Quadrature Set
	att.resize(4);
	wtt.resize(3);
	att[0] = 0.218217890235992381266097;
	att[1] = 0.577350269189625764509149;
	att[2] = 0.786795792469443145800830;
	att[3] = 0.951189731211341853132399;
	wtt[0] = 0.120987654320987654320988;
	wtt[1] = 0.0907407407407407407407407;
	wtt[2] = 0.0925925925925925925925926;
 	size_t wp8[10]  = {0,1,1,0,1,2,1,1,1,0};
	for( size_t m=0;m<=octantOrdinates-1;++m)
	    wt[m] = wtt[wp8[m]];
	break;
    }
	
    case 10: {  // LQn S-10 Quadrature Set
	att.resize(5);
	wtt.resize(4);
	att[0] = 0.189321326478010476671494;
	att[1] = 0.508881755582618974382711;
	att[2] = 0.694318887594384317279217;
	att[3] = 0.839759962236684758403029;
	att[4] = 0.963490981110468484701598;
	wtt[0] = 0.0893031479843567214704325;
	wtt[1] = 0.0725291517123655242296233;
	wtt[2] = 0.0450437674364086390490892;
	wtt[3] = 0.0539281144878369243545650;
 	size_t wp10[15] = {0,1,2,1,0,1,3,3,1,2,3,2,1,1,0};
	for( size_t m=0;m<=octantOrdinates-1;++m)
	    wt[m] = wtt[wp10[m]];
	break;
    }
	
    case 12: { // LQn S-12 Quadrature Set */
	att.resize(6);
	wtt.resize(5);
	att[0] = 0.167212652822713264084504;
	att[1] = 0.459547634642594690016761;
	att[2] = 0.628019096642130901034766;
	att[3] = 0.760021014833664062877138;
	att[4] = 0.872270543025721502340662;
	att[5] = 0.971637719251358378302376;
	wtt[0] = 0.0707625899700910439766549;
	wtt[1] = 0.0558811015648888075828962;
	wtt[2] = 0.0373376737588285824652402;
	wtt[3] = 0.0502819010600571181385765;
	wtt[4] = 0.0258512916557503911218290;
 	size_t wp12[21] = {0,1,2,2,1,0,1,3,4,3,1,2,4,4,2,2,3,2,1,1,0};
	for( size_t m=0;m<=octantOrdinates-1;++m)
	    wt[m] = wtt[wp12[m]];
	break;
    }

    case 14: { // LQn S-14 Quadrature Set
	att.resize(7);
	wtt.resize(7);
	att[0] = 0.151985861461031912404799;
	att[1] = 0.422156982304796966896263;
	att[2] = 0.577350269189625764509149;
	att[3] = 0.698892086775901338963210;
	att[4] = 0.802226255231412057244328;
	att[5] = 0.893691098874356784901111;
	att[6] = 0.976627152925770351762946;
	wtt[0] = 0.0579970408969969964063611;
	wtt[1] = 0.0489007976368104874582568;
	wtt[2] = 0.0227935342411872473257345;
	wtt[3] = 0.0394132005950078294492985;
	wtt[4] = 0.0380990861440121712365891;
	wtt[5] = 0.0258394076418900119611012;
	wtt[6] = 0.00826957997262252825269908;
 	size_t wp14[28] = { 0,1,2,3,2,1,0,1,4,5,5,4,1,2,5,6,5,2,3,5,5,3,2,4,2,1,
			     1,0 }; 
	for( size_t m=0;m<=octantOrdinates-1;++m)
	    wt[m] = wtt[wp14[m]];
	break;
    }
	
    case 16: { // LQn S-16 Quadrature Set
	att.resize(8);
	wtt.resize(8);
	att[0] = 0.138956875067780344591732;
	att[1] = 0.392289261444811712294197;
	att[2] = 0.537096561300879079878296;
	att[3] = 0.650426450628771770509703;
	att[4] = 0.746750573614681064580018;
	att[5] = 0.831996556910044145168291;
	att[6] = 0.909285500943725291652116;
	att[7] = 0.980500879011739882135849;
	wtt[0] = 0.0489872391580385335008367;
	wtt[1] = 0.0413295978698440232405505;
	wtt[2] = 0.0203032007393652080748070;
	wtt[3] = 0.0265500757813498446015484;
	wtt[4] = 0.0379074407956004002099321;
	wtt[5] = 0.0135295047786756344371600;
	wtt[6] = 0.0326369372026850701318409;
	wtt[7] = 0.0103769578385399087825920;
 	size_t wp16[36] = { 0,1,2,3,3,2,1,0,1,4,5,6,5,4,1,2,5,7,7,5,2,3,6,7,6,3,
			     3,5,5,3,2,4,2,1,1,0 };
	for( size_t m=0;m<=octantOrdinates-1;++m)
	    wt[m] = wtt[wp16[m]];
	break;
    }
	
    case 18: { // LQn S-18 Quadrature Set
	att.resize(9);
	wtt.resize(10);
	att[0]  = 0.129344504545924818514086;
	att[1]  = 0.368043816053393605686086;
	att[2]  = 0.504165151725164054411848;
	att[3]  = 0.610662549934881101060239;
	att[4]  = 0.701166884252161909657019;
	att[5]  = 0.781256199495913171286914;
	att[6]  = 0.853866206691488372341858;
	att[7]  = 0.920768021061018932899055;
	att[8]  = 0.983127661236087115272518;
	wtt[0]  = 0.0422646448843821748535825;
	wtt[1]  = 0.0376127473827281471532380;
	wtt[2]  = 0.0122691351637405931037187;
	wtt[3]  = 0.0324188352558815048715646;
	wtt[4]  = 0.00664438614619073823264082;
	wtt[5]  = 0.0312093838436551370068864;
	wtt[6]  = 0.0160127252691940275641645;
	wtt[7]  = 0.0200484595308572875885066;
	wtt[8]  = 0.000111409402059638628382279;
	wtt[9]  = 0.0163797038522425240494567;
 	size_t wp18[45] = { 0,1,2,3,4,3,2,1,0,1,5,6,7,7,6,5,1,2,6,8,9,8,6,2,3,7, 
			     9,9,7,3,4,7,8,7,4,3,6,6,3,2,5,2,1,1,0 };
	for( size_t m=0;m<=octantOrdinates-1;++m)
	    wt[m] = wtt[wp18[m]];
	break;
    }
	
    case 20: { // LQn S-20 Quadrature Set
	att.resize(10);
	wtt.resize(12);
	att[0]  =  0.120603343036693597409418;
	att[1]  =  0.347574292315847257336779;
	att[2]  =  0.476519266143665680817278;
	att[3]  =  0.577350269189625764509149;
	att[4]  =  0.663020403653288019308783;
	att[5]  =  0.738822561910371432904974;
	att[6]  =  0.807540401661143067193530;
	att[7]  =  0.870852583760463975580977;
	att[8]  =  0.929863938955324566667817;
	att[9]  =  0.985347485558646574628509;
	wtt[0]  =  0.0370210490604481342320295;
	wtt[1]  =  0.0332842165376314841003910;
	wtt[2]  =  0.0111738965965092519614021;
	wtt[3]  =  0.0245177476959359285418987;
	wtt[4]  =  0.0135924329650041789567081;
	wtt[5]  =  0.0318029065936585971501960;
	wtt[6]  =  0.00685492401402507781062634;
	wtt[7]  =  0.0308105481755299327227893;
	wtt[8]  = -0.000139484716502602877593527;
	wtt[9]  =  0.00544675187330776223879437;
	wtt[10] =  0.00474564692642379971238396;
	wtt[11] =  0.0277298541009064049325246;
 	size_t wp20[55] = { 0,1,2,3,4,4,3,2,1,0,1,5,6,7,8,7,6,5,1,2,6,9,10,10,9,
			     6,2,3,7,10,11,10,7,3,4,8,10,10,8,4,4,7,9,7,4,3,6,6,3,
			     2,5,2,1,1,0};
	for( size_t m=0;m<=octantOrdinates-1;++m)
	    wt[m] = wtt[wp20[m]];
	break;
    }
	
    case 22: { // LQn S-22 Quadrature Set
	att.resize(11);
	wtt.resize(14);
	att[0]  =  0.113888641383070838173488;
	att[1]  =  0.330271760593086736334651;
	att[2]  =  0.452977095507524183904005;
	att[3]  =  0.548905330875560154226714;
	att[4]  =  0.630401360620980621392149;
	att[5]  =  0.702506006153654989703184;
	att[6]  =  0.767869456282208576047898;
	att[7]  =  0.828089557415325768804621;
	att[8]  =  0.884217805921983001958912;
	att[9]  =  0.936989829997455780115072;
	att[10] =  0.986944149751056870330152;
	wtt[0]  =  0.0329277718552552308051381;
	wtt[1]  =  0.0309569328165031538543025;
	wtt[2]  =  0.00577105953220643022391829;
	wtt[3]  =  0.0316834548379952775919418;
	wtt[4]  = -0.00669350304140992494103696;
	wtt[5]  =  0.0368381622687682466526634;
	wtt[6]  =  0.0273139698006629537455404;
	wtt[7]  =  0.0100962716435030437817055;
	wtt[8]  =  0.0195181067555849392224199;
	wtt[9]  =  0.0117224275470949786864925;
	wtt[10] = -0.00442773155233893239996431;
	wtt[11] =  0.0156214785078803432781324;
	wtt[12] = -0.0101774221315738297143270;
	wtt[13] =  0.0135061258938431808485310;
 	size_t wp22[66] = { 0,1,2,3,4,5,4,3,2,1,0,1,6,7,8,9,9,8,7,6,1,2,7,10,11,
			     12,11,10,7,2,3,8,11,13,13,11,8,3,4,9,12,13,12,9,4,5,
			     9,11,11,9,5,4,8,10,8,4,3,7,7,3,2,6,2,1,1,0 };
	for( size_t m=0;m<=octantOrdinates-1;++m)
	    wt[m] = wtt[wp22[m]];
	break;
    }
	
    case 24: // LQn S-24 Quadrature Set
    { 
	att.resize(12);
	wtt.resize(16);
	att[0]  =  0.107544208775147285552086;
	att[1]  =  0.315151630853896874875332;
	att[2]  =  0.432522073446742487657060;
	att[3]  =  0.524242441631224399254880;
	att[4]  =  0.602150256328323868809286;
	att[5]  =  0.671073561381361944701265;
	att[6]  =  0.733549261041044861004094;
	att[7]  =  0.791106384731321324814121;
	att[8]  =  0.844750913317919895113069;
	att[9]  =  0.895186516397704814461305;
	att[10] =  0.942928254285052510917188;
	att[11] =  0.988366574868785749937406;
	wtt[0]  =  0.0295284942030736546025272;
	wtt[1]  =  0.0281530651743695026834932;
	wtt[2]  =  0.00519730128072174996473824;
	wtt[3]  =  0.0259897467786242920448933;
	wtt[4]  =  0.00146378160153344429844948;
	wtt[5]  =  0.0166609651269037212368055;
	wtt[6]  =  0.0281343344093849194875108;
	wtt[7]  =  0.00214364311909247909952968;
	wtt[8]  =  0.0331943417648083019611294;
	wtt[9]  = -0.0142483904822400753741381;
	wtt[10] =  0.0416812529998231580614934;
	wtt[11] =  0.00323522898964475022578598;
	wtt[12] =  0.000813552611571786631179287;
	wtt[13] =  0.00228403610697848813660369;
	wtt[14] =  0.0338971925236628645848112;
	wtt[15] = -0.00644725595698339499416262;
 	size_t wp24[78] = { 0,1,2,3,4,5,5,4,3,2,1,0,1,6,7,8,9,10,9,8,7,6,1,2,7,
			     11,12,13,13,12,11,7,2,3,8,12,14,15,14,12,8,3,4,9,13,
			     15,15,13,9,4,5,10,13,14,13,10,5,5,9,12,12,9,5,4,8,
			     11,8,4,3,7,7,3,2,6,2,1,1,0 };
	for( size_t m=0;m<=octantOrdinates-1;++m)
	    wt[m] = wtt[wp24[m]];
	break;
    }
    default:
	Insist(false,"Unspported quadrature order selected.");
	break;
    }

    // Evaluate mu and eta for octant 1
    size_t m(0);
    for(size_t i=0;i<=levels/2-1;++i)
	for( size_t j=0;j<=(levels/2)-(i+1);++j)
	{
	    mu[m]  = att[i];
	    eta[m] = att[j];
	    ++m;
	}

    // Evaluate mu and eta for octants 2-4
    for(size_t octant=2; octant<=4; ++octant)
	for(size_t n=0; n<=octantOrdinates-1; ++n) 
	{
	    m = (octant-1)*octantOrdinates+n;
	    switch (octant) {
	    case 2:
		mu[m]  = -mu[n];
		eta[m] =  eta[n];
		wt[m]  =  wt[n];
		break;
		
	    case 3:
		mu[m]  = -mu[n];
		eta[m] = -eta[n];
		wt[m]  =  wt[n];
		break;
		
	    case 4:
		mu[m]  =  mu[n];
		eta[m] = -eta[n];
		wt[m]  =  wt[n];
		break;
	    default:
		Insist(false,"Octant value should only be 2, 3 or 4 in this loop.");
		break;
	    }
	}
    
    // Evaluate mu and eta for octants 5-8
    for( size_t n=0; n<=numOrdinates/2-1; ++n)
    {
	mu[n+numOrdinates/2]  = mu[n];
	eta[n+numOrdinates/2] = eta[n];
	wt[n+numOrdinates/2]  = wt[n];
    }
    
    // Evaluate xi for all octants
    for( size_t n=0;n<=numOrdinates/2-1;++n)
	xi[n] = std::sqrt(1.0-(mu[n]*mu[n]+eta[n]*eta[n]));
    
    for(size_t n=0;n<=numOrdinates/2-1;++n)
	xi[n+numOrdinates/2] = -xi[n];
    
    // Normalize the quadrature set
    double wsum = 0.0;
    for(size_t n=0;n<=numOrdinates-1;++n)
	wsum = wsum + wt[n];
    
    for(size_t n=0; n<=numOrdinates-1; ++n)
	wt[n] = wt[n]*(norm/wsum);
    
    // clear the temporaries.
    att.clear();
    wtt.clear();

    // Verify that the quadrature meets our integration requirements.
    Ensure( soft_equiv(iDomega(),norm) );

#ifdef ENSURE_ON
    // check each component of the vector result
    vector<double> iod = iOmegaDomega();
    Ensure( soft_equiv(iod[0],0.0) );
    Ensure( soft_equiv(iod[1],0.0) );
    Ensure( soft_equiv(iod[2],0.0) );

    // check each component of the tensor result
    vector<double> iood = iOmegaOmegaDomega();
    Ensure( soft_equiv(iood[0],norm/3.0) );  // mu*mu
    Ensure( soft_equiv(iood[1],0.0) ); // mu*eta
    Ensure( soft_equiv(iood[2],0.0) ); // mu*xi
    Ensure( soft_equiv(iood[3],0.0) ); // eta*mu
    Ensure( soft_equiv(iood[4],norm/3.0) ); // eta*eta
    Ensure( soft_equiv(iood[5],0.0) ); // eta*xi
    Ensure( soft_equiv(iood[6],0.0) ); // xi*mu
    Ensure( soft_equiv(iood[7],0.0) ); // xi*eta
    Ensure( soft_equiv(iood[8],norm/3.0) ); // xi*xi
#endif
    // ENSURE_ON

    // Copy quadrature data { mu, eta, xi } into the vector omega.
    omega.resize( numOrdinates );
    size_t ndims = dimensionality();
    for ( size_t i=0; i<numOrdinates; ++i )
    {
	omega[i].resize( ndims );
	omega[i][0] = mu[i];
	omega[i][1] = eta[i];
	omega[i][2] = xi[i];
    }

} // end of Q3LevelSym() constructor.

//---------------------------------------------------------------------------//

void Q3DLevelSym::display() const 
{
    using std::cout;
    using std::endl;
    using std::setprecision;

    cout << endl << "The Quadrature directions and weights are:" 
	 << endl << endl;
    cout << "   m  \t    mu        \t    eta       \t    xi        \t     wt      " << endl;
    cout << "  --- \t------------- \t------------- \t------------- \t-------------" << endl;
    double sum_wt = 0.0;
    for ( size_t ix = 0; ix < mu.size(); ++ix ) {
	cout << "   "
	     << ix << "\t"
	     << setprecision(10) << mu[ix]  << "\t"
	     << setprecision(10) << eta[ix] << "\t"
	     << setprecision(10) << xi[ix]  << "\t"
	     << setprecision(10) << wt[ix]  << endl;
	sum_wt += wt[ix];
    }
    cout << endl << "  The sum of the weights is " << sum_wt << endl;
    cout << endl;
}

} // end namespace rtt_quadrature

//---------------------------------------------------------------------------//
//                 end of Q3DLevelSym.cc
//---------------------------------------------------------------------------//
