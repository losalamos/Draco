//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi_analytic/test/tstAnalytic_EoS.cc
 * \author Thomas M. Evans
 * \date   Thu Oct  4 11:45:19 2001
 * \brief  Analytic_EoS test.
 * \note   Copyright (C) 2001-2012 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "../Analytic_EoS.hh"
#include "../Analytic_Models.hh"
#include "cdi/CDI.hh"
#include "cdi/EoS.hh"
#include "ds++/ScalarUnitTest.hh"
#include "ds++/Assert.hh"
#include "ds++/SP.hh"
#include "ds++/Soft_Equivalence.hh"
#include "ds++/Release.hh"

#include <iostream>
#include <vector>
#include <cmath>
#include <typeinfo>

using namespace std;

using rtt_cdi_analytic::Analytic_EoS;
using rtt_cdi_analytic::Analytic_EoS_Model;
using rtt_cdi_analytic::Polynomial_Specific_Heat_Analytic_EoS_Model;
using rtt_cdi::CDI;
using rtt_cdi::EoS;
using rtt_dsxx::SP;
using rtt_dsxx::soft_equiv;

#define PASSMSG(m) ut.passes(m)
#define FAILMSG(m) ut.failure(m)
#define ITFAILS    ut.failure( __LINE__, __FILE__ )

//---------------------------------------------------------------------------//
// TESTS
//---------------------------------------------------------------------------//

void analytic_eos_test( rtt_dsxx::UnitTest & ut )
{
    typedef Polynomial_Specific_Heat_Analytic_EoS_Model Polynomial_Model;

    // make an analytic model (polynomial specific heats)
    // elec specific heat = a + bT^c
    // ion specific heat  = d + eT^f
    SP<Polynomial_Model> model(new Polynomial_Model(0.0,1.0,3.0,
						    0.2,0.0,0.0));

    if (!model) FAILMSG("Failed to build an Analytic EoS Model!");

    // make an analtyic eos
    Analytic_EoS analytic(model);

    // checks
    {
	double T   = 5.0;
	double rho = 3.0;

	double Cve = T*T*T;
	double Cvi = 0.2;

	double Ue = T*T*T*T/4.0;
	double Ui = 0.2*T;

	// specific heats
	if (analytic.getElectronHeatCapacity(T,rho) != Cve)           ITFAILS;
	if (analytic.getIonHeatCapacity(T,rho) != Cvi)                ITFAILS;

	// specific internal energies
	if (!soft_equiv(analytic.getSpecificElectronInternalEnergy(T,rho), Ue))
	    ITFAILS;

	if (!soft_equiv(analytic.getSpecificIonInternalEnergy(T,rho), Ui))
	    ITFAILS;

	// everything else is zero
	if (analytic.getNumFreeElectronsPerIon(T,rho) != 0.0)         ITFAILS;
	if (analytic.getElectronThermalConductivity(T,rho) != 0.0)    ITFAILS;
    }

    // Check the root finder for new Te, given delta Ue.
    {
        double rho = 3.0;  // not currently used by getElectronTemperature().
        double Ue = 4.0;
        double Tguess = 2.0;
        
        double T_new = analytic.getElectronTemperature( rho, Ue, Tguess );
        if( ! soft_equiv( T_new, 2.0 ) ) ITFAILS;
    }

    // field check
    {
	vector<double> T(6);
	vector<double> rho(6);
	
	T[0] = .993;
	T[1] = .882;
	T[2] = .590;
	T[3] = .112;
	T[4] = .051;
	T[5] = .001;

	std::fill(rho.begin(), rho.end(), 3.0);
	rho[3] = 2.5;

	vector<double> Cve = analytic.getElectronHeatCapacity(T,rho);
	vector<double> Cvi = analytic.getIonHeatCapacity(T,rho);
	vector<double> eie = analytic.getSpecificElectronInternalEnergy(T,rho);
	vector<double> iie = analytic.getSpecificIonInternalEnergy(T,rho);
	vector<double> nfe = analytic.getNumFreeElectronsPerIon(T,rho);
	vector<double> etc = analytic.getElectronThermalConductivity(T,rho);

	if (Cve.size() != 6) ITFAILS;
	if (Cvi.size() != 6) ITFAILS;
	if (eie.size() != 6) ITFAILS;
	if (iie.size() != 6) ITFAILS;
	if (nfe.size() != 6) ITFAILS;
	if (etc.size() != 6) ITFAILS;

	for (int i = 0; i < 6; i++)
	{
	    double cve_ref = T[i] * T[i] * T[i];

	    double cvi_ref = 0.2;

	    double Ue_ref  = T[i] * T[i] * T[i] * T[i] / 4.0;

	    double Ui_ref  = 0.2 * T[i];

	    if (!soft_equiv(Cve[i], cve_ref)) ITFAILS;
	    if (!soft_equiv(Cvi[i], cvi_ref)) ITFAILS;
	    if (!soft_equiv(eie[i], Ue_ref)) ITFAILS;
	    if (!soft_equiv(iie[i], Ui_ref)) ITFAILS;
	    
	    // all else are zero
	    if (nfe[i] != 0.0) ITFAILS;
	    if (etc[i] != 0.0) ITFAILS;
	}
    }

    // Test the get_Analytic_Model() member function.
    {
        SP<Polynomial_Model const> myEoS_model
            = analytic.get_Analytic_Model();
        SP<Polynomial_Model const> expected_model( model );
        
        if( expected_model == myEoS_model )
            PASSMSG("get_Analytic_Model() returned the expected EoS model.");
        else
            FAILMSG("get_Analytic_Model() did not return the expected EoS model.");
    }

    // Test the get_parameters() members function
    {
        std::vector<double> params( model->get_parameters() );
        
        std::vector<double> expectedValue(6,0.0);
        expectedValue[1] = 1.0;
        expectedValue[2] = 3.0;
        expectedValue[3] = 0.2;

        double const tol(1.0e-12);
        
        if( params.size() != expectedValue.size() ) ITFAILS;
        
        if( soft_equiv( params.begin(), params.end(),
                        expectedValue.begin(), expectedValue.end(), tol ) )
            PASSMSG("get_parameters() returned the analytic expression coefficients.");
        else
            FAILMSG("get_parameters() did not return the analytic expression coefficients.");
    }
   
    // make an analytic model like Su-Olson (polynomial specific heats)
    // elec specific heat = a + bT^c
    // ion specific heat  = d + eT^f
    SP<Polynomial_Model> so_model(new Polynomial_Model(0.0,54880.0,3.0,
                                  0.2,0.0,0.0));
    // make an analtyic eos
    Analytic_EoS so_analytic(so_model);

    // do a sanity check on the table lookups.
    {
        double rho = 1.0;  // not currently used by getElectronTemperature().
        double Ue =  1.371999999e-20;
        double Te = 1.0e-6;
       
        // First request the internal energy given the temperature 
	double Ue0  = so_analytic.getSpecificElectronInternalEnergy(Te,rho);
        if( ! soft_equiv( Ue, Ue0 ) ) ITFAILS;
        // Next request the temperature given the internal energy (should match).
        double T_new = so_analytic.getElectronTemperature( rho, Ue0, Te );
        if( ! soft_equiv( T_new, Te ) ) ITFAILS;
    }
 
    return;
}

//---------------------------------------------------------------------------//
 
void CDI_test(rtt_dsxx::UnitTest & ut)
{
    typedef Polynomial_Specific_Heat_Analytic_EoS_Model Polynomial_Model;

    // cdi object
    CDI eosdata;

    // analytic model
    SP<Analytic_EoS_Model> model(new Polynomial_Model(0.0,1.0,3.0,
						      0.0,0.0,0.0));

    // assign the eos object
    SP<Analytic_EoS> analytic_eos(new Analytic_EoS(model));

    // EoS object
    SP<const EoS> eos = analytic_eos;
    if (typeid(*eos) != typeid(Analytic_EoS)) ITFAILS;

    // Assign the object to cdi
    eosdata.setEoS(eos);

    // check
    if (!eosdata.eos()) FAILMSG("Can't reference EoS smart pointer");

    // make temperature and density fields
    vector<double> T(6);
    vector<double> rho(6);
    
    T[0] = .993;
    T[1] = .882;
    T[2] = .590;
    T[3] = .112;
    T[4] = .051;
    T[5] = .001;
    
    std::fill(rho.begin(), rho.end(), 3.0);
    rho[3] = 2.5;
    
    vector<double> Cve;
    vector<double> Cvi;
    vector<double> eie;
    vector<double> iie;
    vector<double> nfe;
    vector<double> etc;

    // test the data
    {
    
	Cve = eosdata.eos()->getElectronHeatCapacity(T,rho);
	Cvi = eosdata.eos()->getIonHeatCapacity(T,rho);
	eie = eosdata.eos()->getSpecificElectronInternalEnergy(T,rho);
	iie = eosdata.eos()->getSpecificIonInternalEnergy(T,rho);
	nfe = eosdata.eos()->getNumFreeElectronsPerIon(T,rho);
	etc = eosdata.eos()->getElectronThermalConductivity(T,rho);

	if (Cve.size() != 6) ITFAILS;
	if (Cvi.size() != 6) ITFAILS;
	if (eie.size() != 6) ITFAILS;
	if (iie.size() != 6) ITFAILS;
	if (nfe.size() != 6) ITFAILS;
	if (etc.size() != 6) ITFAILS;

	for (int i = 0; i < 6; i++)
	{
	    double cve_ref = T[i] * T[i] * T[i];
	    double ue_ref = T[i] * T[i] * T[i] * T[i] / 4.0;

	    if (!soft_equiv(Cve[i], cve_ref)) ITFAILS;
	    if (!soft_equiv(eie[i], ue_ref)) ITFAILS;
	    
	    // all else are zero
	    if (Cvi[i] != 0.0) ITFAILS;
	    if (iie[i] != 0.0) ITFAILS;
	    if (nfe[i] != 0.0) ITFAILS;
	    if (etc[i] != 0.0) ITFAILS;
	}
    }

    // now reset the CDI
    eosdata.reset();

    // should catch this
    bool caught = false;
    try
    {
	eosdata.eos();
    }
    catch(const rtt_dsxx::assertion &ass)
    {
	PASSMSG("Good, caught an unreferenced EoS SP!");
	caught = true;
    }
    if (!caught) FAILMSG("Failed to catch an unreferenced SP<EoS>!");

    // now assign the analytic eos to CDI directly
    eosdata.setEoS(analytic_eos);
    if (!eosdata.eos())                                 ITFAILS;
    if (typeid(*eosdata.eos()) != typeid(Analytic_EoS)) ITFAILS;

    // now test the data again

    // test the data
    {
	Cve = eosdata.eos()->getElectronHeatCapacity(T,rho);
	Cvi = eosdata.eos()->getIonHeatCapacity(T,rho);
	eie = eosdata.eos()->getSpecificElectronInternalEnergy(T,rho);
	iie = eosdata.eos()->getSpecificIonInternalEnergy(T,rho);
	nfe = eosdata.eos()->getNumFreeElectronsPerIon(T,rho);
	etc = eosdata.eos()->getElectronThermalConductivity(T,rho);

	if (Cve.size() != 6) ITFAILS;
	if (Cvi.size() != 6) ITFAILS;
	if (eie.size() != 6) ITFAILS;
	if (iie.size() != 6) ITFAILS;
	if (nfe.size() != 6) ITFAILS;
	if (etc.size() != 6) ITFAILS;

	for (int i = 0; i < 6; i++)
	{
	    double cve_ref = T[i] * T[i] * T[i];
	    double eie_ref = 0.25 * T[i] * T[i] * T[i] * T[i];

	    if (!soft_equiv(Cve[i], cve_ref)) ITFAILS;
	    if (!soft_equiv(eie[i], eie_ref)) ITFAILS;
	    
	    // all else are zero
	    if (Cvi[i] != 0.0) ITFAILS;
	    if (iie[i] != 0.0) ITFAILS;
	    if (nfe[i] != 0.0) ITFAILS;
	    if (etc[i] != 0.0) ITFAILS;
	}
    }
    
    return;
}

//---------------------------------------------------------------------------//

void packing_test(rtt_dsxx::UnitTest & ut)
{
    typedef Polynomial_Specific_Heat_Analytic_EoS_Model Polynomial_Model;

    vector<char> packed;

    {
	// make an analytic model (polynomial specific heats)
	SP<Polynomial_Model> model(new Polynomial_Model(0.0,1.0,3.0,
							0.2,0.0,0.0));
	
	// make an analtyic eos
	SP<EoS> eos(new Analytic_EoS(model));

	packed = eos->pack();
    }
    
    Analytic_EoS neos(packed);

    // checks
    {
	double T   = 5.0;
	double rho = 3.0;

	double Cve = T*T*T;
	double Cvi = 0.2;

	double Ue  = 0.25*T*T*T*T;
	double Ui  = 0.2 * T;

	// specific heats
	if (neos.getElectronHeatCapacity(T,rho) != Cve)           ITFAILS;
	if (neos.getIonHeatCapacity(T,rho) != Cvi)                ITFAILS;

	// specific internal energies
	if (!soft_equiv(neos.getSpecificElectronInternalEnergy(T,rho), Ue))
	    ITFAILS;

	if (!soft_equiv(neos.getSpecificIonInternalEnergy(T,rho), Ui))
	    ITFAILS;


	// everything else is zero
	if (neos.getNumFreeElectronsPerIon(T,rho) != 0.0)         ITFAILS;
	if (neos.getElectronThermalConductivity(T,rho) != 0.0)    ITFAILS;
    }

    if( ut.numFails == 0 )
	PASSMSG("Analytic EoS packing/unpacking test successfull.");

    return;
}

//---------------------------------------------------------------------------//

int main(int argc, char *argv[])
{
    rtt_dsxx::ScalarUnitTest ut( argc, argv, rtt_dsxx::release );
    try
    {
	// >>> UNIT TESTS
	analytic_eos_test(ut);
	CDI_test(ut);
	packing_test(ut);
    }
    catch (rtt_dsxx::assertion &err)
    {
        std::cout << "ERROR: While testing " << argv[0] << ", "
                  << err.what() << std::endl;
        ut.numFails++;
    }
    catch( ... )
    {
        std::cout << "ERROR: While testing " << argv[0] << ", " 
                  << "An unknown exception was thrown on processor "
                  << std::endl;
        ut.numFails++;
    }
    return ut.numFails;
}   

//---------------------------------------------------------------------------//
// end of tstAnalytic_EoS.cc
//---------------------------------------------------------------------------//
