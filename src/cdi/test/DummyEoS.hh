//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi/test/DummyEoS.hh
 * \author Kelly Thompson
 * \date   Mon Jan 8 15:29:17 2001
 * \brief  DummyEoS class header file (derived from ../EoS)
 */
//---------------------------------------------------------------------------//
// $Id$1
//---------------------------------------------------------------------------//

#ifndef __cdi_DummyEoS_hh__
#define __cdi_DummyEoS_hh__

#include "../EoS.hh"

namespace rtt_dummyEoS
{

    //========================================================================
    /*!
     * \class DummyEoS
     *
     * \brief This is an equation of state class that derives its
     *        interface from cdi/EoS and is used for testing purposes
     *        only. 
     *
     * This equation of state class always contains the same data (set by the
     * default constructor).  The data table has the following properties:
     *
     *     internal energy = temperature + 1000*density
     *     heat capacity   = temperature + density/1000
     *     num free electrons = temperature/100
     *     thermal conductivity = 1000*temperature + density
     *
//      * In addition to providing definitions for the member functions
//      * outlined in EoS this class provides three additional 1-D
//      * STL-like accessors for opacity data.
     */
    
    /*!
     * \example cdi/test/tEoS.cc
     * \example cdi/test/tCDI.cc
     */
    //========================================================================

    class DummyEoS : public rtt_cdi::EoS
    {
      public:
	
	// -------------------------- //
	// Constructors & Destructors //
	// -------------------------- //
	
	/*!
	 * \brief Constructor for DummyEoS object.
	 * 
	 * The constructor assigns fixed values for all of the member
	 * data.  Every instance of this object has the same member
	 * data. 
	 */
	DummyEoS();

	/*!
	 * \brief Default DummyEoS() destructor.
	 *
	 * This is required to correctly release memory when a
	 * DummyEoS object is destroyed.
	 */
	~DummyEoS();
	
	// --------- //
	// Accessors //
	// --------- //
	
	/*!
	 * \brief EoS accessor that returns a single specific electron 
	 *        internal energy that corresponds to the provided
	 *        temperature and density. 
	 *
	 *    internal energy = temperature + 1000*density
	 *
	 * \param temperature The temperature value for which an
	 *     opacity value is being requested (Kelvin).
	 * \param density The density value for which an opacity 
	 *     value is being requested (g/cm^3).
	 * \return A specific electron internal energy (kJ/g).
	 */
	double getSpecificElectronInternalEnergy(
	    double temperature, double density ) const;
	
	/*!
	 * \brief EoS accessor that returns a vector of specific
	 *     electron internal energies that
	 *     correspond to the provided vectors of temperatures and 
	 *     densities. 
	 *
	 *    internal energy[i] = temperature[i] + 1000*density[i]
	 *
	 * \param temperature A vector of temperature values for
	 *     which the EoS values are being requested (Kelvin).
	 * \param density A vector of density values for
	 *     which the EoS values are being requested (g/cm^3).
	 * \return A vector of specific electron internal energies (kJ/g).
	 */
	std::vector< double > getSpecificElectronInternalEnergy(
	    const std::vector< double >& vtemperature, 
	    const std::vector< double >& vdensity ) const;

	/*!
	 * \brief Retrieve the electron based heat capacity for this
	 *        material at the provided density and temperature.
	 *
	 *     heat capacity = temperature + density/1000
	 *
	 * \param density Density of the material in g/cm^3
	 * \param temperature Temperature of the material in Kelvin.
	 * \return The electron based heat capacity in kJ/g/K.
	 */
	double getElectronHeatCapacity(
	    double temperature, double density ) const;

	/*!
	 * \brief Retrieve a set of electron based heat capacities for
	 *        this material that correspond to the tuple list of
	 *        provided densities and temperatures. 
	 *
	 *     heat capacity = temperature + density/1000
	 *
	 * \param density Density of the material in g/cm^3
	 * \param temperature Temperature of the material in Kelvin.
	 * \return The electron based heat capacity in kJ/g/K.
	 */
	std::vector< double > getElectronHeatCapacity(
	    const std::vector< double >& vtemperature, 
	    const std::vector< double >& vdensity ) const;

	/*!
	 * \brief Retrieve the specific ion internal energy for this
	 *        material at the provided density and temperature.    
	 *
	 *     internal energy = temperature + 1000*density
	 *
	 * \param density Density of the material in g/cm^3
	 * \param temperature Temperature of the material in Kelvin.
	 * \return The specific ion internal energy in kJ/g.
	 */
	double getSpecificIonInternalEnergy(
	    double temperature, double density ) const;

	/*!
	 * \brief Retrieve a set of specific ion internal energies for
	 *        this material that correspond to the tuple list of
	 *        provided densities and temperatures.      
	 *
	 *     internal energy = temperature + 1000*density
	 *
	 * \param density Density of the material in g/cm^3
	 * \param temperature Temperature of the material in Kelvin.
	 * \return A vector of specific ion internal energies in kJ/g.
	 */
	std::vector< double > getSpecificIonInternalEnergy(
	    const std::vector< double >& vtemperature, 
	    const std::vector< double >& vdensity ) const;

	/*!
	 * \brief Retrieve the ion based heat capacity for this
	 *        material at the provided density and temperature.
	 *
	 *     heat capacity   = temperature + density/1000
	 *
	 * \param density Density of the material in g/cm^3
	 * \param temperature Temperature of the material in Kelvin.
	 * \return The ion based heat capacity in kJ/g/K.
	 */
	double getIonHeatCapacity(
	    double temperature, double density ) const;

	/*!
	 * \brief Retrieve a set of ion based heat capacities for
	 *        this material that correspond to the tuple list of
	 *        provided densities and temperatures. 
	 *
	 *     heat capacity   = temperature + density/1000
	 *
	 * \param density Density of the material in g/cm^3
	 * \param temperature Temperature of the material in Kelvin.
	 * \return A vector of ion based heat capacities in kJ/g/K.
	 */
	std::vector< double > getIonHeatCapacity(
	    const std::vector< double >& vtemperature,
	    const std::vector< double >& vdensity ) const;

	/*!
	 * \brief Retrieve the number of free electrons per ion for this
	 *        material at the provided density and temperature.
	 *
	 *     num free electrons = temperature/100
	 *
	 * \param density Density of the material in g/cm^3
	 * \param temperature Temperature of the material in Kelvin.
	 * \return The number of free electrons per ion.
	 */
	double getNumFreeElectronsPerIon(
	    double temperature, double density ) const;

	/*!
	 * \brief Retrieve a set of free electrons per ion averages for
	 *        this material that correspond to the tuple list of
	 *        provided densities and temperatures. 
	 *
	 * \param density Density of the material in g/cm^3
	 * \param temperature Temperature of the material in Kelvin.
	 * \return A vector of the number of free electrons per ion.
	 */
	std::vector< double > getNumFreeElectronsPerIon(
	    const std::vector< double >& vtemperature,
	    const std::vector< double >& vdensity ) const;

	/*!
	 * \brief Retrieve the electron based thermal conductivity for this
	 *        material at the provided density and temperature.
	 *
	 *     thermal conductivity = 1000*temperature + density
	 *
	 * \param density Density of the material in g/cm^3
	 * \param temperature Temperature of the material in Kelvin.
	 * \return The electron based thermal conductivity in 1/s/cm.
	 */
	double getElectronBasedThermalConductivity(
	    double temperature, double density ) const;

	/*!
	 * \brief Retrieve a set of electron based thermal conductivities for
	 *        this material that correspond to the tuple list of
	 *        provided densities and temperatures. 
	 *
	 *     thermal conductivity = 1000*temperature + density
	 *
	 * \param density Density of the material in g/cm^3
	 * \param temperature Temperature of the material in Kelvin.
	 * \return A vector of electron based thermal conductivities
	 * in 1/s/cm.
	 */
	std::vector< double > getElectronBasedThermalConductivity(
	    const std::vector< double >& vtemperature,
	    const std::vector< double >& vdensity ) const;

    }; // end of class DummyEoS
    
} // end namespace rtt_dummyEoS

#endif // __cdi_DummyEoS_hh__

//---------------------------------------------------------------------------//
// end of cdi/test/DummyEoS.hh
//---------------------------------------------------------------------------//
