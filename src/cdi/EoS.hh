//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi/EoS.hh
 * \author Kelly Thompson
 * \date   Fri Apr 13 16:15:59 2001
 * \brief  EoS class header file (an abstract class)
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __cdi_EoS_hh__
#define __cdi_EoS_hh__

#include <vector>
// #include <string>

namespace rtt_cdi
{
    
    //========================================================================
    /*!
     * \class EoS
     *
     * \brief This is a pure virtual class that defines a standard
     *  interface for all derived EoS objects. 
     *
     * Any derived EoS object must provide as a minumum the
     * functionality outlined in this routine.  This functionality
     * includes access to the data grid and the ability to return
     * interpolated opacity values.
     */
    
//     /*!
//      * \example cdi/test/tDummyOpacity.cc
//      * \example cdi/test/tCDI.cc
//      *
//      */
    //========================================================================
    
    class EoS
    {
	// DATA
	
	// There is no data for a pure virtual object.  This class
	// provides an interface and does not preserve state.
	
      public:
	
	// ---------- //
	// Destructor //
	// ---------- //
	
	/*!
	 * \brief Default EoS() destructor.
	 *
	 * This is required to correctly release memory when any
	 * object derived from EoS is destroyed.
	 */
	virtual ~EoS();
	
	// --------- //
	// Accessors //
	// --------- //
	
	/*!
	 * \brief EoS accessor that returns a single specific electron 
	 *        internal energy that corresponds to the provided
	 *        temperature and density. 
	 *
	 * \param temperature The temperature value for which an
	 *     opacity value is being requested (Kelvin).
	 * \param density The density value for which an opacity 
	 *     value is being requested (g/cm^3).
	 * \return A specific electron internal energy (kJ/g).
	 */
	virtual double getSpecificElectronInternalEnergy(
	    double density, double temperature ) const = 0;
	
	/*!
	 * \brief EoS accessor that returns a vector of specific
	 *     electron internal energies that
	 *     correspond to the provided vectors of temperatures and 
	 *     densities. 
	 *
	 * \param temperature A vector of temperature values for
	 *     which the EoS values are being requested (Kelvin).
	 * \param density A vector of density values for
	 *     which the EoS values are being requested (g/cm^3).
	 * \return A vector of specific electron internal energies (kJ/g).
	 */
	virtual std::vector< double > getSpecificElectronInternalEnergy(
	    const std::vector< double >& vdensity, 
	    const std::vector< double >& vtemperature ) const = 0;

	/*!
	 * \brief Retrieve the electron based heat capacity for this
	 *        material at the provided density and temperature.
	 *
	 * \param density Density of the material in g/cm^3
	 * \param temperature Temperature of the material in Kelvin.
	 * \return The electron based heat capacity in kJ/g/K.
	 */
	virtual double getElectronHeatCapacity(
	    double density, double temperature ) const = 0;

	/*!
	 * \brief Retrieve a set of electron based heat capacities for
	 *        this material that correspond to the tuple list of
	 *        provided densities and temperatures. 
	 *
	 * \param density Density of the material in g/cm^3
	 * \param temperature Temperature of the material in Kelvin.
	 * \return The electron based heat capacity in kJ/g/K.
	 */
	virtual std::vector< double > getElectronHeatCapacity(
	    const std::vector< double >& density, 
	    const std::vector< double >& temperature ) const = 0;

	/*!
	 * \brief Retrieve the specific ion internal energy for this
	 *        material at the provided density and temperature.
	 *
	 * \param density Density of the material in g/cm^3
	 * \param temperature Temperature of the material in Kelvin.
	 * \return The specific ion internal energy in kJ/g.
	 */
	virtual double getSpecificIonInternalEnergy(
	    double density, double temperature ) const = 0; 

	/*!
	 * \brief Retrieve a set of specific ion internal energies for
	 *        this material that correspond to the tuple list of
	 *        provided densities and temperatures. 
	 *
	 * \param density Density of the material in g/cm^3
	 * \param temperature Temperature of the material in Kelvin.
	 * \return A vector of specific ion internal energies in kJ/g.
	 */
	virtual std::vector< double > getSpecificIonInternalEnergy(
	    const std::vector< double >& density, 
	    const std::vector< double >& temperature ) const = 0;

	/*!
	 * \brief Retrieve the ion based heat capacity for this
	 *        material at the provided density and temperature.
	 *
	 * \param density Density of the material in g/cm^3
	 * \param temperature Temperature of the material in Kelvin.
	 * \return The ion based heat capacity in kJ/g/K.
	 */
	virtual double getIonHeatCapacity(
	    double density, double temperature ) const = 0;

	/*!
	 * \brief Retrieve a set of ion based heat capacities for
	 *        this material that correspond to the tuple list of
	 *        provided densities and temperatures. 
	 *
	 * \param density Density of the material in g/cm^3
	 * \param temperature Temperature of the material in Kelvin.
	 * \return A vector of ion based heat capacities in kJ/g/K.
	 */
	virtual std::vector< double > getIonHeatCapacity(
	    const std::vector< double >& density,
	    const std::vector< double >& temperature ) const = 0;

	/*!
	 * \brief Retrieve the number of free electrons per ion for this
	 *        material at the provided density and temperature.
	 *
	 * \param density Density of the material in g/cm^3
	 * \param temperature Temperature of the material in Kelvin.
	 * \return The number of free electrons per ion.
	 */
	virtual double getNumFreeElectronsPerIon(
	    double density, double temperature ) const = 0;

	/*!
	 * \brief Retrieve a set of free electrons per ion averages for
	 *        this material that correspond to the tuple list of
	 *        provided densities and temperatures. 
	 *
	 * \param density Density of the material in g/cm^3
	 * \param temperature Temperature of the material in Kelvin.
	 * \return A vector of the number of free electrons per ion.
	 */
	virtual std::vector< double > getNumFreeElectronsPerIon(
	    const std::vector< double >& density,
	    const std::vector< double >& temperature ) const = 0;

	/*!
	 * \brief Retrieve the electron based thermal conductivity for this
	 *        material at the provided density and temperature.
	 *
	 * \param density Density of the material in g/cm^3
	 * \param temperature Temperature of the material in Kelvin.
	 * \return The electron based thermal conductivity in 1/s/cm.
	 */
	virtual double getElectronBasedThermalConductivity(
	    double density, double temperature ) const = 0; 

	/*!
	 * \brief Retrieve a set of electron based thermal conductivities for
	 *        this material that correspond to the tuple list of
	 *        provided densities and temperatures. 
	 *
	 * \param density Density of the material in g/cm^3
	 * \param temperature Temperature of the material in Kelvin.
	 * \return A vector of electron based thermal conductivities
	 * in 1/s/cm.
	 */
	virtual std::vector< double > getElectronBasedThermalConductivity(
	    const std::vector< double >& density,
	    const std::vector< double >& temperature ) const = 0;
	
    }; // end of class EoS    
    
} // end namespace rtt_cdi

#endif // __cdi_EoS_hh__

//---------------------------------------------------------------------------//
// end of cdi/EoS.hh
//---------------------------------------------------------------------------//
