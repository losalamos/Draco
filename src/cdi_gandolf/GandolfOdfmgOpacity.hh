//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi_gandolf/GandolfOdfmgOpacity.hh
 * \author Kelly Thompson
 * \date   Mon Jan 22 13:56:01 2001
 * \brief  GandolfOdfmgOpacity class header file (derived from
 *         cdi/OdfmgOpacity) 
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __cdi_gandolf_GandolfOdfmgOpacity_hh__
#define __cdi_gandolf_GandolfOdfmgOpacity_hh__

#include <vector>
#include <string>
#include <cmath> // we need to define log(double) and exp(double)

#include "ds++/Assert.hh" 
#include "ds++/SP.hh"
#include "cdi/OdfmgOpacity.hh"
#include "cdi/OpacityCommon.hh"

#include "GandolfWrapper.hh"    // we make calls to the wrapper routines.
#include "GandolfDataTable.hh"  // we have a smart pointer to a
                                // GandolfDataTable object.

namespace rtt_cdi_gandolf
{
// -------------------- //
// Forward declarations //
// -------------------- //

class GandolfFile;
class GandolfDataTable;

//===========================================================================//
/*!
 * \class GandolfOdfmgOpacity
 *
 * \brief GandolfOdfmgOpacity allows the client code to retrieve
 *        opacity data for a particular material.  Each GandolfOpacity
 *        object represents a specific type of data defined by five
 *        attributes: an IPCRESS File (via a GandolfFile object), a
 *        material identifier, an energy model (already selecte since
 *        this is a Odfmg Opacity class), a physics model and a
 *        reaction type.
 *
 *      This is a concrete class derived from cdi/OdfmgOpacity.
 *      This class allows to client to access the data in IPCRESS
 *      files via the Gandolf libraries.
 * <p>
 *      This class is designed to be used in conjuction with the CDI.
 *      The client code will create a GandolfOdfmgOpacity object
 *      and use this object as an argument during the CDI
 *      instantiation.  The purpose of this class is to provide a
 *      mechanism for accessing data in IPCRESS files and works by
 *      calling the Gandolf library provided by X-5.  The
 *      GandolfOdfmgOpacity constructor expects four arguments: a
 *      hook to IPCRESS data file (spGandolfFile), a material
 *      identifier, an opacity model (Rosseland or Plank) and an
 *      opacity reaction specifier (total, scattering or absorption).
 *      Once constructed, this object allows the client to access any
 *      data found in the IPCRESS file for that one material.  The
 *      client code will need to create a separate
 *      GandolfOdfmgOpacity object for each material that it
 *      needs information about. Multiple opacity objects can exist
 *      per IPCRESS file.
 * <p> 
 *      This class only provides access to multigroup opacity data.
 *      If the user needs gray opacity IPCRESS data he/she should use
 *      the cdi_gandolf/GandolfGrayOpacity class.
 * <p>
 *      When instantiated, the GandolfOdfmgOpacity object creates
 *      a GandolfDataTable object.  The IPCRESS data is cached in this
 *      table object.  When the client requests an opacity value at a
 *      specified temperature and density the GandolfOdfmgOpcity
 *      object calls the appropriate GANDOLF library routine, which in
 *      turn, interpolates on the data cached in the GandolfDataTable
 *      object.
 * <p>
 *      When compiling DRACO with support for the IPCRESS file reader
 *      (via Gandolf) you must add the following option on the
 *      configure line:<br><br>
 *      <tt>   --with-gandolf-lib=${VENDORS}/gandolf/IRIX64/lib64</tt><br><br>
 *      Currently, the Gandolf library is only available on 64 bit
 *      IRIX architectures.
 */

/*!
 * \example cdi_gandolf/test/tGandolfOpacity.cc
 *
 * Example of GandolfOdfmgOpacity usage independent of CDI.  In
 *      this example we construct a GandolfOdfmgOpacity object for
 *      the material Aluminum (matID=10001 in our example IPCRESS
 *      file).  We then use the GandolfOdfmgOpacity object to compute
 *      Rosseland Odfmg opacity values for a specified material,
 *      temperature and density.  Other forms of the getOpacity()
 *      accessor are tested along with accessors that return
 *      information about the data set and the cached data table.
 *
 * \example cdi_gandolf/test/tGandolfWithCDI.cc
 * 
 * This example tests and demonstrates how to use the cdi_gandolf
 * package as a plug-in for the CDI class.
 */
//===========================================================================//

class GandolfOdfmgOpacity : public rtt_cdi::OdfmgOpacity
{

	// DATA

	// ----------------------- //
	// Specify unique material //
	// ----------------------- //

	/*!
	 * \brief DS++ Smart Pointer to a GandolfFile object.
	 *     spGandolfFile acts as a hook to link this object to an
	 *     IPCRESS file.
	 */
	rtt_dsxx::SP< const GandolfFile > spGandolfFile;

	/*!
	 * \brief Identification number for one of the materials found in
	 *     the IPCRESS file pointed to by spGandolfFile.
	 */
	int materialID;

	// -------------------- //
	// Available data types //
	// -------------------- //

	// The IPCRESS file only holds specific data for each of its materials.

	/*!
	 * \brief Number of types of data found in the IPCRESS file.
	 */
	int numKeys;

	/*!
	 * \brief A list of keys known by the IPCRESS file.
	 */
	std::vector< std::string > vKnownKeys;

	// --------------- //
	// Data specifiers //
	// --------------- //

	/*!
	 * \brief The physics model that the current data set is based on.
	 *        { Rosseland, Plank }.  This enumeration is defined
	 *        in cdi/OpacityCommon.hh.
	 */
	rtt_cdi::Model opacityModel;

	/*!
	 * \brief The type of reaction rates that the current data set
	 *        represents { Total, Scattering, Absorption }. This
	 *        enumeration is defined in cdi/OpacityCommon.hh.
	 */
	rtt_cdi::Reaction opacityReaction;

	/*!
	 * \brief A string that identifies the energy policy for this
	 *     class.
	 */
	const std::string energyPolicyDescriptor;

	// -------------------- //
	// Opacity lookup table //
	// -------------------- //

	/*!
	 * \brief spGandolfDataTable contains a cached copy of the
	 *        requested IPCRESS opacity lookup table.
	 *
	 * There is a one-to-one relationship between GandolfOpacity and
	 * GandolfDataTable. 
	 */
	rtt_dsxx::SP< const GandolfDataTable > spGandolfDataTable;

	/*!
	 * \brief The group boundaries that we use are not the same as those read
	 * in by Gandolf, because Tops pretends that the ODFs are extra groups. We
	 * therefore have to store our own group and band structure.
	 */
	std::vector<double> groupBoundaries;

	/*!
	 * \brief The band boundaries are derived from the group boundaries read
	 * from the IPCRESS file and the number of bands input by the user.
	 * It should be a vector with size(the number of bands + 1), 0 on the left,
	 * and 1 on the right.
	 */
	std::vector<double> bandBoundaries;

	/*!
	 * \brief Sometimes the opacities in the IPCRESS file are in ascending order
	 * inside each band, sometimes not. Figure it out at the beginning and store
	 * the result so that we know whether to reverse it later.
	 */
	bool reverseBands;

	/*!
	 * \brief This private function translates the group structure stored in the
	 * IPCRESS file to the correct group and band structure, given the
	 * user-input number of bands. This should only need to be called once
	 * during initialization.
	 */
	void loadGroupsAndBands(const int numBands);
public:

	// ------------ //
	// Constructors //
	// ------------ //

	/*!
	 * \brief This is the default GandolfOdfmgOpacity
	 *     constructor.  It requires four arguments plus the energy
	 *     model (this class) to be instantiated.
	 * 
	 *     The combiniation of a data file and a material ID uniquely 
	 *     specifies a material.  If we add the Model, Reaction and
	 *     EnergyPolicy the opacity table is uniquely defined.
	 *
	 * \param spGandolfFile This smart pointer links an IPCRESS
	 *     file (via the GandolfFile object) to a GandolfOpacity
	 *     object. There may be many GandolfOpacity objects per
	 *     GandolfFile object but only one GandolfFile object for each 
	 *     GandolfOpacity object.
	 * \param materialID An identifier that links the
	 *     GandolfOpacity object to a single material found in the
	 *     specified IPCRESS file.
	 * \param opacityModel The physics model that the current
	 *     data set is based on.
	 * \param opacityReaction The type of reaction rate that the
	 *     current data set represents. 
	 */
	GandolfOdfmgOpacity( const rtt_dsxx::SP< const GandolfFile >& spGandolfFile,
			int materialID, 
			rtt_cdi::Model opacityModel,
			rtt_cdi::Reaction opacityReaction,
		  	int numBands);

	/*!
	 * \brief Unpacking constructor.
	 *
	 * This constructor unpacks a GandolfOdfmgOpacity object from a
	 * state attained through the pack function.
	 *
	 * \param packed vector<char> of packed GandolfOdfmgOpacity state;
	 * the packed state is attained by calling pack()
	 */
	explicit GandolfOdfmgOpacity(const std::vector<char> &packed);

	/*!
	 * \brief Default GandolfOpacity() destructor.
	 *
	 *     This is required to correctly release memory when a 
	 *     GandolfOdfmgOpacity is destroyed.  This constructor's
	 *     * definition must be declared in the implementation file so
	 *     that * we can avoid including too many header files
	 */
	~GandolfOdfmgOpacity();

	// --------- //
	// Accessors //
	// --------- //

	/*!
	 * \brief Opacity accessor that utilizes STL-like iterators.  This 
	 *     accessor expects a list of (temperature,density) tuples.
	 *     A set of opacity multigroup values will be returned for
	 *     each tuple.  The temperature and density iterators are
	 *     required to be the same length.  The opacity container
	 *     should have a length equal to the number of tuples times
	 *     the number of energy groups for multigroup data set.
	 * 
	 * \param temperatureFirst The beginning position of a STL
	 *     container that holds a list of temperatures.
	 * \param temperatureLast The end position of a STL
	 *     container that holds a list of temperatures.
	 * \param densityFirst The beginning position of a STL
	 *     container that holds a list of densities.
	 * \param densityLast The end position of a STL
	 *     container that holds a list of temperatures.
	 * \param opacityFirst The beginning position of a STL
	 *     container into which opacity values corresponding to the
	 *     given (temperature,density) tuple will be stored.
	 * \return A list (of type OpacityIterator) of opacities are
	 *     returned.  These opacities correspond to the temperature
	 *     and density values provied in the two InputIterators.  
	 */
	template < class TemperatureIterator, class DensityIterator,
				class OpacityIterator >
					OpacityIterator getOpacity( TemperatureIterator temperatureFirst,
							TemperatureIterator temperatureLast,
							DensityIterator densityFirst, 
							DensityIterator densityLast,
							OpacityIterator opacityFirst ) const;

	/*!
	 * \brief Opacity accessor that utilizes STL-like iterators.  This 
	 *     accessor expects a list of temperatures in an STL
	 *     container.  A set of multigroup opacity values will be
	 *     returned for each temperature provided.  The opacity
	 *     container should have a length equal to the number of
	 *     temperatures times the number of energy groups for
	 *     multigroup data set.
	 *
	 * \param temperatureFirst The beginning position of a STL
	 *     container that holds a list of temperatures.
	 * \param temperatureLast The end position of a STL
	 *     container that holds a list of temperatures.
	 * \param targetDensity The single density value used when
	 *     computing opacities for each given temperature.
	 * \param opacityFirst The beginning position of a STL
	 *     container into which opacity values corresponding to the
	 *     given temperature values will be stored.
	 * \return A list (of type OpacityIterator) of opacities are
	 *     returned.  These opacities correspond to the temperature
	 *     provided in the STL container and the single density value.
	 */
	template < class TemperatureIterator, class OpacityIterator >
		OpacityIterator getOpacity( TemperatureIterator temperatureFirst,
				TemperatureIterator temperatureLast,
				double targetDensity,
				OpacityIterator opacityFirst ) const;

	/*!
	 * \brief Opacity accessor that utilizes STL-like iterators.  This 
	 *     accessor expects a list of densities in an STL
	 *     container.  A set of multigroup opacity values will be
	 *     returned for each density provided.  The opacity
	 *     container should have a length equal to the number of
	 *     density times the number of energy groups for
	 *     multigroup data set.
	 *
	 * \param targetTemperature The single temperature value used when
	 *     computing opacities for each given density.
	 * \param densityFirst The beginning position of a STL
	 *     container that holds a list of densities.
	 * \param densityLast The end position of a STL
	 *     container that holds a list of densities.
	 * \param opacityFirst The beginning position of a STL
	 *     container into which opacity values corresponding to the
	 *     given density values will be stored.
	 * \return A list (of type OpacityIterator) of opacities are
	 *     returned.  These opacities correspond to the density
	 *     provided in the STL container and the single temperature value.
	 */
	template < class DensityIterator, class OpacityIterator >
		OpacityIterator getOpacity( double targetTemperature,
				DensityIterator densityFirst,
				DensityIterator densityLast,
				OpacityIterator opacityFirst ) const;

	/*!
	 * \brief Opacity accessor that returns a 2-D vector of opacities (
	 *     groups * bands ) that correspond to the
	 *     provided temperature and density. 
	 *
	 * \param targetTemperature The temperature value for which an
	 *     opacity value is being requested.
	 * \param targetDensity The density value for which an opacity 
	 *     value is being requested.
	 * \return A vector of opacities.
	 */
	std::vector< std::vector<double> > getOpacity( 
			double targetTemperature,
			double targetDensity ) const; 

	/*!
	 * \brief Opacity accessor that returns a vector of multigroupband
	 *     opacity 2-D vectors that correspond to the provided vector of
	 *     temperatures and a single density value.
	 *
	 * \param targetTemperature A vector of temperature values for
	 *     which opacity values are being requested.
	 * \param targetDensity The density value for which an opacity 
	 *     value is being requested.
	 * \return A vector of vectors of opacities.
	 */
	std::vector< std::vector< std::vector<double> > > getOpacity( 
			const std::vector<double>& targetTemperature,
			double targetDensity ) const; 

	/*!
	 * \brief Opacity accessor that returns a vector of 2-D vectors of
	 *     opacities that correspond to the provided vector of
	 *     densities and a single temperature value.
	 *
	 * \param targetTemperature The temperature value for which an 
	 *     opacity value is being requested.
	 * \param targetDensity A vector of density values for which
	 *     opacity values are being requested.
	 * \return A vector of vectors of opacities.
	 */
	std::vector< std::vector< std::vector<double> > > getOpacity( 
			double targetTemperature,
			const std::vector<double>& targetDensity ) const; 

	/*!
	 * \brief Query whether the data is in tables or functional form.
	 */
	bool data_in_tabular_form() const {
		return true;
	}; 

	/*!
	 * \brief Query to determine the reaction model.
	 */
	rtt_cdi::Reaction getReactionType() const {
		return opacityReaction;
	};

	/*!
	 * \brief Query to determine the physics model.
	 */
	rtt_cdi::Model getModelType() const {
		return opacityModel;
	};

	/*!
	 * \brief Returns a string that describes the templated
	 *     EnergyPolicy.  Currently this will return either "mg" or
	 *     "gray."
	 */ 
	std::string getEnergyPolicyDescriptor() const {
		return energyPolicyDescriptor;
	};

	/*!
	 * \brief Returns a "plain English" description of the opacity
	 *     data that this class references. (e.g. "Odfmg Rosseland
	 *     Scattering".) 
	 *
	 * The definition of this function is not included here to prevent 
	 *     the inclusion of the GandolfFile.hh definitions within this 
	 *     header file.
	 */
	std::string getDataDescriptor() const;

	/*!
	 * \brief Returns the name of the associated IPCRESS file.
	 *
	 * The definition of this function is not included here to prevent 
	 *     the inclusion of the GandolfFile.hh definitions within this 
	 *     header file.
	 */
	std::string getDataFilename() const;

	/*!
	 * \brief Returns a vector of temperatures that define the cached
	 *     opacity data table.
	 * 
	 * We do not return a const reference because this function
	 * must construct this information from more fundamental tables.
	 */
	std::vector< double > getTemperatureGrid() const;

	/*!
	 * \brief Returns a vector of densities that define the cached
	 *     opacity data table.
	 * 
	 * We do not return a const reference because this function
	 * must construct this information from more fundamental tables.
	 */
	std::vector< double > getDensityGrid() const;

	/*!
	 * \brief Returns a vector of energy values (keV) that define the
	 *     energy boundaries of the cached multigroup opacity data
	 *     table. 
	 */
	std::vector< double > getGroupBoundaries() const {
		return groupBoundaries;
	};

	/*!
	 * \brief Returns the size of the temperature grid.
	 */
	int getNumTemperatures() const;

	/*! 
	 * \brief Returns the size of the density grid.
	 */
	int getNumDensities() const;

	/*!
	 * \brief Returns the number of group boundaries found in the
	 *     current multigroup data set.
	 */
	int getNumGroupBoundaries() const {
		return groupBoundaries.size();
	};

	/*!
	 * \brief Returns the number of gruops found in the current
	 *     multigroup data set.
	 */
	int getNumGroups() const {
		return getNumGroupBoundaries() - 1;
	}; 

	/*!
	 * \brief Returns a vector of points along the cumulative opacity 
	 * 		distribution that mark the fraction of each band. First point is
	 * 		always zero, last point is always one.
	 */
	std::vector<double> getBandBoundaries() const {
		return bandBoundaries;
	};

	/*!
	 * \brief Returns the number of group boundaries found in the
	 *     current multigroup data set.
	 */
	int getNumBandBoundaries() const {
		return bandBoundaries.size();
	};

	/*!
	 * \brief Returns the number of band boundaries set in the
	 *     current multigroup data set. The lower boundary is always zero,
	 *     and the upper boundary is always one. Analagous to multigroup:
	 *     numBands = numBandBoundaries - 1.
	 */
	int getNumBands() const {
		return getNumBandBoundaries() - 1;
	};


	/*!
	 * \brief Pack a GandolfMulitgroupOpacity object.
	 *
	 * \return packed state in a vector<char>
	 */ 
	std::vector<char> pack() const;

	/*!
	 * \brief Returns the general opacity model type, defined in OpacityCommon.hh
	 *
	 * Since this is a Gandolf model, return 2 (rtt_cdi::GANDOLF_TYPE)
	 */
	rtt_cdi::OpacityModelType getOpacityModelType() const {
		return rtt_cdi::GANDOLF_TYPE;
	}
}; // end of class GandolfOdfmgOpacity

//---------------------------------------------------------------------------//
// INCLUDE TEMPLATE MEMBER DEFINITIONS FOR AUTOMATIC TEMPLATE INSTANTIATION
//---------------------------------------------------------------------------//

// --------------------------------- //
// STL-like accessors for getOpacity //
// --------------------------------- //

// ------------------------------------------ //
// getOpacity with Tuple of (T,rho) arguments //
// ------------------------------------------ //

template < class TemperatureIterator, class DensityIterator,
			class OpacityIterator >
				OpacityIterator GandolfOdfmgOpacity::getOpacity(
						TemperatureIterator tempIter, 
						TemperatureIterator tempLast,
						DensityIterator densIter, 
						DensityIterator densLast,
						OpacityIterator opIter ) const
{ 
	using std::log;

	// from twix:/scratch/tme/kai/KCC_BASE/include/algorithm

	// assert that the two input iterators have compatible sizes.
	Require( std::distance( tempIter, tempLast )
			== std::distance( densIter, densLast ) );

	// number of groups in this multigroup set.
	const int ng = spGandolfDataTable->getNumGroupBoundaries()-1;

	// temporary opacity vector used by the wrapper.  The returned 
	// data will be copied into the opacityIterator.
	std::vector<double> mgOpacity( ng );

	// loop over the (temperature,density) tuple.
	for ( ; tempIter != tempLast; ++tempIter, ++densIter )
	{
		// Call the Gandolf interpolator.
		// the vector opacity is returned.
		mgOpacity = wrapper::wgintmglog( 
				spGandolfDataTable->getLogTemperatures(),
				spGandolfDataTable->getNumTemperatures(), 
				spGandolfDataTable->getLogDensities(), 
				spGandolfDataTable->getNumDensities(),
				spGandolfDataTable->getNumGroupBoundaries(),
				spGandolfDataTable->getLogOpacities(), 
				spGandolfDataTable->getNumOpacities(),
				log( *tempIter ),
				log( *densIter ) );

		// The opacity vector contains the solution.  Now
		// we copy this solution into the OpacityIterator
		for ( int i=0; i<ng; ++i, ++opIter )
			*opIter = mgOpacity[i];
	}
	return opIter;
}

// ------------------------------------ // 
// getOpacity() with container of temps //
// ------------------------------------ // 

template < class TemperatureIterator, class OpacityIterator >
OpacityIterator GandolfOdfmgOpacity::getOpacity(
		TemperatureIterator tempIter,
		TemperatureIterator tempLast,
		double targetDensity,
		OpacityIterator opIter ) const
{ 
	using std::log;

	// number of groups in this multigroup set.
	const int ng = spGandolfDataTable->getNumGroupBoundaries()-1;

	// temporary opacity vector used by the wrapper.  The returned 
	// data will be copied into the opacityIterator.
	std::vector<double> mgOpacity( ng );

	// loop over the (temperature,density) tuple.
	for ( ; tempIter != tempLast; ++tempIter )
	{
		// Call the Gandolf interpolator.
		// the vector opacity is returned.
		mgOpacity = wrapper::wgintmglog(
				spGandolfDataTable->getLogTemperatures(),
				spGandolfDataTable->getNumTemperatures(), 
				spGandolfDataTable->getLogDensities(), 
				spGandolfDataTable->getNumDensities(),
				spGandolfDataTable->getNumGroupBoundaries(),
				spGandolfDataTable->getLogOpacities(), 
				spGandolfDataTable->getNumOpacities(),
				log( *tempIter ),
				log( targetDensity ) );

		// The opacity vector contains the solution.  Now
		// we copy this solution into the OpacityIterator
		for ( int i=0; i<ng; ++i, ++opIter )
			*opIter = mgOpacity[i];
	}
	return opIter;
}

// ---------------------------------------- // 
// getOpacity() with container of densities //
// ---------------------------------------- //

template < class DensityIterator, class OpacityIterator >
OpacityIterator GandolfOdfmgOpacity::getOpacity(
		double targetTemperature,
		DensityIterator densIter, 
		DensityIterator densLast,
		OpacityIterator opIter ) const
{ 
	using std::log;

	// number of groups in this multigroup set.
	const int ng = spGandolfDataTable->getNumGroupBoundaries()-1;

	// temporary opacity vector used by the wrapper.  The returned 
	// data will be copied into the opacityIterator.
	std::vector<double> mgOpacity( ng );

	// loop over the (temperature,density) tuple.
	for ( ; densIter != densLast; ++densIter )
	{
		// Call the Gandolf interpolator.
		// the vector opacity is returned.
		mgOpacity = wrapper::wgintmglog( 
				spGandolfDataTable->getLogTemperatures(),
				spGandolfDataTable->getNumTemperatures(), 
				spGandolfDataTable->getLogDensities(), 
				spGandolfDataTable->getNumDensities(),
				spGandolfDataTable->getNumGroupBoundaries(),
				spGandolfDataTable->getLogOpacities(), 
				spGandolfDataTable->getNumOpacities(),
				log( targetTemperature ),
				log( *densIter ) );

		// The opacity vector contains the solution.  Now
		// we copy this solution into the OpacityIterator
		for ( int i=0; i<ng; ++i, ++opIter )
			*opIter = mgOpacity[i];
	}
	return opIter;
}

} // end namespace rtt_cdi_gandolf

#endif // __cdi_gandolf_GandolfOdfmgOpacity_hh__

//---------------------------------------------------------------------------//
//                end of cdi_gandolf/GandolfOdfmgOpacity.hh
//---------------------------------------------------------------------------//
