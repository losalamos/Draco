//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi/Opacity.hh
 * \author Kelly Thompson
 * \date   Thu Jun 23 13:55:06 2000
 * \brief  Opacity header file (an abstract class).
 */
//---------------------------------------------------------------------------//
// $Ids: Opacity.hh,v 1.4 2000/07/24 22:49:20 kellyt Exp $
//---------------------------------------------------------------------------//

#ifndef __cdi_Opacity_hh__
#define __cdi_Opacity_hh__

#include <vector>
#include <string>

namespace rtt_cdi
{
 
using std::string;
using std::vector;

//===========================================================================//
/*!
 * \class Opacity
 *
 * \brief This is an abstract base class that CDI opacity objects must be
 *        derived from.
 *
 * \sa Any opacity object that is to be used in conjunction with the CDI
 * must be derived from this abstract base class.  The CDI class must
 * know about the virtual functions found in this class.  Any derived
 * object must define all of these functions.
 *
 */

// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

class Opacity
{

    // NESTED CLASSES AND TYPEDEFS

    // DATA

  public:

    // CREATORS
    
    // defaulted Opacity() 
    // defaulted Opacity(const CDI &rhs);
    
    /*!
     * \brief The opacity destructor must be defined here so that we
     *        do not create a memory leak (the virtual opacity class
     *        must be destroyed when the concrete base class is
     *        destroyed). 
     */
    virtual ~Opacity() {}

    // MANIPULATORS
    
    // defaulted CDI& operator=(const CDI &rhs);

    // ACCESSORS
    
    /*!
     * \breif Returns the opacity data filename.
     */
    virtual string const getDataFilename() = 0;

    /*!
     * \breif Return a vector material ID's found in the opacity data
     *        file.
     */
    virtual vector<int> const getMatIDs() = 0;

    /*!
     * \breif Returns a single opacity value for the prescribed
     *        temperature and density.
     *
     * This opacity object only knows how to access the data for one 
     * material.  The material identification is specified in the
     * construction of the derived concrete opacity object. 
     *
     * \param targetTemperature The temperature (in keV) of the
     *                          material. 
     * \param targetDensity The density (in g/cm^3) of the material.
     *
     * \return Gray opacity value for the current material at
     *         targetTemperature keV and targetDensity g/cm^3.
     */
    virtual double getGray( const double targetTemperature, 
			    const double targetDensity ) = 0;       

    /*!
     * \breif Returns a vector of the opacity values for each energy
     *        group for the prescribed temperature and density.
     *
     * The opacity object only knows how to access the data for one
     * material.  The material identification is specified in the
     * construction of the derived concrete opacity object. 
     *
     * \param targetTemperature The temperature (in keV) of the
     *                          material.
     * \param targetDensity The density (in g/cm^3) of the material.
     *
     * \return A vector of opacity values for the current material at 
     *         targetTemperature keV and targetDensity g/cm^3.  The
     *         vector has ngroups entries.  The number of groups is
     *         specified by the data file. 
     */
    virtual vector<double> getMG( const double targetTemperature, 
				  const double targetDensity ) = 0;       

    //  protected:     

    // DATA

  private:
    
    // IMPLEMENTATION
};

} // end namespace rtt_cdi

#endif // __cdi_Opacity_hh__

//---------------------------------------------------------------------------//
//                              end of cdi/Opacity.hh
//---------------------------------------------------------------------------//
