//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi/CDI.hh
 * \author Kelly Thompson
 * \date   Thu Jun 22 16:22:06 2000
 * \brief  CDI class header file.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __cdi_CDI_hh__
#define __cdi_CDI_hh__

#include "ds++/SP.hh"

namespace rtt_cdi
{

    class GrayOpacity;
    class MultigroupOpacity;

//===========================================================================//
/*!
 * \class CDI
 *
 * \brief This class provides a Common Data Interface (CDI) to Atomic, 
 *        Nuclear and Equation of State (EOS) data.
 *
 * \sa The client must first instantiate concrete Opacity, Nuclear and EOS 
 * classes that are derived from abstrat classes found in the CDI
 * package.  A CDI object is then created using these concrete classes 
 * as constructor parameters.  Each CDI object will provide access to
 * data for one material.
 *
 * Since this header file does not include the definitions for
 * GrayOpacity or MultigroupOpacity, the calling routine must include
 * these header files.  If the calling routine does not make use of
 * one of these classes then it's definition file does not need to be
 * included, however this will result in the compile-time warning:
 *
 *       line 69: warning: delete of pointer to incomplete class
 *  	    delete p;
 *
 * The user should not worry about this warning as long he/she is not
 * trying to instantiate the class specified by the error message.
 */

/*!
 * \example cdi/test/tCDI.cc
 *
 * This test code provides an example of how to use CDI to access an
 * user defined opacity class.  We have created an opacity class
 * called dummyOpacity that is used in the creation of a CDI object.
 * The CDI object is then used to obtain obacity data (via
 * dummyOpacity).
 *
 * The test code also provides a mechanism to test the CDI independent 
 * of any "real" data objects.
 *
 */

// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

class CDI 
{

    // NESTED CLASSES AND TYPEDEFS

    // DATA

    /*!
     * \brief Smart pointer to the opacity object.
     *
     * spOpacity is a smart pointer that links a CDI object to an
     * opacity object (any type of opacity - Gandolf, EOSPAC,
     * Analytic, etc.).  The pointer is established in the CDI
     * constructor. 
     *
     */
    const rtt_dsxx::SP< GrayOpacity > spGrayOpacity;
    const rtt_dsxx::SP< MultigroupOpacity > spMultigroupOpacity;
    
    // future expansion
    // const rtt_dsxx::SP<EOS> spEOS;
    
  public:

    // CREATORS
    
    /*!
     * \brief The CDI object instantiates a CDI object by hooking
     *        itself to Opacity, Nuclear, and EOS Data objects.
     *
     * \sa Currently, CDI only interfaces opacity data (either gray or 
     *     multigroup).
     *
     * \param _spOpacity A smart pointer object to an opacity class.
     *                   The opacity class must be derived from the
     *                   abstract class found in the CDI package.
     * \return A CDI object.  A CDI object will be able to access the 
     *         data for a single material
     *
     * We may need to consider keyword arguments in the constructor to 
     * avoid requiring a large number of constructors.  See
     * B. Stroustrup, "The Design and Evolution of C++," Section 6.5.1.
     */
     CDI( const rtt_dsxx::SP< GrayOpacity > _spGrayOpacity );
     CDI( const rtt_dsxx::SP< MultigroupOpacity > _spMultigroupOpacity );
     CDI( const rtt_dsxx::SP< GrayOpacity > _spGrayOpacity, 
	  const rtt_dsxx::SP< MultigroupOpacity > _spMultigroupOpacity );
    
    /*!
     * \brief Destructor for CDI objects.
     *
     * \sa We include a destructor for the CDI class so that, if
     *     another object inherits from CDI, the derived object
     *     correctly destroys the CDI base class.
     */
    virtual ~CDI() {};

    // ACCESSORS

    rtt_dsxx::SP< GrayOpacity > gray();
    rtt_dsxx::SP< MultigroupOpacity > mg();

  private:
    
    // IMPLEMENTATION
};

} // end namespace rtt_cdi

#endif // __cdi_CDI_hh__

//---------------------------------------------------------------------------//
// end of cdi/CDI.hh
//---------------------------------------------------------------------------//
