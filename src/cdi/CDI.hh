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
 *     classes that are derived from the abstract classes found in the CDI
 *     package.  A CDI object is then created using these concrete classes 
 *     as constructor parameters.  Each CDI object will provide access to
 *     data for <b><i>one<i><b> material.
 * <p>
 *     Since this header file does not include the definitions for
 *     GrayOpacity or MultigroupOpacity, the calling routine must include
 *     these header files.  If the calling routine does not make use of
 *     one of these classes then it's definition file does not need to be
 *     included, however this will result in the compile-time warning:
 * <pre>
 *       line 69: warning: delete of pointer to incomplete class
 *  	    delete p;
 * </pre>
 *     The user should not worry about this warning as long he/she is not
 *     trying to instantiate the class specified by the error message.
 */

/*!
 * \example cdi/test/tCDI.cc
 *
 * \sa This test code provides an example of how to use CDI to access an
 *     user defined opacity class.  We have created an opacity class
 *     called dummyOpacity that is used in the creation of a CDI object.
 *     The CDI object is then used to obtain obacity data (via
 *     dummyOpacity).
 * <p>
 *     The test code also provides a mechanism to test the CDI independent 
 *     of any "real" data objects.
 *
 */
//===========================================================================//

class CDI 
{

    // NESTED CLASSES AND TYPEDEFS

    // DATA

    /*!
     * \brief Smart pointer to the GrayOpacity object.
     *
     * \sa spGrayOpacity is a smart pointer that links a CDI object to
     *     an GrayOpacity object (any type of gray opacity - Gandolf,
     *     EOSPAC, Analytic, etc.).  The pointer is established in the
     *     CDI constructor. 
     */
    const rtt_dsxx::SP< GrayOpacity > spGrayOpacity;

    /*!
     * \brief Smart pointer to the MultigroupOpacity object.
     *
     * \sa spMultigroupOpacity is a smart pointer that links a CDI
     *     object to an MultigroupOpacity object (any type of gray
     *     opacity - Gandolf, EOSPAC, Analytic, etc.).  The pointer is
     *     established in the CDI constructor. 
     */
    const rtt_dsxx::SP< MultigroupOpacity > spMultigroupOpacity;
    
    /*!
     * \brief Smart pointer to the EOS object.
     *
     * \sa spEOS is a smart pointer that links a CDI object to an 
     *     EOS object (any type of gray opacity - Gandolf, EOSPAC,  
     *     Analytic, etc.).  The pointer is established in the CDI
     *     constructor. 
     *
     * EOS objects have not yet been implemented in the CDI.
     */
    // const rtt_dsxx::SP<EOS> spEOS;
    
  public:

    // CREATORS
    
    /*!
     * \brief The CDI object instantiates a CDI object by hooking
     *        itself to Opacity, Nuclear, and EOS Data objects for a
     *        single material.
     *
     * \sa Currently, CDI only interfaces opacity data (either gray or 
     *     multigroup).  There are a number of constructors.  The
     *     variation in constructors allows CDI objects to be
     *     instatiated with different sets of components.
     * <p>
     *     This type of interface is somewhat cumbersome and possible
     *     should be replaced with another model.  One such model is
     *     the use of keyword arguments in the constructor.  (See
     *     B. Stroustrup, "The Design and Evolution of C++," Section
     *     6.5.1.) 
     * <p> 
     *     Currently, there is no mechanism in place to assert that
     *     the materials associated with each component object
     *     (GrayOpacity, MultigroupOpacity, EOS, etc.) are the the
     *     same.  This means it is possible to create a single CDI
     *     object that contains GrayOpacity data for Al,
     *     MultigroupOpacity data for C and EOS data for Be.  The user 
     *     must be careful to ensure that a single CDI object allows
     *     access to data for a single material!
     *
     * \param _spGrayOpacity A smart pointer object to a GrayOpacity
     *        class.  The GrayOpacity class must be derived from the
     *        abstract class found in the CDI package.
     *
     * \param _spMultigroupOpacity A smart pointer object to a
     *        MultigroupOpacity class.  The MultigroupOpacity class
     *        must be derived from the abstract class found in the CDI
     *        package. 
     *
     * \return A CDI object.  A CDI object will be able to access the 
     *         data for a single material
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

    /*!
     * \brief This fuction returns the GrayOpacity object.
     *
     * \sa This provides the CDI with the full functionality of the
     *     interface defined in GrayOpacity.hh.  For example, the host 
     *     code could make the following call:<br>
     * <pre>
     *     double newOp = spCDI1->gray()->getOpacity( 55.3, 27.4 );
     * </pre>
     */
    rtt_dsxx::SP< GrayOpacity > gray();

    /*!
     * \brief This fuction returns the MultigroupOpacity object.
     *
     * \sa This provides the CDI with the full functionality of the
     *     interface defined in MultigroupOpacity.hh.  For example,
     *     the host code could make the following call:<br>
     * <pre>
     *    int numGroups = spCDI1->mg()->getNumGroupBoundaries()
     * </pre>
     */
    rtt_dsxx::SP< MultigroupOpacity > mg();

  private:
    
    // IMPLEMENTATION
};

} // end namespace rtt_cdi

#endif // __cdi_CDI_hh__

//---------------------------------------------------------------------------//
// end of cdi/CDI.hh
//---------------------------------------------------------------------------//
