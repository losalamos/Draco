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

#include "OpacityCommon.hh"
#include "ds++/SP.hh"
#include <vector>
#include <string>

namespace rtt_cdi
{

class GrayOpacity;
class MultigroupOpacity;
class EoS;
    
//===========================================================================//
/*!
 * \class CDI
 *
 * \brief This class provides a Common Data Interface (CDI) to Atomic,
 * Nuclear and Equation of State (EOS) data.
 *
 * \sa The Draco packages cdi_gandolf (opacity data) and cdi_eospac (equation
 * of state data).
 *
 * The client must first instantiate concrete Opacity, Nuclear and EOS 
 * classes that are derived from the abstract classes found in the CDI
 * package.  A CDI object is then created using these concrete classes 
 * as constructor parameters.  Each CDI object will provide access to
 * data for <b><i>one</i></b> material.  This material may be a
 * mixture (e.g. water) if that mixture has been defined in the
 * underlying data tables.  However, CDI will not mix data table entries
 * to create a new material.  This type of mixing should be done by a
 * seperate package or the client code.
 * 
 * Since this header file does not include the definitions for
 * GrayOpacity or MultigroupOpacity, the calling routine must include
 * these header files.  If the calling routine does not make use of
 * one of these classes then it's definition file does not need to be
 * included, however this will result in the compile-time warning:
 *
 * <pre>
 *       line 69: warning: delete of pointer to incomplete class
 *  	    delete p;
 * </pre>
 *
 * The user should not worry about this warning as long he/she is not trying
 * to instantiate the class specified by the error message.
 */
/*!
 * \example cdi/test/tCDI.cc
 *
 * This test code provides an example of how to use CDI to access an user
 * defined opacity class.  We have created an opacity class called
 * dummyOpacity that is used in the creation of a CDI object.  The CDI object
 * is then used to obtain obacity data (via dummyOpacity).
 *
 * The test code also provides a mechanism to test the CDI independent of any
 * "real" data objects.
 */
//===========================================================================//

class CDI 
{
    // NESTED CLASSES AND TYPEDEFS
    typedef rtt_dsxx::SP<const GrayOpacity>       SP_GrayOpacity;
    typedef rtt_dsxx::SP<const MultigroupOpacity> SP_MultigroupOpacity;
    typedef rtt_dsxx::SP<const EoS>               SP_EoS;
    typedef std::vector<SP_GrayOpacity>           SF_GrayOpacity;
    typedef std::vector<SF_GrayOpacity>           VF_GrayOpacity;
    typedef std::vector<SP_MultigroupOpacity>     SF_MultigroupOpacity;
    typedef std::vector<SF_MultigroupOpacity>     VF_MultigroupOpacity;
    typedef std::string                           std_string;
    
    // DATA
	
    /*!
     * \brief Array that stores the matrix of possible GrayOpacity types.
     *
     * gray_opacities contains smart pointers that links a CDI object to a
     * GrayOpacity object (any type of gray opacity - Gandolf, EOSPAC,
     * Analytic, etc.).  The smart pointers is entered in the set functions.
     *
     * grayOpacities is indexed [0,num_Models-1][0,num_Reactions-1].  It is
     * accessed by [rtt_cdi::Model][rtt_cdi::Reaction]
     *
     */
    VF_GrayOpacity grayOpacities;
	
    /*!
     * \brief Array that stores the list of possible MultigroupOpacity types.
     *
     * multigroupOpacities contains a list of smart pointers to
     * MultigroupOpacity objects for different rtt_cdi::Reaction types.  It
     * is indexed [0,num_Reactions-1].  It is accessed by
     * [rtt_cdi::Reaction]. 
     */
    SF_MultigroupOpacity multigroupOpacities;
	
    /*!
     * \brief Smart pointer to the equation of state object.
     *
     * spEoS is a smart pointer that links a CDI object to an equation of
     * state object (any type of EoS - EOSPAC, Analytic, etc.).  The pointer
     * is established in the CDI constructor.
     */	     
    SP_EoS spEoS;

    /*!
     * \brief Material ID.
     */
    std_string matID;
	
  public:
	
    // CREATORS
	
    /*!
     * \brief Construct a CDI object.
     *
     * Builds a CDI object.  The opacity and eos objects that this holds must
     * be loaded using the set functions.  There is no easy way to guarantee
     * that all of the set objects point to the same material.  CDI does do
     * checking that only one of each Model:Reaction pair of opacity objects
     * are assigned; however, the user can "fake" CDI with different
     * materials if he/she is malicious enough.  
     *
     * CDI does allow a string material ID indicator.  It is up to the client
     * to ascribe meaning to the indicator.
     *
     * \param id string material id descriptor, this is defaulted to null
     */
    CDI(const std_string &id = std_string());
	
    /*!
     * \brief Destructor for CDI objects.
     *
     * We include a destructor for the CDI class so that, if another object
     * inherits from CDI, the derived object correctly destroys the CDI base
     * class.
     */
    virtual ~CDI();
	
    // ACCESSORS

    // "set" functions

    // Register a gray opacity (rtt_cdi::GrayOpacity) with CDI.
    void setGrayOpacity(const SP_GrayOpacity &spGOp);

    // Register a multigroup opacity (rtt_cdi::MultigroupOpacity) with CDI.
    void setMultigroupOpacity(const SP_MultigroupOpacity &spMGOp);

    //! Register an EOS (rtt_cdi::Eos) with CDI.
    void setEoS(const SP_EoS &in_spEoS);

    // "get" functions
	
    /*!
     * \brief This fuction returns a GrayOpacity object.
     *
     * This provides the CDI with the full functionality of the interface
     * defined in GrayOpacity.hh.  For example, the host code could make the
     * following call: <tt> double newOp = spCDI1->gray()->getOpacity(
     * 55.3, 27.4 ); </tt>
     *
     * The appropriate gray opacity is returned for the given model and
     * reaction type.
     *
     * \param m rtt_cdi::Model specifying the desired physics model
     * \param r rtt_cdi::Reaction specifying the desired reaction type
     */
    SP_GrayOpacity gray(rtt_cdi::Model m, rtt_cdi::Reaction r) const;

    /*!
     * \brief This fuction returns the MultigroupOpacity object.
     *
     * This provides the CDI with the full functionality of the interface
     * defined in MultigroupOpacity.hh.  For example, the host code could
     * make the following call:<br> <tt> int numGroups =
     * spCDI1->mg()->getNumGroupBoundaries(); </tt>
     *
     * The appropriate multigroup opacity is returned for the given reaction
     * type.
     *
     * \param r rtt_cdi::Reaction specifying the desired reaction type.
     *
     */
    SP_MultigroupOpacity mg(rtt_cdi::Reaction r) const;
	
    /*!
     * \brief This fuction returns the EoS object.
     *
     * This provides the CDI with the full functionality of the interface
     * defined in EoS.hh.  For example, the host code could make the
     * following call:<br> <tt> double Cve =
     * spCDI1->eos()->getElectronHeatCapacity( * density, temperature );
     * </tt>
     */
    SP_EoS eos() const;

    /*!
     * \brief Return material ID string.
     */
    const std_string& getMatID() const { return matID; }

    /*!
     * \brief Reset the CDI object.
     *
     * This function "clears" all data objects (GrayOpacity,
     * MultigroupOpacity, EoS) held by CDI.  After clearing, new objects can
     * be set using the set functions.  
     *
     * As stated in the set functions documentation, you are not allowed to
     * overwrite a data object with the same attributes as one that already
     * has been set.  The only way to "reset" these objects is to call
     * CDI::reset().  Note that CDI::reset() resets \b ALL of the objects
     * stored by CDI.
     */
    void reset();

  private:
	
    // IMPLEMENTATION
};
    
} // end namespace rtt_cdi

#endif // __cdi_CDI_hh__

//---------------------------------------------------------------------------//
// end of cdi/CDI.hh
//---------------------------------------------------------------------------//
