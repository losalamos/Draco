//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   cdi/CDI.hh
 * \author Kelly Thompson
 * \date   Thu Jun 22 16:22:06 2000
 * \brief  CDI class header file.
 * \note   Copyright © 2003 The Regents of the University of California.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef rtt_cdi_CDI_hh
#define rtt_cdi_CDI_hh

#include "GrayOpacity.hh"
#include "MultigroupOpacity.hh"
#include "EoS.hh"
#include "OpacityCommon.hh"
#include "ds++/SP.hh"
#include <vector>
#include <string>

namespace rtt_cdi
{
    
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
 *
 * Including CDI.hh in the client file will automatically include all of the
 * pertinent cdi headers and definitions (rtt_cdi::GrayOpacity,
 * rtt_cdi::MultigroupOpacity, EoS, rtt_cdi::Model, rtt_cdi::Reaction).
 * 
 * CDI also contains static services that allow a client to integrate the
 * (normalized) Planckian and Rosseland (\f$ \partial B/\partial T \f$)
 * integrals.  The overloaded functions that perform these services are:
 *
 * \arg integratePlanckSpectrum(const double lowFreq, const double highFreq,
 * const double T)
 *
 * \arg integratePlanckSpectrum(const double freq, const double T)
 *
 * \arg integrateRosselandSpectrum(const double lowf, const double hif,
 * double T);
 * 
 * \arg integrate_Rosseland_Planckian_Spectrum(const double lowf, const
 * double hif, const double T, double& PL, double& ROSL);
 *
 * \arg integratePlanckSpectrum(const int groupIndex, const double T);
 *
 * \arg integratePlanckSpectrum(const double T);
 *
 * \arg integrateRosselandSpectrum(const int groupIndex,  const double T);
 * 
 * \arg integrate_Rosseland_Planckian_Spectrum(const int groupIndex, const
 * double T, double& PL, double& ROSL);
 *
 * The first four forms can be called (CDI::integratePlanckSpectrum())
 * anytime.  They simply integrate the normalized Planckian or Rosseland over
 * a frequency range.  The next four forms may only be called after the
 * multigroup frequency boundaries have been set.  These boundaries are set
 * after a call, from any CDI object, to setMultigroupOpacity().  The
 * frequency boundaries are stored statically.  After an initial call by any
 * CDI object of setMultigroupOpacity(), the frequency bounds are checked to
 * make sure they do not change.  Changing the boundaries throws an
 * exception.  Thus, clients are allowed to view the group structure through
 * any multigroup data set (cdi.mg()->getGroupBoundaries()) or "globally" by
 * calling CDI::getFrequencyGroupBoundaries().  The context of usage dictates
 * which type of call to make; the result is invariant.  See the test
 * (tCDI.cc) for usage examples of the CDI Plankian and Rosseland integration
 * routines.
 *
 * We detail the two types of integration here, instead of in the individual
 * methods:
 *
 *
 * The Planckian functions integrate the normalized Plankian that is defined:
 *
 * \f[
 *    b(x) = \frac{15}{\pi^4} \frac{x^3}{e^x - 1}
 * \f]
 *
 * where 
 * 
 * \f[
 *    x = \frac{h\nu}{kT}
 * \f]
 *
 * and 
 *
 * \f[ 
 *    B(\nu,T)d\nu = \frac{acT^4}{4\pi} b(x)dx
 * \f]
 * 
 * where \f$B(\nu,T)\f$ is the Plankian and is defined
 *
 * \f[
 *    B(\nu,T) = \frac{2h\nu^3}{c^2} \frac{1}{e^{h\nu/kt} - 1}
 * \f]
 *
 * The normalized Plankian, integrated from 0 to \f$\infty\f$, equals
 * one. However, depending upon the maximum and minimum group boundaries, the
 * normalized Planck function may integrate to something less than one.
 *
 * This function performs the following integration:
 *\f[
 *      \int_{x_low}^{x_high} b(x) dx
 *\f]
 * where \f$x_low\f$ is calculated from the input low frequency bound and
 * \f$x_high\f$ is calculated from the input high frequency bound.  This
 * integration uses the method of B. Clark (JCP (70)/2, 1987).  We use a
 * 10-term Polylogarithmic expansion for the normalized Planckian, except in
 * the low-x limit, where we use a 21-term Taylor series expansion.
 *
 * The user is responsible for applying the appropriate constants to the
 * result of this integration.  For example, to make the result of this
 * integration equivalent to
 * \f[
 *      \int_{\nu_low}^{\nu_high} B(\nu,T) d\nu 
 * \f]
 * then you must multiply by a factor of \f$\frac{acT^4}{4\pi}\f$ where a is
 * the radiation constant.  If you want to evaluate expressions like the
 * following:
 *\f[
 *      \int_{4\pi} \int_{\nu_low}^{\nu_high} B(\nu,T) d\nu d\Omega
 *\f]
 * then you must multiply by \f$acT^4\f$.
 *
 * In the limit of \f$T \rightarrow 0, b(T) \rightarrow 0, therefore we
 * return a hard zero for a temperature equal to a hard zero.
 *
 * The integral is calculated using a polylogarithmic series approximation,
 * except for low frequencies, where the Taylor series approximation is more
 * accurate.  Each of the Taylor and polylogarithmic approximation has a
 * positive trucation error, so they intersect above the correct solution;
 * therefore, we always use the smaller one for a continuous concatenated
 * function.  When both frequency bounds reside above the Planckian peak
 * (above 2.822 T), we skip the Taylor series calculations and use the
 * polylogarithmic series minus one (the minus one is for roundoff control).
 *
 *
 *
 * This Rosseland functions integrate the normalized Rosseland that is defined:
 * \f[
 *    r(x) = \frac{15}{4\pi^4} \frac{x^4 e^x}{(e^x - 1)^2}
 * \f]
 * where 
 * \f[
 *    x = \frac{h\nu}{kT}
 * \f]
 * and 
 * \f[
 *    R(\nu, T)d\nu = \frac{4 acT^3}{4\pi}r(x)dx
 * \f]
 * where \f$R(\nu, T)\f$ is the Rosseland and is defined
 * \f[
 *    R(\nu, T) = \frac{\partial B(\nu, T)}{\partial T}
 * \f]
 * \f[
 *    B(\nu, T) = \frac{2h\nu^3}{c^2} \frac{1}{e^{\frac{h\nu}{kT}} - 1}
 * \f]
 * The normalized Rosseland, integrated from 0 to \f$\infty\f$, equals
 * one. However, depending upon the maximum and minimum group boundaries, the
 * normalized Rosseland function may integrate to something less than one.
 *
 * This function performs the following integration:
 * \f[
 *      \int_{x_1}^{x_N} r(x) dx
 * \f]
 * where \f$x_1\f$ is the low frequency bound and \f$x_N\f$ is the high
 * frequency bound of the multigroup data set.  This integration uses the
 * method of B. Clark (JCP (70)/2, 1987).  We use a 10-term Polylogarithmic
 * expansion for the normalized Planckian, except in the low-x limit, where
 * we use a 21-term Taylor series expansion.
 *
 * For the Rosseland we can relate the group interval integration to the
 * Planckian group interval integration, by equation 27 in B. Clark paper.
 * \f[
 *     \int_{x_1}^{x_N} r(x) dx = \int_{x_1}^{x_N} b(x) dx
 *     - \frac{15}{4\pi^4} \frac{x^4}{e^x - 1}
 * \f]
 * Therefore our Rosslenad group integration function can simply wrap the
 * Planckian function integratePlanckSpectrum(lowFreq, highFreq, T)
 * 
 * The user is responsible for applying the appropriate constants to the
 * result of this integration.  For example, to make the result of this
 * integration equivalent to
 * \f[
 *      \int_{\nu_1}^{\nu_N} R(\nu,T) d\nu
 * \f]
 * then you must multiply by a factor of \f$\frac{4 acT^3}{4\pi}\f$ where a is
 * the radiation constant.  If you want to evaluate expressions like the
 * following:
 *\f[
 *      \int_{4\pi} \int_{\nu_1}^{\nu_N} B(\nu,T) d\nu d\Omega
 *\f]
 * then you must multiply by \f$ 4 acT^3 \f$.
 *
 * In the limit of \f$ T \rightarrow 0, r(T) \rightarrow 0\f$, therefore we
 * return a hard zero for a temperature equal to a hard zero.
 *
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
     * is indexed [0,num_Models-1][0,num_Reactions-1].  It is accessed by
     * [rtt_cdi::Model][rtt_cdi::Reaction]. 
     */
    VF_MultigroupOpacity multigroupOpacities;

    /*!
     * \brief Frequency group boundaries for multigroup data.
     *
     * This is a static vector that contains the frequency boundaries for
     * multigroup data sets.  The number of frequency (energy) groups is the
     * size of the vector minus one.  
     *
     * This data is stored as static so that the same structure is guaranteed
     * for all multigroup data sets.  Thus, each CDI object will have access
     * to the same energy group structure.
     *
     */
    static std::vector<double> frequencyGroupBoundaries;
	
    /*!
     * \brief Smart pointer to the equation of state object.
     *
     * spEoS is a smart pointer that links a CDI object to an equation of
     * state object (any type of EoS - EOSPAC, Analytic, etc.).  The pointer
     * is established in the CDI constructor.
     */	     
    SP_EoS spEoS;

    //! Material ID.
    std_string matID;



    
  public:
	
    // STRUCTORS
    // ---------
	
    CDI(const std_string &id = std_string());
    virtual ~CDI();
	


    // SETTERS
    // -------

    //! Register a gray opacity (rtt_cdi::GrayOpacity) with CDI.
    void setGrayOpacity(const SP_GrayOpacity &spGOp);

    //! Register a multigroup opacity (rtt_cdi::MultigroupOpacity) with CDI.
    void setMultigroupOpacity(const SP_MultigroupOpacity &spMGOp);

    //! Register an EOS (rtt_cdi::Eos) with CDI.
    void setEoS(const SP_EoS &in_spEoS);

    //! Clear all data objects
    void reset();


    // GETTERS
    // -------
	
    SP_GrayOpacity       gray(rtt_cdi::Model m, rtt_cdi::Reaction r) const;
    SP_MultigroupOpacity mg  (rtt_cdi::Model m, rtt_cdi::Reaction r) const;
    SP_EoS eos() const;

    //! Return material ID string.
    const std_string& getMatID() const { return matID; }

    bool isGrayOpacitySet      (rtt_cdi::Model, rtt_cdi::Reaction) const;
    bool isMultigroupOpacitySet(rtt_cdi::Model, rtt_cdi::Reaction) const;
    bool isEoSSet() const;

    static std::vector<double> getFrequencyGroupBoundaries();

    static int getNumberFrequencyGroups();




    // INTEGRATORS:
    // -----------

    //! Integrate the normalized Planckian from 0 to x (hnu/kT).
    static double integratePlanckSpectrum(const double frequency, 
					  const double T); 

    //! Integrate the normalized Rosseland from 0 to x (hnu/kT).
    static double integrateRosselandSpectrum(const double frequency,
                                             const double T);

    //! Integrate the normalized Planckian and Rosseland from 0 to x (hnu/kT)
    static void integratePlanckRosselandSpectrum(const double frequency,
                                                 const double T,
                                                 double& placnk,
                                                 double& rosseland);




    //! Integrate the normalized Planckian over a frequency range.
    static double integratePlanckSpectrum(const double lowf,
                                          const double hif, 
					  const double T); 
    
    //! Integrate the normalized Planckian spectrum over a frequency group.
    static double integratePlanckSpectrum(const int groupIndex,
                                          const double T);

    //! Integrate the normalized Planckian spectrum over all frequency groups.
    static double integratePlanckSpectrum(const double T);
    
    //! Integrate the normalized Rosseland spectrum over a frequency group
    static double integrateRosselandSpectrum(const int groupIndex, 
					     const double T);

    //! Integrate the normalized Rosseland over a frequency range.  
    static double integrateRosselandSpectrum(const double lowf,
					     const double hif, 
					     const double T);

    //! Integrate the Planckian and Rosseland over a frequency range.
    static void integrate_Rosseland_Planckian_Spectrum(const double lowf,
						       const double hif,
						       const double T, 
						       double& planck, 
						       double& rosseland);

    //! Integrate the Planckian and Rosseland over a frequency group.
    static void integrate_Rosseland_Planckian_Spectrum(const int groupIndex,
						       const double T,
						       double& planck,
						       double& rosseland);

    //! Integrate the Planckian over all frequency groups
    static void integrate_Planckian_Spectrum(const std::vector<double>& bounds,
                                             const double T,
                                             std::vector<double>& planck);


    //! Integrate the Planckian and Rosseland over all frequency groups
    static void integrate_Rosseland_Planckian_Spectrum(const std::vector<double>& bounds,
                                                       const double T,
                                                       std::vector<double>& planck,
                                                       std::vector<double>& rosseland);

};


//---------------------------------------------------------------------------//
// INLINE FUNCTIONS
//---------------------------------------------------------------------------//


    
} // end namespace rtt_cdi

#endif // rtt_cdi_CDI_hh

//---------------------------------------------------------------------------//
// end of cdi/CDI.hh
//---------------------------------------------------------------------------//
