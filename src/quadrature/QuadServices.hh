//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   quadrature/QuadServices.hh
 * \author Kelly Thompson
 * \date   Mon Nov 8 11:17:12 2004
 * \brief  Provide Moment-to-Discrete and Discrete-to-Moment operations.
 * \note   Copyright 2004 The Regents of the University of California.
 *
 * Long description.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef quadrature_QuadServices_hh
#define quadrature_QuadServices_hh

#include "ds++/SP.hh"
#include "Quadrature.hh"

namespace rtt_quadrature
{

//===========================================================================//
/*!
 * \class QuadServices
 * \brief Provide Moment-to-Discrete and Discrete-to-Moment operations.
 *
 * This class provides the Moment-to-Discrete and Discrete-to-Moment matrices
 * as described by Jim Morel in "A Hybrid Collocation-Galerkin-Sn Method for
 * Solving the Boltzmann Transport Equation," Nuclear Science and
 * Engineering, Vol. 101, pp. 72-87 (2001).
 *
 * Some code was derived from the Partisn source file "galerkq.f."  Other
 * pieces were derived from the Snac2 soure file "sncYnm.f90."
 *
 * Background:
 *
 * Start with the steady state transport equation for neutrons,
 *
 * \f[
 * \mathbf\Omega\cdot\mathbf\nabla\psi(\mathbf{r},\mathbf{\Omega},E)
 * + \sigma_t \psi - Q(\mathbf{r},\mathbf{\Omega},E) =
 * \int\limits^\infty_0{dE' \int\limits_{4\pi} {d\mathbf\Omega'
 * \sigma_s(\mathbf{r},\mathbf\Omega'\rightarrow\mathbf\Omega,E'\rightarrow
 * E) \psi(\mathbf{r},\mathbf\Omega',E') }}.
 * \f]
 *
 * Here, \f$ Q(\mathbf{r},\mathbf{\Omega},E) \f$ represents external and
 * fission sources.  Now, approximate the integral over angles (the
 * scattering term) using a spherical harmonics treatment.
 * 
 * \f[
 * \mathbf\Omega\cdot\mathbf\nabla\psi(\mathbf{r},\mathbf{\Omega},E)
 * + (\sigma_a + \sigma_{s0}) \psi - Q(\mathbf{r},\mathbf{\Omega},E) 
 * = \int\limits^\infty_0{dE' \sum\limits_{k=0}^\infty{
 * \sum\limits_{\ell=-k}^k{\frac{2k+1}{4\pi} c_{\ell,k} Y_{\ell,k}(\mathbf\Omega)
 * \phi_{\ell,k}(\mathbf{r},E')\sigma_{sk}(\mathbf{r},E'\rightarrow E)}}}.
 * \f]
 *
 * We have defined the moment based scattering cross section and the flux
 * moments in terms of Legendre polynomials and spherical harmonics,
 *
 * \f[
 * \sigma_{sk}(\mathbf{r},E'\rightarrow E) = \int\limits_{-1}^1{d\mu_0
 * P_k(\mu_0)\sigma_s(\mathbf{r},\mu_0,E'\rightarrow E)}
 * \f]
 * and
 * \f[
 * \phi_{\ell,k}(\mathbf{r},E) = \int\limits_{4\pi}{d\mathbf\Omega' c_{\ell,k}
 * Y_{\ell,k}^*(\mathbf\Omega')\psi(\mathbf{r},\mathbf\Omega',E)}.
 * \f]
 *
 * In this equation, \f$ \sigma_s(\mu_0) \f$ is the cross section for
 * scattering through an angle whose cosine is \f$ \mu_0 \f$.  The parameter
 * \f$ c_{\ell,k} \f$ is a normalization factor that depends on the index
 * pair \f$ (\ell, k) \f$ and the normalization chosen for the definition of
 * \f$ Y_{\ell,k}(\mathbf\Omega) \f$.  The analytic scattering source
 * expanded with spherical harmonics can be expressed as,
 *
 * \f[
 * S(\mathbf{r},\mathbf\Omega,E) =
 * \int\limits_0^\infty{dE'\sum\limits_{\ell=0}^L{\sum\limits_{k=-\ell}^{\ell}{
 * \frac{2\ell+1}{4\pi} c_{\ell,k} Y_{\ell,k}(\mathbf\Omega)
 * \sigma_{s,\ell}(\mathbf{r},E'\rightarrow E) \int\limits_{4\pi}{d\mathbf\Omega'
 * c_{\ell,k} Y_{\ell,k}^*(\mathbf\Omega')\psi(\mathbf{r},\mathbf\Omega',E')}}}} 
 * \f]
 *
 * In this equation, we have also truncated the spherical harmonics expansion
 * at L, the order of the scattering matrix.  We can formulate a discrete
 * scattering source by using the discrete ordinate approximation,
 *
 * \f[
 * \mathbf{S}(\mathbf{r},E) = \int\limits_0^\infty{dE' \sum\limits_{\ell=0}^L{ 
 * \sum\limits_{k=-\ell}^{\ell}{
 * \frac{2\ell+1}{4\pi} c_{\ell,k} Y_{\ell,k}(\mathbf\Omega_m)
 * \sigma_{s,\ell}(\mathbf{r},E'\rightarrow E) 
 * \sum\limits_{m'=1}^M{w_{m'} c_{\ell,k} Y_{\ell,k}^*(\mathbf\Omega'_m)
 * \psi_m(\mathbf{r},E')}}}} 
 * \f]
 * or
 * \f[
 * \mathbf{S}(\mathbf{r},E) = 
 * \int\limits_0^\infty { dE' 
 * \left[ 
 * \mathbf{M}_{m,(\ell, k)} \cdot \mathbf\Sigma_{(\ell,k)}(E' \rightarrow E)
 * \cdot \mathbf{D}_{m,(\ell,k)} \cdot \psi_m(\mathbf{r},E')
 * \right]} \,\,\,\,\, \forall \,\,\,\,\, m = 1, \ldots, M, \,\,\,\,\, \ell = 0,
 * \ldots, L,\,\,\,\,\, and \,\,\,\,\, k = -\ell, \ldots, \ell.
 * \f]
 * 
 * \f$ \mathbf{M} \f$ is the moments-to-discrete operator, \f$ \mathbf{D} \f$
 * is the discrete-to-moments operator and \f$ \mathbf\Sigma \f$ is a matrix
 * containing the scattering cross section moments for \f$ L^{th} \f$ order
 * anisotropic scattering.  In this representation, M is the total number of
 * discrete angles.  
 *
 * Now, we define the moments-to-discrete operator to be
 *
 * \f[
 * \mathbf{M}_{m,(\ell,k)} = \frac{2\ell+1}{\sum\limits_{m'=1}^{M}{w_{m'}}} c_{\ell,k}
 * Y_{\ell,k}(\mathbf\Omega_m),
 * \,\,\,\,\, \forall \,\,\,\,\, m = 1, \ldots, M, \,\,\,\,\, \ell = 0,
 * \ldots, L,\,\,\,\,\, and \,\,\,\,\,k = -\ell, \ldots, \ell.
 * \f]
 *
 * This moment-to-discrete operator can be used to perform the following
 * operation:
 *
 * \f[
 * \psi_m = \sum\limits_{\ell=0}^L
 * \sum\limits_{k=-\ell}^\ell\mathbf{M}_{m,(\ell,k)}\Phi_\ell^k.
 * \f]
 * 
 * \section galerkin_formulation Galerkin Formulation
 *
 * For this class, we always take the discrete-to-moment operator
 * to be the inverse of the moment-to-discrete operator.  That is,
 *
 * \f[
 * \mathbf{D}_{m,(\ell,k)} = \mathbf{M}_{m,(\ell,k)}^{-1}.
 * \f]
 * 
 * Morel suggests the following formula for computing \f$ c_{\ell,k} \f$,
 *
 * \f[
 * c_{\ell,k} = \sqrt{(2-\delta_{k,0})\frac{ (\ell - |k|)!}{(\ell + |k|)!}},
 * \f]
 * where \f$ \delta_x \f$ is the Kronecker delta function.
 *
 * In Morel's Galerkin treatment, specific \f$ (\ell, k) \f$ moments are
 * chosen to produce a square matrix, \f$ \mathbf{M} \f$.  
 *
 * For 3D quadrature sets, the following rules are used:
 *
 * \li Use all moments defined by \f$
 * k=-\ell,\ldots,\ell\,\,\,\,\,\forall\,\,\,\,\,\ell=0,\ldots,N-1. \f$
 * \li Augment this list with all moments defined by \f$ k=-\ell,\ldots, -1
 * \,\,\,\,\,\forall\,\,\,\,\, \ell=N. \f$
 * \li Augment this list with all moments defined by \f$ k=1,\ldots, \ell, \,\,\,\,\,
 * k \f$ odd \f$ \,\,\,\,\,\forall\,\,\,\,\, \ell=N. \f$
 * \li Augment this list with all moments defined by \f$ k=-\ell,\ldots, -1, \,\,\,\,\,
 * k \f$ even \f$ \,\,\,\,\,\forall\,\,\,\,\, \ell=N+1. \f$
 *
 * For 2D quadrature sets, the following rules are used:
 *
 * \li Use all moments defined by \f$
 * k=0,\ldots,\ell \,\,\,\,\,\forall\,\,\,\,\,\ell=0,\ldots,N-1. \f$
 * \li Augment this list with all moments defined by \f$ k=1, \ldots, \ell,
 * \,\,\,\,\, k \f$ odd \f$
 * \,\,\,\,\,\forall\,\,\,\,\, \ell=N. \f$
 *
 * For 1D quadrature sets, the following rules are used:
 *
 * \li Use all moments defined by \f$
 * k=0 \,\,\,\,\,\forall\,\,\,\,\,\ell=0,\ldots,N-1. \f$
 *
 * In all of these algorithms, \f$ N \f$ is the order of the \f$ S_N \f$
 * quadrature set. 
 *
 *
 *
 * \section sph_and_leg_polys Spherical Harmonics and Legendre Polynomials
 * 
 * For our uses, we define the spherical harmonic function to have the
 * following form,
 *
 * \f[
 * Y_{\ell,k}(\mathbf\Omega_m) = Y_{\ell,k}(\mu_m,\omega_m)
 * \f]
 * and
 * \f[
 * Y_{\ell,k}(\mu_m,\omega_m) = P_{\ell,k}(\mu_m) cos(k\omega_m), \,\,\,\,\,
 * \forall \,\,\,\,\, k \ge 0,
 * \f]
 * \f[
 * Y_{\ell,k}(\mu_m,\omega_m) = P_{\ell,-k}(\mu_m) sin(-k\omega_m), \,\,\,\,\,
 * \forall \,\,\,\,\, k < 0,
 * \f]
 * where \f$ P_{\ell,k} \f$ are the associated Legendre functions of degree
 * \f$ \ell \f$ and order k.  Note, for 1D quadrature sets, \f$ k=0 \f$, and
 * the spherical harmonics functions become associated Legendre
 * polynomials.
 *
 * \f$\omega_m\f$ is the azimuthal angle
 * defined as
 * \f[
 * \omega_m = tan^{-1}\frac{\xi_m}{\eta_m}, \,\,\,\,\, \forall
 * \,\,\,\,\,\eta_m \ge 0,\,\,\,\,\,\xi_m \ge 0,
 * \f]
 * \f[
 * \omega_m = \pi - tan^{-1}\frac{\xi_m}{\eta_m}, \,\,\,\,\, \forall
 * \,\,\,\,\,\eta_m < 0,\,\,\,\,\,\xi_m \ge 0,
 * \f]
 * \f[
 * \omega_m = tan^{-1}\frac{\xi_m}{\eta_m} + \pi, \,\,\,\,\, \forall
 * \,\,\,\,\,\eta_m, < 0\,\,\,\,\,\xi_m<0,
 * \f]
 * \f[
 * \omega_m = 2\pi - tan^{-1}\frac{\xi_m}{\eta_m}, \,\,\,\,\, \forall
 * \,\,\,\,\,\eta_m \ge 0,\,\,\,\,\,\xi_m < 0.
 * \f]
 * For 2D quadrature sets, \f$ \xi \f$ can be computed from the other two
 * ordinates as
 * \f[
 * \xi = \sqrt{1 - \mu^2 - \eta^2}.
 * \f]
 *
 * \sa QuadServices.cc for detailed descriptions.
 *
 * Code Sample:
 * \code
 * \\ 
 * SP<Quadrature> spQuad = QuadCreator().quadCreate( 
 *                               QuadCreator::LevelSym, 4 );
 *
 * \\ Create a QuadServices object that is attached to a specific
 * \\ quadrature set. 
 * QuadServices myQS( spQuad );
 *
 * \endcode
 */
/*! 
 * \example quadrature/test/tQuadServices.cc 
 * 
 * description of example
 */
//===========================================================================//

class QuadServices 
{
  public:

    // NESTED CLASSES AND TYPEDEFS

    typedef std::pair< unsigned, int > lk_index;

    // CREATORS
    
    //! Default constructor assumes that only isotropic scattering is used. 
    explicit QuadServices( rtt_dsxx::SP< const Quadrature > const spQuad_ );
    
    //! Constructor that allows the user to pick the (k,l) moments to use.
    //! \todo This still needs to be defined.
    QuadServices( rtt_dsxx::SP< const Quadrature > const   spQuad_,
		  std::vector< lk_index >          const & lkMoments_ );

    //! Copy constructor (the long doxygen description is in the .cc file).
    QuadServices( QuadServices const & rhs );

    //! Destructor.
    virtual ~QuadServices() { /* empty */ }

    // MANIPULATORS
    
    //! Assignment operator for QuadServices.
    QuadServices& operator=( QuadServices const & rhs );

    // ACCESSORS

    //! \brief Return the moment-to-discrete operator.
    std::vector< double > getM() const { return Mmatrix; }

    //! \brief Return the discrete-to-moment operator.
    std::vector< double > getD() const { return Dmatrix; }

  private:

    // NESTED CLASSES AND TYPEDEFS

    // IMPLEMENTATION

    //! \bug move to ds++
    template< typename T > T kronecker_delta( T const, T const ) const;
    
    //! \bug move to ds++
    template< typename T > T factorial( T const ) const;

    template< typename T > 
    void print_matrix( std::string           const & matrix_name,
		       std::vector<T>        const & x,
		       std::vector<unsigned> const & dims ) const;
    
    //! Compute the (k,ell) spherical harmonic evaluated at Omega_m.
    double spherical_harmonic( unsigned const m,
			       unsigned const k,
			       int      const ell ) const;    

    //! \brief constuct maps between moment index n and the tuple (k,l).
    // This can optionally be provided by the user.
    std::vector< lk_index > compute_n2lk()    const;
    std::vector< lk_index > compute_n2lk_1D() const;
    std::vector< lk_index > compute_n2lk_2D() const;
    std::vector< lk_index > compute_n2lk_3D() const;

    //! Build the Mmatrix.
    std::vector< double > computeM() const;
    std::vector< double > computeD() const;

    //! Helper functions to compute coefficients
   inline double compute_clk( unsigned k, int ell ) const;
   inline double compute_azimuthalAngle( double const mu,
	  			         double const eta,
				         double const xi ) const;

    //! Checks
    bool D_equals_M_inverse() const;

    // DATA
    rtt_dsxx::SP< const Quadrature > const spQuad;
    unsigned                         const numMoments;
    std::vector< lk_index >          const n2lk;
    std::vector< double >            const Mmatrix;
    std::vector< double >            const Dmatrix;

};

} // end namespace rtt_quadrature

#include "QuadServices.i.hh"

#endif // quadrature_QuadServices_hh

//---------------------------------------------------------------------------//
//              end of quadrature/QuadServices.hh
//---------------------------------------------------------------------------//
