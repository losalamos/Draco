//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   quadrature/Quadrature.hh
 * \author Kelly Thompson
 * \date   Tue Feb 22 10:21:50 2000
 * \brief  Quadrature class header file.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __quadrature_Quadrature_hh__
#define __quadrature_Quadrature_hh__

#include <vector>
#include <string>

#include "ds++/Assert.hh"
// Require(cond) -- check input vals
// Ensure(cond)  -- check output vals
// Check(cond) or Assert(cond) -- other checks
// Insist(cond,msg) -- Always checks, others can be turned off with compile
//                     option.

namespace rtt_quadrature
{

using std::vector;
using std::string;

//! The numerical constant PI.
const double PI = 3.14159265358979323846;

//! Tol specifies how accurately we compare two double values.
const double TOL = 1.0e-13;

//===========================================================================//
/*!
 * \class Quadrature
 *
 * \brief A class to encapsulate the angular discretization.
 *
 * The Quadrature class provides services related to the angular
 * discretization scheme.  It creates a set of quadrature directions
 * (abscissas) and weights associated with a particular quadrature scheme
 * specified by the calling routine.
 *
 */

/*! 
 * \example quadrature/test/tQuadrature.cc
 * 
 * Example of Quadrature usage and testing algorithm.  This test code
 * generates several different quadrature sets and demonstrates access to the 
 * member data and functions.
 *
 * A QuadCreator object must be instatiated.  This object is
 * responsible for returning a Smart Pointer to the new Quadrature
 * object. 
 *
 */
// revision history:
// -----------------
// 0) original
// 1) Added lots of comments.
//    Changed "getmu()" to "getMu()" because Tycho was already using
//       the 2nd convention. 
//    Added several accessors: getEta(), getXi(), getWt(), getMu(int),
//       getEta(int), getXi(int), getWt(int), getOmega(int), name(),
//       dimensionality(),getSnOrder(). 
//    Implemented clients ability to specify norm == sumwt.
//    Added checks to verify that 
//       Integral(dOmega) = norm
//       Integral(Omega*dOmega) = (0,0,0)
//       Integral(Omega*Omega*dOmega) = norm/3 * ((1,0,0),(0,1,0),(0,0,1))
//    Added use of the ds++/Assert class.
// 
//===========================================================================//

class Quadrature 
{
  public:

    // CREATORS

    /*!
     * \brief The default constructor for the quadrature class.
     *
     * The default constructor sets the SN Order and normalization values.
     * This constructor complements the concrecte constructor for the
     * quadrature class being instantiated.  The concrete class will also
     * initialize the variable numAngles to an appropriate value.
     *
     * \param snOrder_ Integer specifying the order of the SN set to be
     *                 constructed.  Number of angles = (snOrder+2)*snOrder.
     * \param norm_    A normalization constant.  The sum of the quadrature
     *                 weights will be equal to this value.  Its default
     *                 value is set in QuadCreator.
     */

    Quadrature( int snOrder_, double norm_ )
	: snOrder( snOrder_ ), norm( norm_ ) { }

    // ACCESSORS

    /*!
     * \brief Return the mu vector.
     *
     * Retrieves a vector containing the first (mu) component of the
     * direction vector.  The direction vector, Omega, in general has three
     * components (mu, eta, xi).  
     *
     * Omega = mu * ex + eta * ey + xi * ey.
     *
     * mu  = cos(theta)*sin(phi)
     * eta = sin(theta)*sin(phi)
     * xi  = cos(phi)
     *
     * theta is the azimuthal angle.
     * phi is the polar angle.
     */
    const vector<double>& getMu() const { return mu; }

    /*!
     * \brief Return the eta vector.
     * 
     * Retrieves a vector whose elements contain the second (eta) components
     * of the direction vector.  If this function is called for a 1D set an
     * error will be issued.
     *
     * See comments for getMu().
     */
    const vector<double>& getEta() const 
    {
	// The quadrature set must have at least 2 dimensions to return eta.
	Require( dimensionality() >= 2 );
	return eta;
    }

    /*!
     * \brief Return the xi vector.
     *
     * Retrieves a vector whose elements contain the third (xi) components of
     * the direction vector.  If this function is called for a 1D or a 2D set 
     * an error will be issued.
     *
     * See comments for getMu().
     */
    const vector<double>& getXi() const 
    {
	Require( dimensionality() >= 3 );
	return xi;
    }

    /*!
     * \brief Return the wt vector.
     *
     * Retrieves a vector whose elements contain the weights associated with
     * each direction omega_m. 
     *
     * See comments for getMu().
     */
    const vector<double>& getWt() { return wt; }

    /*!
     * \brief Return the mu component of the direction Omega_m.
     *
     * Retrieves the m-th element of the mu component of the direction
     * vector. 
     *
     * See comments for getMu().
     *
     * \param m The direction index must be contained in the range
     *          (0,numAngles). 
     */
    // const double getMu ????
    double getMu( const int m ) const
    {
	// Angle index m must be greater than zero and less than numAngles.
	Require( m >= 0 && m < getNumAngles() );             
	// Die if the vector mu appears to be the wrong size.
	Check( mu.size() >= m );       
	return mu[m];
    }

    /*!
     * \brief Return the eta component of the direction Omega_m.
     *
     * Retrieves the m-th element of the eta component of the direction
     * vector.  If this accessor is called for a 1D set an error will be
     * issued. 
     *
     * See comments for getMu().
     *
     * \param m The direction index must be contained in the range
     *          (0,numAngles). 
     */
    double getEta( const int m ) const
    {
	// The quadrature set must have at least 2 dimensions to return eta.
	Require( dimensionality() >= 2 );
	// Angle index m must be greater than zero and less than numAngles.
	Require( m >= 0 && m < getNumAngles() ); 
	return eta[m];
    }

    /*!
     * \brief Return the xi component of the direction Omega_m.
     * 
     * Retrieves the m-th element of the xi component of the direction
     * vector.  If this accessor is called for a 1D or a 2D set an error will 
     * be issued.
     *
     * See comments for getMu().
     *
     * \param m The direction index must be contained in the range
     *          (0,numAngles). 
     */
    double getXi( const int m ) const
    {
	// The quadrature set must have at least 3 dimensions to return xi.
	Require( dimensionality() >= 3 );
	// Angle index m must be greater than zero and less than numAngles.
	Require( m >= 0 && m < getNumAngles() ); 
	return xi[m];
    }

    /*!
     * \brief Return the weight associated with the direction Omega_m.
     * 
     * Retrieves the weight associated with the m-th element of the direction
     * vector. 
     *
     * See comments for getMu().
     *
     * \param m The direction index must be contained in the range
     *          (0,numAngles). 
     */
    double getWt( const int m ) const
    {
	// Angle index m must be greater than zero and less than numAngles.
	Require( m >= 0 && m < getNumAngles() ); 
	return wt[m];
    }

    /*!
     * \brief Returns the Omega vector for all directions.
     *
     * Returns a vector of length numAngles.  Each entry is a length three 
     * vector of doubles that together represent the m-th discrete
     * direction. 
     */
    const vector< vector<double> > &getOmega() const 
    {
	return omega;
    }

    /*!
     * \brief Returns the Omega vector for a single direction.
     *
     * Returns the Omega direction vector associated with the client
     * specified index.
     *
     * \param m The direction index must be contained in the range
     * (0,numAngles). 
     */
    const vector<double> &getOmega( const int m ) const
    {
	Require ( m>=0 && m < getNumAngles() );
	return omega[m];
    }

    /*!
     * \brief Returns the number of directions in the current quadrature set.
     */
    virtual int getNumAngles() const = 0;

    /*!
     * \brief Prints a table containing all quadrature directions and weights.
     */
    virtual void display() const = 0;

    /*!
     * \brief The sum of the quadrature weights will be normalized so
     *        that they sum to this value.
     */
    double getNorm() const { return norm; }

    /*!
     * \brief Returns a string containing the name of the quadrature set.
     */
    virtual string name() const = 0;

    /*!
     * \brief Returns an integer containing the dimensionality of the quadrature set.
     */
    virtual int dimensionality() const = 0;

    /*!
     * \brief Returns an integer containing the Sn order of the quadrature set.
     */
    virtual int getSnOrder() const = 0;

    /*!
     * \brief Integrates dOmega over the unit sphere. (The sum of quadrature weights.)
     */
    double iDomega() const;

    /*!
     * \brief Integrates the vector Omega over the unit sphere. 
     *
     * The solution to this integral is a vector with length equal to the
     * number of dimensions of the quadrature set.  The solution should be
     * the zero vector. 
     *
     * The integral is actually calculated as a quadrature sum over all
     * directions of the quantity: 
     *
     *     Omega(m) * wt(m)
     *
     * Omega is a vector quantity.
     */
    vector<double> iOmegaDomega() const;

    /*!
     * \brief Integrates the tensor (Omega Omega) over the unit sphere. 
     *
     * The solution to this integral is a tensor vector with ndims^2 elements
     * The off-diagonal elements should be zero.  The diagonal elements
     * should have the value sumwt/3.  
     *
     * We actually return a 1D vector whose length is ndims^2.  The "diagonal"
     * elements are 0, 4 and 8.
     *
     * The integral is actually calculated as a quadrature sum over all
     * directions of the quantity:
     *
     *     Omega(m) Omega(m) wt(m)
     *
     * The quantity ( Omega Omega ) is a tensor quantity.
     */
    vector<double> iOmegaOmegaDomega() const;


    // Other accessors that may be needed in the future.
    //
    // virtual int getNumAnglesPerOctant() = 0;

  protected:

    // DATA

    const int snOrder; // defaults to 4.
    const double norm; // 1D: defaults to 2.0.
                       // 2D: defaults to 2*pi.
                       // 3D: defaults to 4*pi.

    // Quadrature directions and weights.
    vector<double> mu;
    vector<double> eta; // will be an empty vector for all 1D sets.
    vector<double> xi;  // will be an empty vector for all 1D and 2D sets.
    vector<double> wt;
    vector< vector< double > > omega;
};

//===========================================================================//
/*!
 * \class Q1DGaussLeg
 *
 * \brief A class to encapsulate a 1D Gauss Legendre Quadrature set.
 *
 * The client only needs to call QuadCreator::QuadCreate with the requested
 * quadrature set specified.  Since this quadrature set is inheireted from
 * the class Quadrature the client never needs to access this class directly.
 * The class overrides the virutal member functions contained in the class
 * Quadrature and contains data members that define the quadrature set.
 */
//===========================================================================//

class Q1DGaussLeg : public Quadrature
{

    // NESTED CLASSES and TYPEDEFS

    // DATA

    int numAngles;  // == snOrder

  public:

    // CREATORS

    /*!
     * \brief Constructs a 1D Gauss Legendre quadrature object.
     *
     * \param snOrder_ Integer specifying the order of the SN set to be
     *                 constructed. 
     * \param norm_    A normalization constant.  The sum of the quadrature
     *                 weights will be equal to this value (default = 2.0).
     */
    // The default values for snOrder_ and norm_ were set in QuadCreator.
    Q1DGaussLeg( int snOrder_, double norm_ );
    // Defaulted:    Q1DGaussLeg(const Q1DGaussLeg &rhs);
    // Defaulted:   ~Q1DGaussLeg();

    // MANIPULATORS
    
    // Defaulted:    Q1DGaussLeg& operator=(const Q1DGaussLeg &rhs);

    // ACCESSORS

    // These functions override the virtual member functions specifed in the
    // parent class Quadrature.
    
    int getNumAngles()   const { return numAngles; }
    void display()       const;
    string name()        const { return "1D Gauss Legendre"; }
    int dimensionality() const { return 1; }
    int getSnOrder()     const { return snOrder; }

  private:
    
    // IMPLEMENTATION

};


//===========================================================================//
/*!
 * \class Q2DLevelSym
 *
 * \brief A class to encapsulate a 2D Level Symmetric quadrature set.
 *
 * The client only needs to call QuadCreator::QuadCreate with the requested
 * quadrature set specified.  Since this quadrature set is inherited from
 * the class Quadrature the client never needs to access this class directly.
 * The class overrides the virtual member functions contained in the class
 * Quadrature and contains data members that define the quadrature set.
 */
//===========================================================================//

class Q2DLevelSym : public Quadrature
{

    // NESTED CLASSES and TYPEDEFS

    // DATA

    int numAngles; // defaults to 12.

  public:

    // CREATORS
    
    /*!
     * \brief Constructs a 2D Level Symmetric quadrature object.
     *
     * \param snOrder_ Integer specifying the order of the SN set to be
     *                 constructed.  Number of angles = (snOrder+2)*snOrder/2.
     * \param norm_    A normalization constant.  The sum of the quadrature
     *                 weights will be equal to this value (default = 2*PI).
     */
    // The default values for snOrder_ and norm_ were set in QuadCreator.
    Q2DLevelSym( int snOrder_, double norm_ );
    // Defaulted:   Q2DLevelSym(const Q2DLevelSym &rhs);
    // Defaulted:  ~Q2DLevelSym();

    // MANIPULATORS
    
    // Defaulted:   Q2DLevelSym& operator=(const Q2DLevelSym &rhs);

    // ACCESSORS

    // These functions override the virtual member functions specifed in the
    // parent class Quadrature.

    //! Returns the number of angles in the current quadrature set.
    int getNumAngles()   const { return numAngles; }
    //! Prints a short table containing the quadrature directions and weights.
    void display()       const;
    //! Returns the official name of the current quadrature set.
    string name()        const { return "2D Level Symmetric"; }
    //! Returns the number of dimensions in the current quadrature set.
    int dimensionality() const { return 2; }
    //! Returns the order of the SN set.
    int getSnOrder()     const { return snOrder; }
    //! Returns the number of eta levels in the quadrature set.
    int getLevels()      const { return snOrder; }

  private:
    
    // IMPLEMENTATION

};


//===========================================================================//
/*!
 * \class Q3DLevelSym
 *
 * \brief A class to encapsulate a 3D Level Symmetric quadrature set.
 *
 * The client only needs to call QuadCreator::QuadCreate with the requested
 * quadrature set specified.  Since this quadrature set is inheireted from
 * the class Quadrature the client never needs to access this class directly.
 * The class overrides the virutal member functions contained in the class
 * Quadrature and contains data members that define the quadrature set.
 */
//===========================================================================//

class Q3DLevelSym : public Quadrature
{

    // NESTED CLASSES and TYPEDEFS

    // DATA

    int numAngles; // defaults to 24.

  public:

    // CREATORS
    
    /*!
     * \brief Constructs a 3D Level Symmetric quadrature object.
     *
     * \param snOrder_ Integer specifying the order of the SN set to be
     *                 constructed.  Number of angles = (snOrder+2)*snOrder.
     * \param norm_    A normalization constant.  The sum of the quadrature
     *                 weights will be equal to this value (default = 4*PI).
     */
    // The default values for snOrder_ and norm_ were set in QuadCreator.
    Q3DLevelSym( int snOrder_, double norm_ );
    // Defaulted:   Q3DLevelSym(const Q3DLevelSym &rhs);
    // Defaulted:  ~Q3DLevelSym();

    // MANIPULATORS
    
    // Defaulted:   Q3DLevelSym& operator=(const Q3DLevelSym &rhs);

    // ACCESSORS

    // These functions override the virtual member functions specifed in the
    // parent class Quadrature.

    //! Returns the number of angles in the current quadrature set.
    int getNumAngles()   const { return numAngles; }
    //! Prints a short table containing the quadrature directions and weights.
    void display()       const;
    //! Returns the official name of the current quadrature set.
    string name()        const { return "3D Level Symmetric"; }
    //! Returns the number of dimensions in the current quadrature set.
    int dimensionality() const { return 3; }
    //! Returns the order of the SN set.
    int getSnOrder()     const { return snOrder; }
    //! Returns the number of xi levels in the quadrature set.
    int getLevels()      const { return snOrder; }

  private:
    
    // IMPLEMENTATION

};


} // end namespace rtt_quadrature

#endif                          // __quadrature_Quadrature_hh__

//---------------------------------------------------------------------------//
//                              end of quadrature/Quadrature.hh
//---------------------------------------------------------------------------//
