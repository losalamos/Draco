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

namespace rtt_quadrature
{

using std::vector;
using std::string;

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
 * \example quadrature/test/tQuadrature.cc
 * 
 * Example of Quadrature usage and testing algorithm.  This test code
 * generates several different quadrature sets and demonstrates access to the 
 * member data and functions.
 *
 */
// revision history:
// -----------------
// 0) original
// 1) 3-17-2000 a) Added lots of comments.
//              b) Changed "getmu()" to "getMu()" because Tycho was already 
//                 using the 2nd convention.
//              c) Added several accessors: getEta(), getXi(), getwt(), 
//                 getMu(int), getEta(int), getXi(int), getwt(int),
//                 getOmega(int), name(), dimensionality(),getSnOrder().
//              d) Implemented clients ability to specify norm == sumwt.
//              e) Added checks to verify that 
//                 Integral(dOmega) = sumwt
//                 Integral(Omega*dOmega) = (0,0,0)
//                 Integral(Omega*Omega*dOmega) 
//                    = norm/3 * ((1,0,0),(0,1,0),(0,0,1))
//              f) Added use of the ds++/Assert class.
//
// 
//===========================================================================//

class Quadrature 
{
  public:

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
    virtual const vector<double>& getMu() = 0;

    /*!
     * \brief Return the eta vector.
     * 
     * Retrieves a vector whose elements contain the second (eta) components
     * of the direction vector.  If this function is called for a 1D set an
     * error will be issued.
     *
     * See comments for getMu().
     */
    virtual const vector<double>& getEta() = 0;

    /*!
     * \brief Return the xi vector.
     *
     * Retrieves a vector whose elements contain the third (xi) components of
     * the direction vector.  If this function is called for a 1D or a 2D set 
     * an error will be issued.
     *
     * See comments for getMu().
     */
    virtual const vector<double>& getXi() = 0;

    /*!
     * \brief Return the wt vector.
     *
     * Retrieves a vector whose elements contain the weights associated with
     * each direction omega_m. 
     *
     * See comments for getMu().
     */
    virtual const vector<double>& getwt() = 0;

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
    virtual double getMu( const int ) const = 0;

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
    virtual double getEta( const int ) const = 0;

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
    virtual double getXi( const int ) const = 0;

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
    virtual double getwt( const int ) const = 0;

    /*!
     * \brief Returns the Omega vector for all directions.
     *
     * Returns a vector of length numAngles.  Each entry is a length three 
     * vector of doubles that together represent the m-th discrete
     * direction. 
     */
//     virtual const vector< vector<double> >& getOmega() const = 0;

    /*!
     * \brief Returns the Omega vector for a single direction.
     *
     * Returns the Omega direction vector associated with the client
     * specified index.
     *
     * \param m The direction index must be contained in the range
     * (0,numAngles). 
     */
    virtual vector<double> getOmega( const int ) const = 0;

    /*!
     * \brief Returns the number of directions in the current quadrature set.
     */
    virtual int getNumAngles() const = 0;

    /*!
     * \brief Prints a table containing all quadrature directions and weights.
     */
    virtual void display() const = 0;

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

    // Other accessors that may be needed:
    //
    // virtual int getLevels() = 0;
    // virtual int getNumAnglesPerOctant() = 0;
};

//===========================================================================//
/*!
 * \class Q1DGaussLeg
 *
 * \brienf A class to encapsulate a 1D Gauss Legendre Quadrature set.
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

    int snOrder;    // defaults to 4.
    int numAngles;  // == snOrder
    double norm;    // defaults to 2.0

    // Quadrature directions and weights.
    vector<double> mu;
    vector<double> eta; // will be an empty vector for this 1D set.
    vector<double> xi;  // will be an empty vector for this 1D set.
    vector<double> wt;

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
    
    const vector<double>& getMu()  { return mu;  };
    const vector<double>& getEta() { return eta; }; 
    const vector<double>& getXi()  { return xi;  }; 
    const vector<double>& getwt()  { return wt;  };
    //    const vector< vector<double> >& getOmega();
    double getMu ( const int ) const;
    double getEta( const int ) const;
    double getXi ( const int ) const;
    double getwt ( const int m ) const { return wt[m]; };
    vector<double> getOmega( const int ) const;
    int getNumAngles()   const { return numAngles; };
    void display()       const;
    string name()        const { return "1D Gauss Legendre"; };
    int dimensionality() const { return 1; };
    int getSnOrder()     const { return snOrder; };

  private:
    
    // IMPLEMENTATION

    // These functions integrate the direction vector over the unit sphere.
    double iDomega() const;
    double iOmegaDomega() const;
    double iOmegaOmegaDomega() const;
};


//===========================================================================//
/*!
 * \class Q3DLevelSym
 *
 * \brienf A class to encapsulate a 3D Level Symmetric quadrature set.
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

    int snOrder; // Quadrature order, numAngles = (snOrder+2)*snOrder
                 // number of levels = snOrder.
                 // defaults to 4.
    int numAngles; // defaults to 24.
    double norm; // Defaults to 4*PI

    // Quadrature directions and weights.
    vector<double> mu;
    vector<double> wt;
    vector<double> eta;
    vector<double> xi;

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

    const vector<double>& getMu()  { return mu;  };
    const vector<double>& getEta() { return eta; };
    const vector<double>& getXi()  { return xi;  };
    const vector<double>& getwt()  { return wt;  };
    //    const vector< vector<double> >& getOmega() const;
    double getMu ( const int ) const;
    double getEta( const int ) const;
    double getXi ( const int ) const;
    double getwt ( const int m ) const { return wt[m]; };
    vector<double> getOmega( const int ) const;
    int getNumAngles()   const { return numAngles; };
    void display()       const;
    string name()        const { return "3D Level Symmetric"; };
    int dimensionality() const { return 3; };
    int getSnOrder()     const { return snOrder; };

  private:
    
    // IMPLEMENTATION

    // These functions integrate the direction vector over the unit sphere.
    double iDomega() const;
    vector<double> iOmegaDomega() const;
    vector<double> iOmegaOmegaDomega() const;

};


} // end namespace rtt_quadrature

#endif                          // __quadrature_Quadrature_hh__

//---------------------------------------------------------------------------//
//                              end of quadrature/Quadrature.hh
//---------------------------------------------------------------------------//
