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

//#include <iostream>
#include <vector>

namespace rtt_quadrature
{

using std::vector;

//using std::cout;
//using std::cin;
 
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
 * \example quadrature/test/it_quad.cc
 * 
 * Example of Quadrature usage.  In this example the test code requests input 
 * from the user via stdin, creates a quadrature set, displays some
 * information about the set and then destroys itself.
 *
 * Do I need a creator, destructor?
 */
// revision history:
// -----------------
// 0) original
// 
//===========================================================================//

class Quadrature 
{
  public:

    // ACCESSORS

    /*!
     * \brief Retrieves a vector containing the first (mu) component of the
     * direction vector.
     *
     * The direction vector, Omega, in general has three components 
     * (mu, eta, xi).  
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

    virtual const vector<double> getmu() = 0;

    /*!
     * \brief Retrieves a vector whose elements contain the weights
     * associated with each direction omega_m.
     *
     * See comments for getmu().
     */
    virtual const vector<double> getwt() = 0;

    /*!
     * \brief Returns the number of directions in the current quadrature set.
     */
    virtual int getNumAngles() = 0;

    /*!
     * \brief Prints a table containing all quadrature directions and weights.
     */
    virtual void display() = 0;

};



class Q1DGaussLeg : public Quadrature
{

    // NESTED CLASSES and TYPEDEFS

    // DATA

    int snOrder;
    int numAngles;

    // Quadrature directions and weights.
    vector<double> mu;
    vector<double> wt;

  public:

    // CREATORS
    
    Q1DGaussLeg( int snOrder_ );
    Q1DGaussLeg(const Q1DGaussLeg &rhs);
    ~Q1DGaussLeg();

    // MANIPULATORS
    
    Q1DGaussLeg& operator=(const Q1DGaussLeg &rhs);

    // ACCESSORS

    const vector<double> getmu() { return mu; };
    const vector<double> getwt() { return wt; };
    int getNumAngles() { return numAngles; };
    void display();

  private:
    
    // IMPLEMENTATION

    // Generate a 1D Gauss-Legendre set with specified number of angles.

};

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
    
    Q3DLevelSym( int snOrder_, double norm_);
    Q3DLevelSym(const Q3DLevelSym &rhs);
    ~Q3DLevelSym();

    // MANIPULATORS
    
    Q3DLevelSym& operator=(const Q3DLevelSym &rhs);

    // ACCESSORS

    const vector<double> getmu()  { return mu;  };
    const vector<double> geteta() { return eta; };
    const vector<double> getxi()  { return xi;  };
    const vector<double> getwt()  { return wt;  };
    int getNumAngles() { return numAngles; };
    void display();

  private:
    
    // IMPLEMENTATION

    bool valid_levels( int );
    bool valid_norm( double );

};


} // end namespace rtt_quadrature

#endif                          // __quadrature_Quadrature_hh__

//---------------------------------------------------------------------------//
//                              end of quadrature/Quadrature.hh
//---------------------------------------------------------------------------//
