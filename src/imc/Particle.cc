//----------------------------------*-C++-*----------------------------------//
// Particle.cc
// Thomas M. Evans
// Fri Jan 30 17:04:24 1998
//---------------------------------------------------------------------------//
// @> source file for Particle class
//---------------------------------------------------------------------------//

#include "imctest/Particle.hh"
#include "imctest/Global.hh"
#include <iostream>
#include <cmath>

IMCSPACE

//---------------------------------------------------------------------------//
// constructors
//---------------------------------------------------------------------------//
// default constructor, explicit definition of Particle required to
// guarantee that the Particle is defined by a Mesh;
// vectors handle their own dynamic memory
template<class MT>
Particle<MT>::Particle(const MT &mesh, int seed, double weight_, double energy_)
    :weight(weight_), energy(energy_),
     r(mesh.Coord().getDim(), 0.0), omega(mesh.Coord().getDim(), 0.0),
     cell(0), random(seed)
{
  // the syntax for initialization of the vector is equivalent to:
  // vector<double> r(dimension, 0.0),
  // which means, "make a vector of type double of dimension
  // initialized with all zeroes--the vectors must be
  // initialized by 0.0 because if the arguments match the iterator
  // constructor will be called instead, see Stroustrup pg. 450
}

// no copy constructor is required because vectors can do assignment
// therefore, we use the default copy constructor (memberwise copy)

//---------------------------------------------------------------------------//
// public member functions
//---------------------------------------------------------------------------//
// calculate source from a given point
template<class MT>
void Particle<MT>::source(vector<double> &r_, const MT &mesh)
{
  // initial location
    r = r_;

  // sample initial direction (isotropic)
    mesh.Coord().setOmega(omega, random);

  // find particle cell
    cell = mesh.getCell(r);
}

// transport particle through mesh
template<class MT>
void Particle<MT>::transport(const MT &mesh, const Opacity<MT> &xs)
{
  // explicit calls from standard library namespace
    using std::log;
  // explicit calls from Global namespace
    using Global::pi;
    
  // end transport conditional
    bool alive = true;

  // transport loop, ended when transport = false
    while (alive)
    {
      // dist-to-collision and dist-to-boundary definitions
        double d_collision, d_boundary;
      // cell face definition
        int face = 0;
        
      // sample distance-to-collision
        d_collision = -log(random.ran()) / sigma;

      // get distance-to-boundary and cell face
        d_boundary  = mesh.getDb(r, omega, cell, face);

      // streaming
        if (d_collision <= d_boundary)
        {
            stream(d_collision);
            alive = collide(mesh, xs);
        }
        else
        {
            stream(d_boundary);
            alive = mesh.Layout()(cell, face);
        }
    }
}

//---------------------------------------------------------------------------//
// private member functions
//---------------------------------------------------------------------------//
// stream particle a distance d
template<class MT>
void Particle<MT>::stream(double distance)
{
  // calculate new location
    for (int i = 0; i <= r.size()-1; i++)
        r[i] = r[i] + distance * omega[i];
}

// calculate everything about a collision
template<class MT>
bool Particle<MT>::collide(const MT &mesh, const Opacity<MT> &xs)
{   
    bool status;

  // determine absorption or collision
    if (random.ran() <= xs.Sigma_a(cell) / xs.Sigma(cell))
        status = false;
    else
    {
        status = true;
        
      // calculate theta and phi (isotropic)
        double costheta, phi;
        costheta = 1 - 2 * random.ran();
        phi      = 2 * Global::pi * random.ran();

      // get new direction cosines
        mesh.Coord().calcOmega(costheta, phi, omega);
    }

    return status;
}

//---------------------------------------------------------------------------//
// overloaded operators
//---------------------------------------------------------------------------//
// no overload assignment operator needed because vector can handle
// its own memory management

CSPACE

//---------------------------------------------------------------------------//
//                              end of Particle.cc
//---------------------------------------------------------------------------//
