//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   mc/ThreeVector.hh
 * \author H. Grady Hughes
 * \date   Tue Jan 25 15:18:13 MST 2000
 * \brief  Cartesian 3-D vector class header file.
 */
//---------------------------------------------------------------------------//

#ifndef __mc_ThreeVector_hh__
#define __mc_ThreeVector_hh__

#include "ds++/Assert.hh"
#include "rng/Sprng.hh"
#include <vector>
#include <cmath>
#include <iostream>

namespace rtt_mc
{

using rtt_rng::Sprng;

//___________________________________________________________________________//
/*!
 * \class ThreeVector
 * \brief Cartesian 3-D vector class with cross product, dot product, etc.
 */
// revision history:
// -----------------
//  0) original : Committed 2000-01-27.
//___________________________________________________________________________//

class ThreeVector
{

 private:

    //! Scalar field of doubles, used for returning an ordinary vector.
    typedef std::vector<double> SF_DOUBLE;

    //! Cartesian coordinates of the vector.
    double x, y, z;

 public:

    /*!
     * \brief Explicit default constructor.
     *
     * For the moment, no invisible type conversion is allowed.
     *
     * Compiler will generate default copy constructor and assignment operator.
     */
    explicit ThreeVector(double x_ = 0.0, double y_ = 0.0, double z_ = 0.0)
        : x(x_), y(y_), z(z_) {}

    /*!
     * \brief Explicit constructor from ordinary vector<double>.
     *
     * Again, no invisible type conversion is allowed.
     */
    explicit ThreeVector(const SF_DOUBLE &v)
    { Require ( v.size() == 3 ); x = v[0]; y = v[1]; z = v[2]; }

    //! \brief Conversion function for returning an ordinary vector.
    SF_DOUBLE convert() const
    {
        SF_DOUBLE v;
        v.push_back(x); v.push_back(y); v.push_back(z);
        return v;
    }

    // Accessor functions.
    double get_x() const { return x; }
    double get_y() const { return y; }
    double get_z() const { return z; }

    /*!
     * \brief Cross-product of two vectors.
     *
     * V.cross(W) == V x W in the usual notation.
     *
     * Order is important:  V.cross(W) == -W.cross(V)
     */
    const ThreeVector cross(const ThreeVector &r) const
    { return ThreeVector(y*r.z - z*r.y, z*r.x - x*r.z, x*r.y - y*r.x); }

    /*!
     * \brief Dot-product of two vectors.
     *
     * V.dot(W) == V . W in the usual notation.
     *
     * Order is not important:  V.dot(W) == W.dot(V)
     */
    double dot(const ThreeVector &r) const { return x*r.x + y*r.y + z*r.z; }

    //! Calculate the magnitude of a vector.
    double get_norm() const { return std::sqrt(x*x + y*y + z*z); }

    //! \brief Normalize a vector. Non-const: alters the object's data.
    void normalize() { double s = get_norm(); x /= s; y /= s; z /= s; }

    //! \brief Add two vectors.
    const ThreeVector operator+(const ThreeVector &r) const
    { return ThreeVector(x + r.x, y + r.y, z + r.z); }

    //! \brief Subtract one vector from another.
    const ThreeVector operator-(const ThreeVector &r) const
    { return ThreeVector(x - r.x, y - r.y, z - r.z); }

    //! \brief Overloaded equality operator.
    bool operator==(const ThreeVector &rhs) const
    { return (x==rhs.x) && (y==rhs.y) && (z==rhs.z); }

    //! \brief Overloaded inequality operator.
    bool operator!=(const ThreeVector &rhs) const { return !(*this == rhs); }

};  // end class ThreeVector

//    //_ _brief Average of two vectors.
//    inline
//    const ThreeVector center(const ThreeVector &A, const ThreeVector &B)
//    {
//        return ThreeVector(
//            0.5*(A.get_x() + B.get_x()),
//            0.5*(A.get_y() + B.get_y()),
//            0.5*(A.get_z() + B.get_z()));
//    }
//
//    //_ _brief Average of three vectors.
//    inline
//    const ThreeVector center(const ThreeVector &A, const ThreeVector &B,
//        const ThreeVector &C)
//    {
//        return ThreeVector(
//            (A.get_x() + B.get_x() + C.get_x())/3.0,
//            (A.get_y() + B.get_y() + C.get_y())/3.0,
//            (A.get_z() + B.get_z() + C.get_z())/3.0);
//    }
//
//    //_ _brief Average of four vectors.
//    inline
//    const ThreeVector center(const ThreeVector &A, const ThreeVector &B,
//        const ThreeVector &C, const ThreeVector &D)
//    {
//        return ThreeVector(
//            0.25*(A.get_x() + B.get_x() + C.get_x() + D.get_x()),
//            0.25*(A.get_y() + B.get_y() + C.get_y() + D.get_y()),
//            0.25*(A.get_z() + B.get_z() + C.get_z() + D.get_z()));
//    }

    /*!
     * \brief Linear combination of two vectors.
     * \param A        The starting vector.
     * \param B        The ending vector.
     * \param fraction The fractional distance from A to B.
     * \return         fraction*A + (1.0 - fraction)*B.
     *
     * The fraction need not be bounded by (0, 1).  Extrapolation is allowed.
     */
    inline
    const ThreeVector lin_comb(const ThreeVector &A, const ThreeVector &B,
        double fraction)
    {
        return ThreeVector(
            fraction*A.get_x() + (1.0 - fraction)*B.get_x(),
            fraction*A.get_y() + (1.0 - fraction)*B.get_y(),
            fraction*A.get_z() + (1.0 - fraction)*B.get_z());
    }

    //! \brief Sample uniformly within a triangle specified by its vertices.
    inline
    const ThreeVector sample_in_triangle(const ThreeVector &A,
        const ThreeVector &B, const ThreeVector &C, Sprng &random)
    {
        double b, c;

        do
            {
                b = random.ran(); Check ( b > 0.0 && b < 1.0 );
                c = random.ran(); Check ( c > 0.0 && c < 1.0 );
            }
            while (b + c == 1.0);

        if (b + c > 1.0)
            { double temp = b; b = 1.0 - c; c = 1.0 - temp; }

        double a = 1.0 - (b + c);
        return ThreeVector(
            a*A.get_x() + b*B.get_x() + c*C.get_x(),
            a*A.get_y() + b*B.get_y() + c*C.get_y(),
            a*A.get_z() + b*B.get_z() + c*C.get_z());
    }

}   // end namespace rtt_mc

#endif  // __mc_ThreeVector_hh__

//---------------------------------------------------------------------------//
//                              end of mc/ThreeVector.hh
//---------------------------------------------------------------------------//
