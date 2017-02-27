//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   parser/Constant_Expression.hh
 * \author Kent Budge
 * \date   Wed Jul 26 07:53:32 2006
 * \brief  Definition of class Constant_Expression
 * \note   Copyright (C) 2016-2017 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef parser_Constant_Expression_hh
#define parser_Constant_Expression_hh

#include "utilities.hh"

namespace rtt_parser {
//---------------------------------------------------------------------------//
DLL_PUBLIC_parser void write_c(Unit const &units, ostream &out);

//===========================================================================//
/*!
 * \class Constant_Expression
 * \brief Represents a constant.
 *
 * This is the only concrete type of Expression that is currently available
 * to clients. The need has never arisen to make any other concrete type of
 * Expression directly available.
 *
 * If you need access to another type of concrete Expression, contact the
 * Capsaicin team (preferably Kent G. Budge if he has not been RIFfed yet.)
 */
//===========================================================================//

//---------------------------------------------------------------------------//
class Constant_Expression : public Expression {
public:
  // CREATORS

  /*! Create a dimensionless constant.
     *
     * \param number_of_variables Number of indepedent variables in the
     * expression.
     *
     * \param value Value of the constant
     */
  Constant_Expression(unsigned const number_of_variables, double const value)
      : Expression(number_of_variables, value * dimensionless) {}

  /*! Create a dimensioned constant.
     *
     * \param number_of_variables Number of indepedent variables in the
     * expression.
     *
     * \param value Dimensions and value of the constant. The value is stored
     * in the \c conv member of this argument.
     */
  Constant_Expression(unsigned const number_of_variables, Unit const &value)
      : Expression(number_of_variables, value) {}

  // ACCESSORS

  virtual bool is_constant() const { return true; }

private:
  // IMPLEMENTATION

  virtual double evaluate_(double const *const) const { return units().conv; }

  virtual void write_(Precedence const, vector<string> const &,
                      ostream &out) const {
    if (is_compatible(units(), dimensionless)) {
      out << units().conv;
    } else {
      write_c(units(), out);
    }
  }

  virtual bool is_constant_(unsigned) const { return true; }
};

} // end namespace rtt_parser

#endif // parser_Constant_Expression_hh

//---------------------------------------------------------------------------//
// end of parser/Constant_Expression.hh
//---------------------------------------------------------------------------//
