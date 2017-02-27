//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   parser/Constant_Expression.cc
 * \author Kent Budge
 * \date   Wed Jul 26 07:53:32 2006
 * \brief  Definition of methods of class Constant_Expression
 * \note   Copyright © 2016-2017 Los Alamos National Security, LLC
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include "Constant_Expression.hh"

namespace {
using namespace std;

//---------------------------------------------------------------------------//
void upper(ostream &out, char const *const label, double const power,
           bool &first, unsigned &icount) {
  if (power > 0.0) {
    if (!first) {
      out << '*';
    }
    first = false;
    if (power == static_cast<unsigned>(power)) {
      out << label;
      for (unsigned i = 1; i < power; ++i) {
        out << '*' << label;
      }
    } else {
      out << "pow(" << label << "," << power << ")";
    }
  } else if (power < 0.0) {
    icount++;
  }
}

//---------------------------------------------------------------------------//
void lower(ostream &out, char const *const label, double const power,
           bool &first) {
  if (power < 0.0) {
    if (!first) {
      out << '*';
    }
    first = false;
    unsigned ipower = static_cast<unsigned>(-power);
    if (power == ipower) {
      out << label;
      for (unsigned i = 1; i < ipower; ++i) {
        out << '*' << label;
      }
    } else {
      out << "pow(" << label << "," << ipower << ")";
    }
  }
}

} // namespace anonymous

namespace rtt_parser {
using namespace rtt_parser;
//---------------------------------------------------------------------------//
void write_c(Unit const &u, ostream &out) {
  unsigned count = u.conv != 1.0;
  double p(0.0);
  if (u.m != 0) {
    count++;
    p = u.m;
  }
  if (u.kg != 0) {
    count++;
    p = u.kg;
  }
  if (u.s != 0) {
    count++;
    p = u.s;
  }
  if (u.A != 0) {
    count++;
    p = u.A;
  }
  if (u.K != 0) {
    count++;
    p = u.K;
  }
  if (u.mol != 0) {
    count++;
    p = u.mol;
  }
  if (u.cd != 0) {
    count++;
    p = u.cd;
  }
  if (u.rad != 0) {
    count++;
    p = u.rad;
  }
  if (u.sr != 0) {
    count++;
    p = u.sr;
  }

  Require(count != 0); // should not come here if dimensionless

  if (count == 1) {
    if (p < 0 && p == static_cast<int>(p)) {
      out << '(';
    }
  }

  bool first = true;
  if (u.conv != 1.0) {
    out << u.conv;
    first = false;
  }
  unsigned icount = 0;
  upper(out, "m", u.m, first, icount);
  upper(out, "kg", u.kg, first, icount);
  upper(out, "s", u.s, first, icount);
  upper(out, "A", u.A, first, icount);
  upper(out, "K", u.K, first, icount);
  upper(out, "mol", u.mol, first, icount);
  upper(out, "cd", u.cd, first, icount);
  upper(out, "rad", u.rad, first, icount);
  upper(out, "sr", u.sr, first, icount);
  if (first) {
    out << '1';
  }
  first = false;
  if (icount > 0) {
    out << '/';
    if (icount > 1) {
      out << '(';
    }
    lower(out, "m", u.m, first);
    lower(out, "kg", u.kg, first);
    lower(out, "s", u.s, first);
    lower(out, "A", u.A, first);
    lower(out, "K", u.K, first);
    lower(out, "mol", u.mol, first);
    lower(out, "cd", u.cd, first);
    lower(out, "rad", u.rad, first);
    lower(out, "sr", u.sr, first);
    if (icount > 1) {
      out << ')';
    }
  }
  if (count == 1) {
    if (p < 0 && p == static_cast<int>(p)) {
      out << ')';
    }
  }
}

} // end namespace rtt_parser

//---------------------------------------------------------------------------//
//              end of parser/Constant_Expression.cc
//---------------------------------------------------------------------------//
