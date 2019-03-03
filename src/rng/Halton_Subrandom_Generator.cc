//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   rng/Halton_Subrandom_Generator.cc
 * \author Kent Budge
 * \brief  Define methods of class Halton_Subrandom_Generator
 * \note   Copyright (C) 2016-2019 Los Alamos National Laboratory,
 *         All rights reserved. */
//---------------------------------------------------------------------------//

#include "Halton_Subrandom_Generator.hh"
#include "ds++/Assert.hh"

namespace rtt_rng {

//---------------------------------------------------------------------------//
/*!
 * \param count Dimension of the subrandom vector generated by this object.
 */
Halton_Subrandom_Generator::Halton_Subrandom_Generator(unsigned const count)
    : Subrandom_Generator(count), sequences_() {
  /* empty */
}

//---------------------------------------------------------------------------//
void Halton_Subrandom_Generator::shift_vector() {
  ++count_;
  size_t const N = sequences_.size();
  for (size_t i = 0; i < N; ++i) {
    sequences_[i].shift();
  }
  element_ = 0;
}

//---------------------------------------------------------------------------//
double Halton_Subrandom_Generator::shift() {
  if (element_ == sequences_.size()) {
    Check(sequences_.size() < UINT_MAX);
    sequences_.push_back(
        Halton_Sequence(static_cast<unsigned>(sequences_.size()), count_));
  }
  double const Result = sequences_[element_].lookahead();
  ++element_;
  return Result;
}

} // end namespace rtt_rng

//---------------------------------------------------------------------------//
// end of Halton_Subrandom_Generator.cc
//---------------------------------------------------------------------------//
