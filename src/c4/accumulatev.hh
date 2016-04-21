//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   c4/accumulatev.hh
 * \author Kelly Thompson, Thomas M. Evans, Bob Webster
 * \date   Monday, Nov 05, 2012, 13:42 pm
 * \brief  Data accumulatev functions
 * \note   Copyright (C) 2016 Los Alamos National Security, LLC.
 *         All rights reserved.
 *
 * This file contains the declarations for determinate and indeterminate
 * variable length accumulate functions.
 */
//---------------------------------------------------------------------------//
// $Id: accumulatev.hh 6288 2011-12-04 03:43:52Z kellyt $
//---------------------------------------------------------------------------//

#ifndef c4_accumulatev_hh
#define c4_accumulatev_hh

#include <vector>

namespace rtt_c4
{

//---------------------------------------------------------------------------// 
/*!
 * \brief Accumulates an array across all processors using an arbitrary
 * functor.
 *
 * Accumulate is a generalization of summation: it computes the sum (or some
 * other binary operation) of init and all of the elements in the across the
 * processor for each container position..
 *
 * The container "local" given by the iterators [localBegin, localEnd) must be
 * the same size and type on all processors.  This function then computes the
 * accumulation row by row of the value of that row, accumulated across
 * processors.
 *
 * Otherwise.  this function is very similar to the STL function
 * std::accumulate.  with two important differences in the requirements on the
 * binary function and the application of the "initial" value init.
 *
 * In the STL, the function object binary_op is not required to be either
 * commutative or associative: the order of all of accumulate's operations is
 * specified. The result is first initialized to init. Then, for each iterator
 * i in [first, last), in order from beginning to end, it is updated by result
 * = result + *i (in the first version) or result = binary_op(result, *i) (in
 * the second version).
 *    
 * This parallel version first performs the accumulation on the individual
 * processors and then performs an accumulation and the sends the result up a
 * binary tree for further accumulation. This can have an unintended
 * consequence of repeating the inclusion of the init value. So to be safe,
 * arithmetic (+,-) operations should be provided with an init of 0,
 * arithmetic (*,/) operations should be provided with an init of 1, max/min
 * functions with appropriate lower/upper bounds. It also leads to the
 * requirement that the binary function be associative and commutative..
 *
 * HISTORY
 * 
 * This implemenation was developed by Bob Webster (Tigs communication
 * library) and then adapted for use in Jayenne by Tom Evan (wedgehog_gs).
 * Kelly Thompson updated and moved the routine from wedgehog_gs to c4 in
 * Nov., 2012.
 
 *     \param localBegin itertor pointing to beginning of the data
 *     \param localEnd itertor pointing to beginning of the data
 *     \param init initial value for the reduction
 *     \param op the binary operation to apply.
 *   
 *     \return The reduced value.
 */
template<typename T, typename Tciter, typename BinaryOp>
void accumulatev( Tciter   localBegin, Tciter   localEnd,
                  T        init,       BinaryOp op );

} // end namespace rtt_c4

#include "accumulatev.i.hh"

#endif // c4_accumulatev_hh

//---------------------------------------------------------------------------//
// end of c4/accumulatev.hh
//---------------------------------------------------------------------------//
