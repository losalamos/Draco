//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   meshTest/compose.hh
 * \author Randy M. Roberts
 * \date   Tue Dec 21 09:26:22 1999
 * \brief  A replacement for the non-standard C++ std::compose suit
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __meshTest_compose_hh__
#define __meshTest_compose_hh__

#include <functional>

namespace rtt_meshTest
{
 
template <class Operation1, class Operation2>
class unary_compose 
  : public std::unary_function<typename Operation2::argument_type,
                               typename Operation1::result_type> {
protected:
    Operation1 op1;
    Operation2 op2;
public:
    typedef typename std::unary_function<typename Operation2::argument_type,
            typename Operation1::result_type>::argument_type argument_type;
    unary_compose (const Operation1& x, const Operation2& y)
    : op1 (x), op2 (y) {}

    typename Operation1::result_type 
    operator () (const argument_type& x) const
    {
        return op1 (op2 (x));
    }
};

template <class Operation1, class Operation2>
inline
unary_compose<Operation1, Operation2> compose1 (const Operation1& op1, 
                                               const Operation2& op2)
{
    return unary_compose<Operation1, Operation2> (op1, op2);
}

template <class Operation1, class Operation2, class Operation3>
class binary_compose 
  : public std::unary_function<typename Operation2::argument_type,
                               typename Operation1::result_type> {
protected:
    Operation1 op1;
    Operation2 op2;
    Operation3 op3;
public:
    typedef typename std::unary_function<typename Operation2::argument_type,
            typename Operation1::result_type>::argument_type argument_type;
    binary_compose (const Operation1& x, const Operation2& y, 
                   const Operation3& z) : op1 (x), op2 (y), op3 (z) { }
    typename Operation1::result_type
    operator () (const argument_type& x) const
    {
        return op1 (op2 (x), op3 (x));
    }
};

template <class Operation1, class Operation2, class Operation3>
inline
binary_compose<Operation1, Operation2, Operation3> 
compose2 (const Operation1& op1, const Operation2& op2,
         const Operation3& op3)
{
    return binary_compose<Operation1, Operation2,
                          Operation3> (op1, op2, op3);
}

} // end namespace rtt_meshTest

#endif                          // __meshTest_compose_hh__

//---------------------------------------------------------------------------//
//                              end of meshTest/compose.hh
//---------------------------------------------------------------------------//
