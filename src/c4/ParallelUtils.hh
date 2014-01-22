//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   c4/ParallelUtils.hh
 * \author Randy M. Roberts
 * \date   Tue Feb 22 16:38:51 2000
 * \brief  
 * \note   Copyright (C) 2000-2014 Los Alamos National Security, LLC.
 *         All rights reserved.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#ifndef __c4_ParallelUtils_hh__
#define __c4_ParallelUtils_hh__

#include "C4_Functions.hh"

namespace C4
{

/*!
 * Distribute the range [first1, last1) from processor 0,
 * to the corresponding range [first2, first2+size2) on all of the
 * other processors, such that...
 *    processor 0 gets the first chunk of
 *    [first1, last1) into its [first2, first2+size2),
 *    processor 1 gets the next chunk of
 *    [first1, last1) into its [first2, first2+size2), etc.
 */

template<class ForwardIterator1, class OutputIterator2, class Distance2>
void distribute_n(ForwardIterator1 first1, ForwardIterator1 last1,
		  OutputIterator2 first2, Distance2 size2_in,
		  std::forward_iterator_tag, std::output_iterator_tag)
{
    typedef typename std::iterator_traits<ForwardIterator1>::value_type value_type1;

    int size2 = size2_in;
    int totalSize2 = size2;
    C4::gsum(totalSize2);

    int sizesOK = 0;

    if (C4::node() == 0)
	sizesOK = (totalSize2 == std::distance(first1, last1));

    C4::gsum(sizesOK);

    Assert(sizesOK == 1);

    // Processor 0, marches through the other processors
    // finding the size of their output ranges,
    // and sending the next chunk of input, of that size,
    // to the processor.
    
    if (C4::node() == 0)
    {
	// First take care of doing proc 0's **own**
	// distribution of its input into its output.
	
	ForwardIterator1 tmp1 = first1;
	std::advance(tmp1, size2);

	// Copy the chunk of data from the input to the output range.
	
	std::copy(first1, tmp1, first2);
	first1 = tmp1;
	
	for (int i=1; i<C4::nodes(); ++i)
	{
	    // Get the size of proc i's output range.
	    
	    C4::Recv<int>(size2, i);
	
	    value_type1 *buf = new value_type1[size2];

	    // Send the next size2 chunk of data.
	    
	    std::advance(tmp1, size2);
	    std::copy(first1, tmp1, buf);
	    C4::Send<value_type1>(buf, size2, i);

	    // Advance where the next chunk of data will come from.
	    
	    first1 = tmp1;

	    delete [] buf;
	}

	Assert(tmp1 == last1);
    }
    else
    {
	// Send proc 0, this proc's output range size.
	
	C4::Send<int>(size2, 0);

	// Receive the correct size of data from proc 0.
	
	value_type1 *buf = new value_type1[size2];
	C4::Recv<value_type1>(buf, size2, 0);

	// Copy the data from the buffer to the output range.
	
	std::copy(buf, buf+size2, first2);

	delete [] buf;
    }
}

/*!
 * Distribute the range [first1, last1) from processor 0,
 * to the corresponding range [first2, first2+size2) on all of the
 * other processors, such that...
 *    processor 0 gets the first chunk of
 *    [first1, last1) into its [first2, first2+size2),
 *    processor 1 gets the next chunk of
 *    [first1, last1) into its [first2, first2+size2), etc.
 */

template<class ForwardIterator1, class OutputIterator2, class Distance2>
inline void distribute_n(ForwardIterator1 first1, ForwardIterator1 last1,
			 OutputIterator2 first2, Distance2 size2)
{
    distribute_n(first1, last1, first2, size2, 
		 std::iterator_traits<ForwardIterator1>::iterator_category(),
		 std::iterator_traits<OutputIterator2>::iterator_category());
}


/*!
 * Distribute the range [first1, last1) from processor 0,
 * to the corresponding range [first2, last2) on all of the
 * other processors, such that...
 *    processor 0 gets the first chunk of
 *    [first1, last1) into its [first2, last2),
 *    processor 1 gets the next chunk of
 *    [first1, last1) into its [first2, last2), etc.
 */

template<class ForwardIterator1, class ForwardIterator2>
inline void distribute(ForwardIterator1 first1, ForwardIterator1 last1,
		       ForwardIterator2 first2, ForwardIterator2 last2,
		       std::forward_iterator_tag, std::forward_iterator_tag)
{
    distribute_n(first1, last1, first2, std::distance(first2, last2));
}

/*!
 * Distribute the range [first1, last1) from processor 0,
 * to the corresponding range [first2, last2) on all of the
 * other processors, such that...
 *    processor 0 gets the first chunk of
 *    [first1, last1) into its [first2, last2),
 *    processor 1 gets the next chunk of
 *    [first1, last1) into its [first2, last2), etc.
 */

template<class ForwardIterator1, class ForwardIterator2>
inline void distribute(ForwardIterator1 first1, ForwardIterator1 last1,
		       ForwardIterator2 first2, ForwardIterator2 last2)
{
    distribute(first1, last1, first2, last2,
	       std::iterator_traits<ForwardIterator1>::iterator_category(),
	       std::iterator_traits<ForwardIterator2>::iterator_category());
}

/*!
 * Collate the range [first1, last1) from all of the processors, in order,
 * into [result, ...) on processor 0.
 */

template<class InputIterator, class OutputIterator>
void collate(InputIterator first, InputIterator last, OutputIterator result,
	     std::input_iterator_tag, std::output_iterator_tag)
{
    typedef typename std::iterator_traits<InputIterator>::value_type value_type;
    // typedef std::iterator_traits<InputIterator>::difference_type diff_type;
    typedef int diff_type;

    // Only Proc 0 copies any data into the result iterator.
    
    if (C4::node() == 0)
    {
	result = std::copy(first, last, result);
	
	for (int i=1; i<C4::nodes(); ++i)
	{
	    diff_type size;
	    C4::Recv<diff_type>(size, i);
	    value_type *buf = new value_type[size];
	    C4::Recv<value_type>(buf, size, i);
	    
	    result = std::copy(buf, buf+size, result);

	    delete [] buf;
	}
    }
    else
    {
	diff_type size = std::distance(first, last);
	C4::Send<diff_type>(size, 0);
	value_type *buf = new value_type[size];
	std::copy(first, last, buf);
	C4::Send<value_type>(buf, size, 0);

	delete [] buf;
    }
}

/*!
 * Collate the range [first1, last1) from all of the processors, in order,
 * into [result, ...) on processor 0.
 */

template<class InputIterator, class OutputIterator>
inline void collate(InputIterator first, InputIterator last,
		    OutputIterator result)
{
    collate(first, last, result,
	    std::iterator_traits<InputIterator>::iterator_category(),
	    std::iterator_traits<OutputIterator>::iterator_category());
}

}

} // end namespace C4

#endif // __c4_ParallelUtils_hh__

//---------------------------------------------------------------------------//
// end of c4/ParallelUtils.hh
//---------------------------------------------------------------------------//
