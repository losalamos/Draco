/*-----------------------------------*-C++-*---------------------------------*/
/* MatrixFactoryTraits.hh */
/* Randy M. Roberts */
/* Thu May 27 13:48:53 1999 */
/*---------------------------------------------------------------------------*/
/* @> */
/*---------------------------------------------------------------------------*/

#ifndef __MatrixFactory_MatrixFactoryTraits_hh__
#define __MatrixFactory_MatrixFactoryTraits_hh__

namespace rtt_MatrixFactory
{

template<class Matrix>
struct MatrixFactoryTraits
{
    template<class T>
    static Matrix *create(const T &rep)
    {
	// You should be specializing this class.
	// BogusMethod is being used to trigger a compilation
	// error.
	
	return MatrixFactoryTraits<Matrix>::BogusMethod(rep);
    }
};

} // end namespace rtt_MatrixFactory

#endif    /* __MatrixFactory_MatrixFactoryTraits_hh__ */

/*---------------------------------------------------------------------------*/
/*    end of MatrixFactoryTraits.hh */
/*---------------------------------------------------------------------------*/
