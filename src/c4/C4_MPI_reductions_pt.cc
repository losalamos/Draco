//----------------------------------*-C++-*----------------------------------//
/*!
 * \file   c4/C4_MPI_reductions_pt.cc
 * \author Thomas M. Evans
 * \date   Mon Mar 25 11:12:35 2002
 * \brief  C4 MPI global reduction instantiations.
 */
//---------------------------------------------------------------------------//
// $Id$
//---------------------------------------------------------------------------//

#include <c4/config.h>

#ifdef C4_MPI

#include "C4_MPI.t.hh"

namespace rtt_c4
{

//---------------------------------------------------------------------------//
// EXPLICIT INSTANTIATIONS OF GLOBAL REDUCTIONS
//---------------------------------------------------------------------------//

template void global_sum< short          >( short          &);
template void global_sum< unsigned short >( unsigned short &);
template void global_sum< int            >( int            &);
template void global_sum< unsigned int   >( unsigned int   &);
template void global_sum< long           >( long           &);
template void global_sum< unsigned long  >( unsigned long  &);
template void global_sum< float          >( float          &);
template void global_sum< double         >( double         &);
template void global_sum< long double    >( long double    &);
template void global_sum< unsigned long long >( unsigned long long &);


template void global_prod< short          >( short          &);
template void global_prod< unsigned short >( unsigned short &);
template void global_prod< int            >( int            &);
template void global_prod< unsigned int   >( unsigned int   &);
template void global_prod< long           >( long           &);
template void global_prod< unsigned long  >( unsigned long  &);
template void global_prod< float          >( float          &);
template void global_prod< double         >( double         &);
template void global_prod< long double    >( long double    &);
template void global_prod< unsigned long long >( unsigned long long &);


template void global_max< short          >( short          &);
template void global_max< unsigned short >( unsigned short &);
template void global_max< int            >( int            &);
template void global_max< unsigned int   >( unsigned int   &);
template void global_max< long           >( long           &);
template void global_max< unsigned long  >( unsigned long  &);
template void global_max< float          >( float          &);
template void global_max< double         >( double         &);
template void global_max< long double    >( long double    &);
template void global_max< unsigned long long >( unsigned long long &);


template void global_min< short          >( short          &);
template void global_min< unsigned short >( unsigned short &);
template void global_min< int            >( int            &);
template void global_min< unsigned int   >( unsigned int   &);
template void global_min< long           >( long           &);
template void global_min< unsigned long  >( unsigned long  &);
template void global_min< float          >( float          &);
template void global_min< double         >( double         &);
template void global_min< long double    >( long double    &);
template void global_min< unsigned long long >( unsigned long long &);


template void global_sum< short          >( short          *, int);
template void global_sum< unsigned short >( unsigned short *, int);
template void global_sum< int            >( int            *, int);
template void global_sum< unsigned int   >( unsigned int   *, int);
template void global_sum< long           >( long           *, int);
template void global_sum< unsigned long  >( unsigned long  *, int);
template void global_sum< float          >( float          *, int);
template void global_sum< double         >( double         *, int);
template void global_sum< long double    >( long double    *, int);

template void global_prod< short          >( short          *, int);
template void global_prod< unsigned short >( unsigned short *, int);
template void global_prod< int            >( int            *, int);
template void global_prod< unsigned int   >( unsigned int   *, int);
template void global_prod< long           >( long           *, int);
template void global_prod< unsigned long  >( unsigned long  *, int);
template void global_prod< float          >( float          *, int);
template void global_prod< double         >( double         *, int);
template void global_prod< long double    >( long double    *, int);
template void global_prod< unsigned long long >( unsigned long long *, int);


template void global_max< short          >( short          *, int);
template void global_max< unsigned short >( unsigned short *, int);
template void global_max< int            >( int            *, int);
template void global_max< unsigned int   >( unsigned int   *, int);
template void global_max< long           >( long           *, int);
template void global_max< unsigned long  >( unsigned long  *, int);
template void global_max< float          >( float          *, int);
template void global_max< double         >( double         *, int);
template void global_max< long double    >( long double    *, int);
template void global_max< unsigned long long >( unsigned long long *, int);


template void global_min< short          >( short          *, int);
template void global_min< unsigned short >( unsigned short *, int);
template void global_min< int            >( int            *, int);
template void global_min< unsigned int   >( unsigned int   *, int);
template void global_min< long           >( long           *, int);
template void global_min< unsigned long  >( unsigned long  *, int);
template void global_min< float          >( float          *, int);
template void global_min< double         >( double         *, int);
template void global_min< long double    >( long double    *, int);
template void global_min< unsigned long long >( unsigned long long *, int);


} // end namespace rtt_c4

#endif // C4_MPI

//---------------------------------------------------------------------------//
//                              end of C4_MPI_reductions_pt.cc
//---------------------------------------------------------------------------//
