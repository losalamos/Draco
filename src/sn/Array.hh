//----------------------------------*-C++-*----------------------------------//
// Array.hh
// Scott Turner
// 17 April 1998
//---------------------------------------------------------------------------//
// @> Simple array classes.
//---------------------------------------------------------------------------//

#ifndef __sn_Array_hh__
#define __sn_Array_hh__

//===========================================================================//
// class Array*D -
// These array classes allow the use of Fortran-style arrays
//   i.e. sigtot(i,j,k) = 0.5
// and memory is also Fortran-style, with (i,0,0) values being
// stored sequentially first, then (i,1,0), etc...
//===========================================================================//

class Array2D
{

    private:

        const int i1_pv, i2_pv;
        REAL *addr;

        // Declare, but do not define, the copy constructor and assignment op.
        // This prevents clients from calling it and compilers from generating
        // it. If functionallity is needed in the future, it will need to be
        // defined.

        Array2D( const Array2D &rhs );
        Array2D &operator=( const Array2D &rhs );

    public:

        // use default constructor for Array2D()

        // use the following constructor when parameters are passed

        Array2D( int i1, int i2 ) : i1_pv(i1), i2_pv(i2)
        {
            addr = new REAL [ i1 * i2 ];

            for ( int i=0 ; i < i2 ; i++ )
                for ( int j=0 ; j < i1 ; j++ )
                    addr[ i*i1 + j ] = 0.0;
        }

        // define a destructor

        ~Array2D()
        {
            delete [] addr;
        }

        // The following function can be called to re-initialize to any value

        void Array2D_reinit( const REAL init_value )
        {
            for ( int i=0 ; i < i2_pv ; i++ )
                for ( int j=0 ; j < i1_pv ; j++ )
                    addr[ i*i1_pv + j ] = init_value;
        }

        // overload () so that Fortran array style can be used (non-const)

        REAL &operator ()( int in1, int in2 )
        {
            return addr[ in2*i1_pv + in1 ];
        }

        // overload () so that Fortran array style can be used (const)

        const REAL &operator ()( int in1, int in2 ) const
        {
            return addr[ in2*i1_pv + in1 ];
        }

};

class Array3D
{

    private:

        const int i1_pv, i2_pv, i3_pv;
        REAL *addr;

        // Declare, but do not define, the copy constructor and assignment op.
        // This prevents clients from calling it and compilers from generating
        // it. If functionallity is needed in the future, it will need to be
        // defined.

        Array3D( const Array3D &rhs );
        Array3D &operator=( const Array3D &rhs );

    public:

        // use default constructor for Array3D()

        // use the following constructor when parameters are passed

        Array3D( int i1, int i2, int i3 ) : i1_pv(i1), i2_pv(i2), i3_pv(i3)
        {
            addr = new REAL [ i1 * i2 * i3 ];

            for ( int i=0 ; i < i3 ; i++ )
                for ( int k=0 ; k < i2 ; k++ )
                    for ( int j=0 ; j < i1 ; j++ )
                        addr[ i*i1*i2 + k*i1 + j ] = 0.0;
        }

        // define a destructor

        ~Array3D()
        {
            delete [] addr;
        }

        // The following function can be called to re-initialize to any value

        void Array3D_reinit( const REAL init_value )
        {
            for ( int i=0 ; i < i3_pv ; i++ )
                for ( int k=0 ; k < i2_pv ; k++ )
                    for ( int j=0 ; j < i1_pv ; j++ )
                        addr[ i*i1_pv*i2_pv + k*i1_pv + j ] = init_value;
        }

        // overload () so that Fortran array style can be used (non-const)

        REAL &operator ()( int in1, int in2, int in3 )
        {
            return addr[ in3*i1_pv*i2_pv + in2*i1_pv + in1 ];
        }

        // overload () so that Fortran array style can be used (const)

        const REAL &operator ()( int in1, int in2, int in3 ) const
        {
            return addr[ in3*i1_pv*i2_pv + in2*i1_pv + in1 ];
        }

};

class Array4D
{

    private:

        const int i1_pv, i2_pv, i3_pv, i4_pv;
        REAL *addr;

        // Declare, but do not define, the copy constructor and assignment op.
        // This prevents clients from calling it and compilers from generating
        // it. If functionallity is needed in the future, it will need to be
        // defined.

        Array4D( const Array4D &rhs );
        Array4D &operator=( const Array4D &rhs );

    public:

        // use default constructor for Array4D()

        // use the following constructor when parameters are passed

        Array4D( int i1, int i2, int i3, int i4 ) : i1_pv(i1), i2_pv(i2),
                                                    i3_pv(i3), i4_pv(i4)
        {
          addr = new REAL [ i1 * i2 * i3 * i4 ];

          for ( int l=0 ; l < i4 ; l++ )
              for ( int i=0 ; i < i3 ; i++ )
                  for ( int k=0 ; k < i2 ; k++ )
                      for ( int j=0 ; j < i1 ; j++ )
                          addr[ l*i1*i2*i3 + i*i1*i2 + k*i1 + j ] = 0.0;
        }

        // define a destructor

        ~Array4D()
        {
            delete [] addr;
        }

        // The following function can be called to re-initialize to any value

        void Array4D_reinit( const REAL init_value )
        {
          for ( int l=0 ; l < i4_pv ; l++ )
              for ( int i=0 ; i < i3_pv ; i++ )
                  for ( int k=0 ; k < i2_pv ; k++ )
                      for ( int j=0 ; j < i1_pv ; j++ )
                          addr[ l*i1_pv*i2_pv*i3_pv +
                                      i*i1_pv*i2_pv +
                                            k*i1_pv + j ] = init_value;
        }

        // overload () so that Fortran array style can be used (non-const)

        REAL &operator ()( int in1, int in2, int in3, int in4 )
        {
            return addr[ in4*i1_pv*i2_pv*i3_pv +
                               in3*i1_pv*i2_pv +
                                     in2*i1_pv + in1 ];
        }

        // overload () so that Fortran array style can be used (const)

        const REAL &operator ()( int in1, int in2, int in3, int in4 ) const
        {
          return addr[ in4*i1_pv*i2_pv*i3_pv +
                             in3*i1_pv*i2_pv +
                                   in2*i1_pv + in1 ];
        }

};

#endif                          // __sn_Array_hh__

//---------------------------------------------------------------------------//
//                              end of Array.hh
//---------------------------------------------------------------------------//

