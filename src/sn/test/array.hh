//----------------------------------*-C++-*----------------------------------//
// array.hh
// Scott Turner
// 19 February 1998
//---------------------------------------------------------------------------//
// @> Simple array classes.
//---------------------------------------------------------------------------//

#ifndef __sn_test_array_hh__
#define __sn_test_array_hh__

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

    const int i1_perm, i2_perm;
    REAL *data;

    // Declare, but do not define, the assignment operator and copy constructor.
    // This prevents clients from calling them and compilers from generating
    // them. If functionallity is needed in the future, they will need to be
    // defined.

    Array2D &operator=( const Array2D &rhs );
    Array2D( const Array2D &rhs );

  public:

    // use default constructor for Array2D()
    // use the following constructor when parameters are passed

    Array2D( int i1, int i2 ): i1_perm(i1), i2_perm(i2)
    {
      data = new REAL [ i1 * i2 ];

      for ( int i=0 ; i < i1 ; i++ ) {
      for ( int j=0 ; j < i2 ; j++ )
        data[ i*i2 + j ] = 0.0;
      }
    }

    // define a correct destructor, the default must be equivalent
    // to "delete data;", since it does not work?

    ~Array2D()
    {
      delete [] data;
    }

    // The following function can be called to re-initialize an array to any
    // value.

    void Array2D_reinit( const REAL init_value )
    {
      for ( int i=0 ; i < i1_perm ; i++ ) {
      for ( int j=0 ; j < i2_perm ; j++ )
        data[ i*i2_perm + j ] = init_value;
      }
    }

    // overload () so that Fortran array style can be used (non-const)

    REAL &operator ()( int in1, int in2 )
    {
      return data[ in2*i1_perm + in1 ];
    }

    // overload () so that Fortran array style can be used (const)

    const REAL &operator ()( int in1, int in2 ) const
    {
      return data[ in2*i1_perm + in1 ];
    }

};

class Array3D
{

  private:

    const int i1_perm, i2_perm, i3_perm;
    REAL *data;

    // Declare, but do not define, the assignment operator and copy constructor.
    // This prevents clients from calling them and compilers from generating
    // them. If functionallity is needed in the future, they will need to be
    // defined.

    Array3D &operator=( const Array3D &rhs );
    Array3D( const Array3D &rhs );

  public:

    // use default constructor for Array3D()
    // use the following constructor when parameters are passed

    Array3D( int i1, int i2, int i3 ) : i1_perm(i1), i2_perm(i2), i3_perm(i3)
    {
      data = new REAL [ i1 * i2 * i3 ];

      for ( int i=0 ; i < i1 ; i++ ) {
      for ( int j=0 ; j < i2 ; j++ ) {
      for ( int k=0 ; k < i3 ; k++ )
        data[ i*i2*i3 + j*i3 + k ] = 0.0;
      } }
    }

    // define a correct destructor, the default must be equivalent
    // to "delete data;", since it does not work?

    ~Array3D()
    {
      delete [] data;
    }

    // The following function can be called to re-initialize an array to any
    // value.

    void Array3D_reinit( const REAL init_value )
    {
      for ( int i=0 ; i < i1_perm ; i++ ) {
      for ( int j=0 ; j < i2_perm ; j++ ) {
      for ( int k=0 ; k < i3_perm ; k++ )
        data[ i*i2_perm*i3_perm + j*i3_perm + k ] = init_value;
      } }
    }

    // overload () so that Fortran array style can be used (non-const)

    REAL &operator ()( int in1, int in2, int in3 )
    {
      return data[ in3*i1_perm*i2_perm + in2*i1_perm + in1 ];
    }

    // overload () so that Fortran array style can be used (const)

    const REAL &operator ()( int in1, int in2, int in3 ) const
    {
      return data[ in3*i1_perm*i2_perm + in2*i1_perm + in1 ];
    }

};

class Array4D
{

  private:

    const int i1_perm, i2_perm, i3_perm, i4_perm;
    REAL *data;

    // Declare, but do not define, the assignment operator and copy constructor.
    // This prevents clients from calling them and compilers from generating
    // them. If functionallity is needed in the future, they will need to be
    // defined.

    Array4D &operator=( const Array4D &rhs );
    Array4D( const Array4D &rhs );

  public:

    // use default constructor for Array4D()
    // use the following constructor when parameters are passed

    Array4D( int i1, int i2, int i3, int i4 ) : i1_perm(i1), i2_perm(i2),
                                                i3_perm(i3), i4_perm(i4)
    {
      data = new REAL [ i1 * i2 * i3 * i4 ];

      for ( int i=0 ; i < i1 ; i++ ) {
      for ( int j=0 ; j < i2 ; j++ ) {
      for ( int k=0 ; k < i3 ; k++ ) {
      for ( int l=0 ; l < i4 ; l++ )
        data[ i*i2*i3*i4 + j*i3*i4 + k*i4 + l ] = 0.0;
      } } }
    }

    // define a correct destructor, the default must be equivalent
    // to "delete data;", since it does not work?

    ~Array4D()
    {
      delete [] data;
    }

    // The following function can be called to re-initialize an array to any
    // value.

    void Array4D_reinit( const REAL init_value )
    {
      for ( int i=0 ; i < i1_perm ; i++ ) {
      for ( int j=0 ; j < i2_perm ; j++ ) {
      for ( int k=0 ; k < i3_perm ; k++ ) {
      for ( int l=0 ; l < i4_perm ; l++ )
        data[ i*i2_perm*i3_perm*i4_perm + j*i3_perm*i4_perm + k*i4_perm + l ] =
                                                                     init_value;
      } } }
    }

    // overload () so that Fortran array style can be used (non-const)

    REAL &operator ()( int in1, int in2, int in3, int in4 )
    {
      return data[ in4*i1_perm*i2_perm*i3_perm + in3*i1_perm*i2_perm +
                   in2*i1_perm + in1 ];
    }

    // overload () so that Fortran array style can be used (const)

    const REAL &operator ()( int in1, int in2, int in3, int in4 ) const
    {
      return data[ in4*i1_perm*i2_perm*i3_perm + in3*i1_perm*i2_perm +
                   in2*i1_perm + in1 ];
    }

};

#endif                          // __sn_test_array_hh__

//---------------------------------------------------------------------------//
//                              end of sn/test/array.hh
//---------------------------------------------------------------------------//

