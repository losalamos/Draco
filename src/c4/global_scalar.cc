//----------------------------------*-C++-*----------------------------------//
// global_scalar.cc
// Maurice LeBrun
// Wed Feb  1 16:01:58 1995
//---------------------------------------------------------------------------//
// @> Global C4 functions for a scalar architecture
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
// Miscellaneous

C4_NAMESPACE_BEG

void Init( int& argc, char **& argv ) {}

void Finalize() {}

int  node()
{
    return 0;
}

int  nodes()
{
    return 1;
}

int  group()
{
    return 0;
}

void gsync() {}

//---------------------------------------------------------------------------//
// Send/receives
/*
int Send( void *buf, int size, int dest, int tag, int group )
{
    return C4_SUCCESS;
}

int Recv( void *buf, int size, int source, int tag, int group )
{
    return 0;
}
*/
C4_Req SendAsync( void *buf, int size, int dest, int tag, int group )
{
    C4_Req r;
    return r;
}

C4_Req RecvAsync( void *buf, int size, int source, int tag, int group )
{
    C4_Req r;
    return r;
}

void SendAsync( C4_Req& r, void *buf, int size, int dest, int tag, 
		int group /*=0*/ )
{
}

void RecvAsync( C4_Req& r, void *buf, int size, int source, int tag, 
		int group /*=0*/ )
{
}

int Send( int *buf, int nels, int dest, int group /*=0*/ )
{
    return C4_SUCCESS;
}

int Recv( int *buf, int nels, int source, int group /*=0*/ )
{
    return 0;
}

int Send( float *buf, int nels, int dest, int group /*=0*/ )
{
    return C4_SUCCESS;
}

int Recv( float *buf, int nels, int source, int group /*=0*/ )
{
    return 0;
}

int Send( double *buf, int nels, int dest, int group /*=0*/ )
{
    return C4_SUCCESS;
}

int Recv( double *buf, int nels, int source, int group /*=0*/ )
{
    return 0;
}

//---------------------------------------------------------------------------//
// Global reductions

void gsum( int& x )                  {}
void gsum( long& x )                 {}
void gsum( float& x )                {}
void gsum( double& x )               {}

void gsum( int *px, int n )          {}
void gsum( long *px, int n )         {}
void gsum( float *px, int n )        {}
void gsum( double *px, int n )       {}

void gmin( int& x )                  {}
void gmin( long& x )                 {}
void gmin( float& x )                {}
void gmin( double& x )               {}

void gmax( int& x )                  {}
void gmax( long& x )                 {}
void gmax( float& x )                {}
void gmax( double& x )               {}

C4_NAMESPACE_END

//---------------------------------------------------------------------------//
//                              end of global_scalar.cc
//---------------------------------------------------------------------------//
