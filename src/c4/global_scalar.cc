//----------------------------------*-C++-*----------------------------------//
// global_scalar.cc
// Maurice LeBrun
// Wed Feb  1 16:01:58 1995
//---------------------------------------------------------------------------//
// @> Global C4 functions for a scalar architecture
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
// Miscellaneous

void C4_Init( int& argc, char **& argv ) {}

void C4_Finalize() {}

int  C4_node()
{
    return 0;
}

int  C4_nodes()
{
    return 1;
}

int  C4_group()
{
    return 0;
}

void C4_gsync() {}

//---------------------------------------------------------------------------//
// Send/receives

int C4_Send( void *buf, int size, int dest, int tag, int group )
{
    return C4_SUCCESS;
}

int C4_Recv( void *buf, int size, int source, int tag, int group )
{
    return 0;
}

C4_Req C4_SendAsync( void *buf, int size, int dest, int tag, int group )
{
    C4_Req r;
    return r;
}

C4_Req C4_RecvAsync( void *buf, int size, int source, int tag, int group )
{
    C4_Req r;
    return r;
}

void C4_SendAsync( C4_Req& r, void *buf, int size, int dest, int tag, 
		   int group /*=0*/ )
{
}

void C4_RecvAsync( C4_Req& r, void *buf, int size, int source, int tag, 
		   int group /*=0*/ )
{
}

int C4_Send( int *buf, int nels, int dest, int group /*=0*/ )
{
    return C4_SUCCESS;
}

int C4_Recv( int *buf, int nels, int source, int group /*=0*/ )
{
    return 0;
}

int C4_Send( float *buf, int nels, int dest, int group /*=0*/ )
{
    return C4_SUCCESS;
}

int C4_Recv( float *buf, int nels, int source, int group /*=0*/ )
{
    return 0;
}

int C4_Send( double *buf, int nels, int dest, int group /*=0*/ )
{
    return C4_SUCCESS;
}

int C4_Recv( double *buf, int nels, int source, int group /*=0*/ )
{
    return 0;
}

//---------------------------------------------------------------------------//
// Global reductions

void C4_gsum( int& x )                  {}
void C4_gsum( long& x )                 {}
void C4_gsum( float& x )                {}
void C4_gsum( double& x )               {}

void C4_gsum( int *px, int n )          {}
void C4_gsum( long *px, int n )         {}
void C4_gsum( float *px, int n )        {}
void C4_gsum( double *px, int n )       {}

void C4_gmin( int& x )                  {}
void C4_gmin( long& x )                 {}
void C4_gmin( float& x )                {}
void C4_gmin( double& x )               {}

void C4_gmax( int& x )                  {}
void C4_gmax( long& x )                 {}
void C4_gmax( float& x )                {}
void C4_gmax( double& x )               {}

//---------------------------------------------------------------------------//
//                              end of global_scalar.cc
//---------------------------------------------------------------------------//
