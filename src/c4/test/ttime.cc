#include "../global.hh"

#include <iostream>

#define MAX 1000000

#define double float

void f( double *z )
{
    for( int i=0; i < MAX; i++ )
        z[i] = 0.;
}

int main( int argc, char *argv[] )
{
    double x[ MAX ], y[ MAX ], z[ MAX ], u[ MAX ], v[ MAX ];

    for( int i=0; i < MAX; i++ )
    {
        x[i] = 0.;
        y[i] = 2.;
        z[i] = 3.;
        u[i] = 4.;
        v[i] = 5.;
    }

    C4::Init( argc, argv );

    double tr = C4::Wtick();
    double ts = C4::Wtime();

    for( int i=0; i < MAX; i++ )
        x[i] = 4.*y[i] + 5.*z[i] + 6.*u[i] + 7.*v[i];

    double te = C4::Wtime();

    double flop_rate = 7.*MAX / (te - ts);

    std::cout << "It took " << te - ts << " +/- " << tr << " seconds.\n";
    std::cout << "Flop rate was " << flop_rate << "\n";
    std::cout << "which is " << flop_rate * 1.0e-6 << " megaflops.\n";

    f( x );

    if (tr > 0.)
        std::cout << "test passed\n";
    else
        std::cout << " test failed\n";

    C4::Finalize();
}
