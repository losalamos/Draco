/* Minimal main program -- everything is loaded from the library */

extern "C" int Py_Main( int argc, char *argv[] );

int main( int argc, char *argv[] )
{
    return Py_Main(argc, argv);}
