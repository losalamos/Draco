//----------------------------------*-C++-*----------------------------------//
// c4run.cc
// Geoffrey M. Furnish
// Thu Feb 19 09:54:06 1998
//---------------------------------------------------------------------------//
// @> Program to manage execution of a C4 SHMEM multi-box program.
//---------------------------------------------------------------------------//

#include <iostream>
#include <string>
#include <vector>

#include <unistd.h>
#include <netdb.h>
#include <netinet/in.h>
#include <sys/socket.h>

#include "ds++/Assert.hh"

using namespace std;

int verbose = 0;

int sock_fe;
int fe_port = 5678;
struct sockaddr_in fe_addr;

struct host_entry {
    string name_;
    int npes_;

    host_entry() {}
    host_entry( string& name, int pes ) : name_(name), npes_(pes) {}
    string name() const { return name_; }
    int npes() const { return npes_; }
};

vector<host_entry> hosts;
int npes = 0;

struct servent *sp;

string command;

//---------------------------------------------------------------------------//
// Parse the ocmmand line options for things we recognize.
//---------------------------------------------------------------------------//

void process_cli( int argc, char ** argv )
{
    argc--, argv++;		// Skip over program name.
    while( argc ) {
	string opt = *argv;

	if (opt == "-v") {
	    verbose++;
	    argc--, argv++;
	    continue;
	}

	if (opt == "-h") {
	    Assert( argc > 2 );
	    string hn = argv[1];
	    int npes = atoi(argv[2]);
	    host_entry h( hn, npes );
	    hosts.push_back( h );
	    argc -= 3, argv += 3;
	    if (verbose)
		cout << "host: " << hn << " nodes: " << npes << endl;
	    continue;
	}

    // Other options here.

//     // Finalize--we have no idea what he's saying...

// 	if (verbose) {
// 	    cout << "Unrecognized option: " << opt << endl;
// 	    argc--, argv++;
// 	}

    // Unrecognized option, so it must be the program the user wants to
    // invoke.  Stuff it all into command so it can be rsh'd.

	if (opt == "--")
	    argc--, argv++;

	while( argc )
	{
	    command += *argv;
	    if (argc > 1)
		command += ' ';
	    argc--, argv++;
	}
    }

    if (command == "") throw "Must specify program to run.";

}

//---------------------------------------------------------------------------//
// Launch the server process on the identified host.
//---------------------------------------------------------------------------//

void launch( const host_entry& hent )
{
    if (verbose)
	cout << "Launching on " << hent.name() << endl;

    string cmd;
    cmd += "rsh ";
    cmd += hent.name();
    cmd += ' ';

// Now form the command to execute on the remote host.

// First, determine current directory.
    cmd += "'(cd /home/furnish/devel/draco/src/c4; ";

    cmd += command;
    cmd += " -npes ";
    char buf[40];
    sprintf( buf, "%d", hent.npes() );
    cmd += buf;
    cmd += " )'";

    cout << "Preparing to execute system command: " << cmd << endl;

    system( cmd.data() );

    cout << "Donce executing command." << endl;

//     const char *ah = hent.name().data();
//     char *ahost = new char[ strlen(ah)+1 ];
//     strcpy( ahost, ah );
//     cout << "ahost = " << ahost << endl;
// //    int fd = rexec( &ahost, sp->s_port, "furnish", NULL, "ls", NULL );
//     int fd = rcmd( &ahost, sp->s_port, "furnish", "furnish", "ls", NULL );
//     cout << "Back from rexec, fd=" << fd << endl;
//     close(fd);
}

//---------------------------------------------------------------------------//
// The purpose of this routine is to bind the port that we will be listening
// on for the connection requests of our client processes.
//---------------------------------------------------------------------------//

void bind_fe_port()
{
    sock_fe = socket( AF_INET, SOCK_STREAM, 0 );
    if (!sock_fe) throw "Can't open socket!";

    memset( &fe_addr, 0, sizeof(fe_addr) );
    fe_addr.sin_family = AF_INET;
    fe_addr.sin_addr.s_addr = htonl( INADDR_ANY );

    int nports = 10;
    while( nports ) {
	fe_addr.sin_port = htons( fe_port );

    // If we can successfully bind, jump out of the loop.
	if ( bind( sock_fe, (struct sockaddr *) &fe_addr,
		   sizeof(fe_addr) ) == 0 )
	    break;

    // That port was in use, try another.
	fe_port++;
	nports--;
    }
    Assert( nports );

// Okay, by this point we have bound to port # fe_port.
    if (verbose)
	cout << "Proxy will listen on port " << fe_port << endl;
}

int main( int argc, char *argv[] )
{
    try {
	process_cli( argc, argv );

	bind_fe_port();

	if (hosts.size() == 0)
	{
	// Just running on the same box we're on.  We can just exec.
	}
	else
	{
    // Doing a multibox run, fire 'em up!
	    for( vector<host_entry>::iterator phost = hosts.begin();
		 phost != hosts.end(); phost++ )
		launch( *phost );
	}
    }
    catch( const char *msg )
    {
	cout << "Failed b/c: " << msg << endl;
    }
    catch(...)
    {
	cout << "Unknown failure!\n";
    }
}

//---------------------------------------------------------------------------//
//                              end of c4run.cc
//---------------------------------------------------------------------------//
