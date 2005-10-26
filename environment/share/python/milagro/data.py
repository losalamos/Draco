#!/usr/bin/env python

import Utils, array, string, re, sys

class Cycle_Data:
    """Contains one cycle of Milagro output data. An object of class 
    Cycle_Data is constructed with a file object from which the data is read.
    """
    cycle_header = 'RESULTS FOR CYCLE[ ]+(\d+)'
    time_header  = 'Total time :[ ]+(\d+)'
    cell_row     = '^[ ]+(\d+)'
    cycle_footer = 'END OF CYCLE'

    def __init__( self, file_handle ):
        
        # Initialize data arrays
        self.cell = array.array( 'i' )
        self.tMat = array.array( 'd' )
        self.eMat = array.array( 'd' )
        self.tRad = array.array( 'd' )
        self.eRad = array.array( 'd' )
        self.cycle = 0
        self.clock = 0
        self.found_cycle = True;

        self.__read_next_cycle( file_handle )
    

    def __read_next_cycle( self, file_handle ):
        """Read the next cycle from the file object and store the
        results in the internal arrays. If no cycle is found, set
        self.found_cycle to False."""

        file_position = file_handle.tell()

        cycle_position = Utils.ScanTo( file_handle, self.__class__.cycle_header )

        if cycle_position:

            self.cycle = int( cycle_position.group( 1 ) )
            cell_position = Utils.ScanTo( file_handle, self.__class__.cell_row )

            if cell_position:
                line = cell_position.string
                while cell_position:
                    numbers = string.split( line )
                    self.cell.append( int(   numbers[0] ) )
                    self.tMat.append( float( numbers[1] ) )
                    self.eMat.append( float( numbers[2] ) )
                    self.tRad.append( float( numbers[3] ) )
                    self.eRad.append( float( numbers[4] ) )

                    line = file_handle.readline()
                    cell_position = re.search( self.cell_row, line )

                # Scan to end of cycle
                cycle_position = Utils.ScanTo( file_handle, self.__class__.cycle_footer )
                if not cycle_position:
                    file_handle.seek( file_position )
                    raise RunTimeError, "Missing end-of-cycle footer"

            else:
                # Failed to find a cell row
                file_handle.seek( file_position )
                raise RunTimeError, "No cell row information found"

        else:
            # Didn't find a cycle. This happens after the last one, so
            # it's no biggie. We reset the file position because we
            # still need to scan for the execution time.
            file_handle.seek( file_position )
            self.found_cycle = False;


    def gnuplot_output( self, file_out ):
        """ Write the cycle data to the file object in a form suitable
        for reading  with gnuplot."""

        print >> file_out, "# Data for cycle %i" % ( self.cycle, )
        for i in range( len( self.cell ) ):
            print >> file_out, "%i %f %f %f %f" % ( self.cell[i], 
                                                    self.tMat[i], 
                                                    self.eMat[i], 
                                                    self.tRad[i], 
                                                    self.eRad[i] ) 


##---------------------------------------------------------------------------##


##---------------------------------------------------------------------------##
## class Milagro_Data
##---------------------------------------------------------------------------##

class Milagro_Data:

    """Parses a milagro output file. If scan_cycles is true, store the results as a list of
    Cycle_Data objects. Otherwise, just record the execution time."""
    
    # Class data:
    timing_line  = ' \*\* We ran for[ ]+(\S+) seconds'
        

    def __init__( self, input, scan_cycles = True ):
        self.scan_cycles = scan_cycles
        self.input       = input

        self.execution_time = 0
        self.cycles = []

        self.__extract_data()


    def __extract_cycle_data( self, file ):

        # Read first cycle
        cycle_data = Cycle_Data( file )

        while cycle_data.found_cycle is True:

            # Add the cycle data to the list
            self.cycles.append( cycle_data )

            # Get the next cycle
            cycle_data = Cycle_Data( file )



    def __extract_data( self ):

        file = Utils.openFileOrString( self.input )

        if self.scan_cycles:
            self.__extract_cycle_data( file )

        # Scan for the timing data
        time_position = Utils.ScanTo( file, self.timing_line )

        if time_position:
            self.execution_time = float( time_position.group( 1 ) )
        else:
            raise RunTimeError, "Could not find timing data"
        


    def gnuplot_output( self, file_out ):

        """Write all of the cycle data in a form suitable for plotting
        with gnuplot. Cycles are stored as sperate "data blocks" in
        the file."""

        for cycle in self.cycles:
            print "Printing data for cycle %i" % cycle.cycle
            cycle.gnuplot_output( file_out )
            print >> file_out

        # Extra line to seperate data sets
        print >> file_out

##---------------------------------------------------------------------------##
## import_data
##---------------------------------------------------------------------------##

def import_file( filename, scan_cycles = True ):

    """Create and return a Milagro_Data object based on the given
    file"""

    return Milagro_Data( filename, scan_cycles )


##---------------------------------------------------------------------------##
## gnuplot
##---------------------------------------------------------------------------##
def gnuplot_output( data_name, filename ):

    """Create a gnuplot file from a Milagro_Data object, Cycle_Data
    object or a milagro output file"""

    # If we can't open it as a filename, assume it is already a file:
    try:
        file = open( filename, 'w' )
    except TypeError:
        file = filename

    # Try calling gnuplot_output method of the data argument
    try:
        data_name.gnuplot_output( file )
        return
    except:
        # Try importing data, treating data_name as a filename
        try:
            data = import_file( data_name )
            gnuplot_output( data, filename )
            return
        except:
            raise 'ERROR: Could not process data object'


if __name__=="__main__":
    import apihelper
    classes = [Milagro_Data, Cycle_Data]
    for item in classes:
        print item.__name__, ":"
        print item.__doc__
        print "Methods:"
        apihelper.help( item, 10, False )
    functions = [import_file, gnuplot_output]
    for item in functions:
        print item.__name__, ":"
        print item.__doc__





