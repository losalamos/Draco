# Geoffrey Furnish
# 28 October 1996

"""Amorphous Data File eXaminer

Classes:

- ADFX

"""

from Tkinter import *
from FileDialog import *

from gnl import *

#from K_Mesh import *

class ADFX(Frame):

    "Standard key list browser for amorphous data files."

    def __init__( s, master=None ):
	Frame.__init__( s, master )

	s.create_widgets()

    def create_widgets(s):
	"Build the megawidget interface for browsing ADFile's."

	s.title = Label( s, text="HED Amorphous Data File eXaminer" )
	s.title.pack()

	s.file = Label(s)
	s.file.pack( fill=X )

	# Construct the frame which will hold the listbox and
	# scrollbars and the so forth.

	s.f = Frame(s)
	s.f.rowconfigure( 1, weight=1, minsize=0 )
	s.f.columnconfigure( 0, weight=1, minsize=0 )

	s.f.key_list = Label( s.f, text="Key List:" )

	s.f.lb = Listbox( s.f )

	s.f.vb = Scrollbar( s.f, orient=VERTICAL )
	s.f.hb = Scrollbar( s.f, orient=HORIZONTAL )

	s.f.key_list.grid( row=0, columnspan=2, sticky='ew' )
	s.f.lb.grid( row=1, column=0, sticky='news' )
	s.f.vb.grid( row=1, column=1, sticky='ns' )
	s.f.hb.grid( row=2, column=0, sticky='ew' )

	s.f.pack( expand=1, fill='both' )

	# Now let's configure the listbox and scrollbars to work
	# together.

	s.f.lb['yscrollcommand'] = s.f.vb.set
	s.f.lb['xscrollcommand'] = s.f.hb.set

	s.f.vb['command'] = s.f.lb.yview
	s.f.hb['command'] = s.f.lb.xview

	# Now buid the bottom row control buttons.

	s.brow = Frame(s)
	s.brow.open = Button( s.brow, text="Open File", command=s.open )
	s.brow.open.pack( side=LEFT, expand=1, fill=X )
	s.brow.clone = Button( s.brow, text="Clone", command=s.clone )
	s.brow.clone.pack( side=LEFT, expand=1, fill=X )
	s.brow.quit = Button( s.brow, text="Dismiss", command=s.dismiss )
	s.brow.quit.pack( side=LEFT, expand=1, fill=X )
	s.brow.pack(fill=X )

	s.pack( expand=1, fill=BOTH )

	# Useful bindings.

	s.f.lb.bind( '<Double-Button-1>', s.autofire )

    def open(s):
	"Open a new ADFile for browsing."

	fd = FileDialog(s)
	file = fd.go( pattern="*.dat" )
	if file is not None:
	    print "Preparing to open ", file

	    # Should clear the listbox.
	    lb = s.f.lb
	    lb.delete( 0, END )

	    # Now open the file.

	    s.adfile = new_adfile( file, 0, 0 )

	    # Now fetch the keys.

	    keys = s.adfile.Get_keys()

	    # Now insert the keys into the listbox.

	    for key in keys:
		lb.insert( END, key )

	    # If we got this far, things must be cool, so adorn
	    # ourselves.

	    s.file['text'] = file

    def clone(s):
	"Bring up another browser megawidget."

	return ADFX( Toplevel() )

    def dismiss(s):
	"Destroy this browser megawidget."

	s.master.destroy()

    def autofire( s, event ):
	"""Attempt to fire up a viewer for the record with the
	selected key.  This is accomplished by regexp matching against
	the key itself.  Typically the first few characters of the key
	will identify the data type for the record, although more
	sophisticated nomenclatures are certainly possible."""

	key = s.selection_get()

	if key[:5] == "text:":
	    return ADF_Text_Viewer( s.adfile, key, Toplevel() )

	if key[:5] == "Text:":
	    return ADF_Text_Viewer( s.adfile, key, Toplevel() )

##	if key[:10] == "K_M::zcsf:":
##	    return K_M_zcsf( s.adfile, key, Toplevel() )

##	if key[:10] == "K_M::ncsf:":
##	    return K_M_ncsf( s.adfile, key, Toplevel() )

	print 'Unrecognized, trying to autofire on key "%s".' % key

class ADF_Text_Viewer(Frame):

    "A viewer widget for ADFile text records"

    def __init__( s, adfile=None, key=None, master=None ):
	Frame.__init__( s, master )

	if adfile is None:
	    print "You must give me an ADFile!"
	    return

	if key is None:
	    print "You must give me a key!"
	    return

	s.adfile = adfile
	s.key = key

	s.create_widgets()

	# Okay, now we need to get the adfile to find the record with
	# this key, and then we need to read the data out.

	entry = s.adfile.Locate( key, 1 )
	len = s.adfile.len(entry)
	data = s.adfile.binread(len)
	s.txtpane.insert( END, data )

    def create_widgets(s):
	"Build the text record visualizer megawidget."

	s.lab = Label( s, text=s.key )
	s.lab.pack( fill=X )

	s.txtpane = Text( s )
	s.txtpane.pack( expand=1, fill=BOTH )

	s.quit = Button( s, text="Dismiss",
			 command=lambda o=s: o.master.destroy() )
	s.quit.pack( fill=X )

	s.pack()
