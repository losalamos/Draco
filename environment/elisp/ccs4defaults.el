;; ========================================
;; CCS-4 default settings
;; Kelly Thompson
;; 10 June 2002
;; ========================================

;;========================================
;; HOW TO USE:
;;
;; In the file ${HOME}/.xemacs/init.el add one of the following lines
;; (use the line that has the correct path for your machine)
;;
;; CCS2/4 Linux:
;; (load-library "/codes/radtran/vendors/environment/elisp/ccs4defaults.el")
;;
;; ICN Machines: 
;; (load-library "/usr/projects/draco/environment/elisp/ccs4defaults.el")
;;
;; Use personal checkout of Draco:
;; (load-library "/path/to/draco/environment/elisp/ccs4defaults.el")
;;
;; NOTE: ALL defvar values may be overriden in users's init.el.
;;========================================


;; Defaults for personal directories
;; ========================================

(defvar my-home-dir (concat (getenv "HOME") "/"))
(defvar ccs4-env-dirs
      (list (concat my-home-dir "draco/environment/")
	    (concat my-home-dir ".xemacs/")
	    "/usr/projects/draco/environment/"
	    "/codes/radtran/vendors/environment/" )
"\nList of directories that will be prepended to the load-path
if they exist.  The directories <dir>/elisp will also be 
examined and prepended to the the load-path if they exist.
\nAdd to this list by using the following command in personal 
elisp files:
\n\t(setq ccs4-env-dirs (cons \"/path/to/extra/dir/\"))")

;; Update load path using ccs4-env-dirs. Also, set ccs4-elisp-dir
(let ((dirlist ccs4-env-dirs))
  (progn
    ;; process least-imporant directories first.  This ensures that
    ;; the first directory in the list "ccs4-env-dirs" become the
    ;; first directory in the load-path.
    (setq dirlist (reverse dirlist))
    ;; Loop over all directories in the list provided.  Add valid
    ;; elisp directories to the load path (prepend).
    (while dirlist
      (setq ldir (car dirlist))
      (if (file-accessible-directory-p ldir)
	  (setq load-path (cons ldir load-path)))
      (if (file-accessible-directory-p (concat ldir "elisp/"))
	  (setq load-path (cons (concat ldir "elisp/") load-path)))
      ;; If we find the ccs4defauls.el file in this directory then set
      ;; this directory to be the ccs4-elisp-dir
      (if (file-readable-p (concat ldir "elisp/ccs4defaults.el"))
	  (progn	
	    (setq ccs4-env-dir ldir)
	    (setq ccs4-elisp-dir (concat ldir "elisp/"))
	    (setq ccs4-templates-dir (concat ldir "templates/"))
	    ))
      ;; Set extra Info paths:
      (if (file-accessible-directory-p (concat ldir "info/"))
	  (setq Info-directory-list (cons (concat ldir "info/")
					  Info-directory-list )))
		
      ;; remove ldir from dirlist and continue the while-loop.
      (setq dirlist (cdr-safe dirlist)))
    ))

;; Default items for setup
;; To override these defaults add code to ~/.xemacs/init.el that has
;; the form "(setq want-feature nil)".  Values can be nil or t
;; (non-nil). 
;; ========================================
; Draco configuration (nil is the default for all options)
(defvar want-pooma-style-by-default  nil)

(defvar want-rtt-menu                t
"\nIf non-nil, an RTT menu will be placed in the XEmacs menubar.
This menu includes selections to create new files from RTT templates
and selections for frequently used XEmacs features (CVS, compile,
etc.)" )

; rtt-config-modes (nil is the default for all options)

(defvar config-pkg-verobse   nil
"\nIf non-nil, XEmacs will be verbose as it applies
RTT customizations to XEmacs modes.")

(defvar config-pkg-colorize-modeline t
"\nIf non-nil, the modeline background will change color
based on the current editing mode.  For example, an elisp 
buffer will have an orange modeline but a C++ buffer would
have a blue modeline.")
(defvar want-mppl-mode       nil)
(defvar want-tcl-mode        nil)
(defvar want-python-mode     t)
(defvar want-makefile-mode   t)
(defvar want-cc-mode         t)
(defvar want-auctex-pkg      t)
(defvar want-f90-mode        t)
(defvar want-python-mode     t)
(defvar want-change-log-mode t)
(defvar want-fortran-mode    t)
(defvar want-emacs-lisp-mode t)
(defvar want-shell-mode      t)
(defvar want-sgml-mode       t)
(defvar want-autoconf-mode   t)
(defvar want-sh-mode         t)
(defvar want-cvs-mode        t
"\nApply RTT customizations to cvs-mode:\n
   - Loads cvs-mode during XEmacs initialization.
   - Turns on auto-fill-mode for cvs-edit-mode.
   - Turns on modeline coloring if config-pkg-colorize-modeline is
     non-nil.
   - sets cvs-inhibit-copyright-message to t.
   - sets cvs-erase-input-buffer to nil.
   - sets cvs-auto-remove-handled to t.
   - sets cvs-auto-remove-handled-directories to t.")
(defvar want-doxymacs-mode t
"\nThis will add the Doxygen keywords to c-mode and c++-mode only.

Default key bindings are:

   - C-c d ? will look up documentation for the symbol under the point.
   - C-c d r will rescan your Doxygen tags file.
   - C-c d f will insert a Doxygen comment for the next function.
   - C-c d i will insert a Doxygen comment for the current file.
   - C-c d ; will insert a Doxygen comment for the current member.
   - C-c d m will insert a blank multiline Doxygen comment.
   - C-c d s will insert a blank singleline Doxygen comment.")

;; if we can't find doxymacs.el then turn this mode off.
(if (file-installed-p "doxymacs.el")
    (setq want-doxymacs-mode t)
  (setq want-doxymacs-mode nil))

; fonts
(defvar want-rtt-font-lock-keywords t
"\nIf non-nil, the default set of keywords used by font-lock
will be overridden for some modes\n
   c++-mode:
   Additional keywords will be highlighted with color.
   (e.g.: Insist, Ensure, ...).")


;; CCS-4 directories
;; ========================================

;; The standard Draco environment files should be installed at 
;; /codes/radtran/vendors/environment in the directories:
;;
;;    elisp       XEmacs elisp scripts and settings
;;    templates   Templates for new documents
;;    bibfiles
;;    bibtex
;;    latex
;;    tex
;;    
(defvar ccs4-env-dir "/codes/radtran/vendors/environment/"
"\nDirectory that contains Draco environment files.
   - elisp       Subdirectory that contains XEmacs elisp scripts that
                 support the Draco development environment.
   - bibfiles    Subdirectory that contains LaTeX bibfiles for the
                 Draco development environment.
   - bibtex      Subdirectory that contains LaTeX style sheets for the
                 Draco development environment.
   - bin
   - doc
   - latex       Subdirectory that contains LaTeX style sheets for the
                 Draco development environment.
   - share
   - templates   Subdirectory that contains templates that support
                 rapid development of C++ source files.  These
                 templates are used when the user selects 
                 \"New file ...\" from the XEmacs RTT menu.
   - tex         currently empty.")

;; If the user specified value or the default value does not exist then
;; 1) look at the default location.
;; 2) print an error message.
(if (not (file-accessible-directory-p ccs4-env-dir))
    (if (file-accessible-directory-p "/usr/projects/draco/environment/")
	(setq ccs4-env-dir "/usr/projects/draco/environment/")
      (error (concat "Unable to find ccs4-env-dir = " ccs4-env-dir))))

(defvar ccs4-elisp-dir (concat ccs4-env-dir "elisp/")
"\nDirectory containing standard CC4 or Draco elisp code.")
(defvar ccs4-templates-dir (concat ccs4-env-dir "templates/")
"\nDirectory containing source code templates that are to be
used with the Draco elisp code (RTT Menu).")

;; Define a verbose load-library function
;; ========================================

(defun if-exist-load-library (libname)
  "If libname exists in load-path then load it.
If libname does not exist then print a warning message
and continue."
  (interactive)
    (if (load libname)
	(message (concat "Done loading file " libname ".el(c)"))
      (message (concat "Warning: " libname ".el not found"))))

;; Setup default Draco environment
;; ========================================

(cond ((string-match "XEmacs\\|Lucid" emacs-version)
; load tme's settings first then over ride some of them with my own.
; all of these files live in ccs4-elisp-dir

       ;; Draco setup
       (if (file-accessible-directory-p ccs4-elisp-dir)
	   (progn
             ; kill the default list and then set only the types we want.
	     (setq auto-mode-alist nil)

	     ; Load the Draco-rtt elisp packages:
	     ;
	     ; 1. This will load fl-keywords, rtt-hacks, tme-hacks, and
	     ;    draco-hacks. (rtt-hacks will load pooma-hacks and
             ;    infer-cc-style) 
	     ; 2. If want-rtt-menu is t then we load rtt-menu, 
             ;    Config-pkg (You must have all of of the config-pkg
	     ;    values alredy set!) and dante-macros.
	     (require 'draco-rtt)

	     ;; Set auto-mode-alist and autoload font-lock modes.
	     ;; Uses want-pkg-mode variables to determine what modes
	     ;; need to be configured using RTT defaults.
	     (require 'Config-pkg)
		 
	     ;; Custom key-bindings.
	     (load-library "Config-key")

	     ;; Lots of misc configuration setup including
	     ;; new-XXX-file commands
	     (if-exist-load-library "xemacs-setup")
	     )
	 (message "Warning: CCS-4 elisp directory not found no custom code loaded.")
	 ) ; end file-accessible-directory-p (ccs4-elisp-dir)
       (message "CCS-4 files loaded successfully.")))

(provide 'ccs4defaults)

;; ========================================
;; end of setup
;; ========================================
