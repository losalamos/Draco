;; ========================================
;; CCS-4 default settings
;; Kelly Thompson
;; 10 June 2002
;; ========================================

;; NOTE: ALL defvar values may be overriden in users's init.el.


;; Defaults for personal directories
;; ========================================

(defvar my-home-dir (concat (getenv "HOME") "/"))
(defvar my-elisp-dir (concat my-home-dir "lib/elisp/"))
(defvar my-templates-dir (concat my-elisp-dir "templates/"))
(defvar my-info-dir (concat my-elisp-dir "info/"))

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

; fonts
(defvar want-rtt-font-lock-keywords t
"\nIf non-nil, the default set of keywords used by font-lock
will be overridden for some modes\n
   c++-mode:
   Additional keywords will be highlighted with color.
   (e.g.: Insist, Ensure, ...).")


;; CCS-4 directories
;; ========================================

;; The standard CCS-4 elisp files should be installed at 
;; /codes/radtran/vendors/xemacs/elisp.
(defvar ccs4-elisp-dir "/codes/radtran/vendors/xemacs/elisp/")
;;(defvar ccs4-elisp-dir "/usr/projects/draco/elisp-head/")

;; If the user specified value or the default value does not exist then
;; 1) look at the default location.
;; 2) print an error message.
(if (not (file-accessible-directory-p ccs4-elisp-dir))
    (if (file-accessible-directory-p "/usr/projects/draco/vendors/xemacs/elisp/")
	(setq ccs4-elisp-dir "/usr/projects/draco/xemacs/elisp/")
      (error (concat "Unable to find ccs4-elisp-dir = " ccs4-elisp-dir))))

(setq ccs4-templates-dir (concat ccs4-elisp-dir "templates/"))
(setq ccs4-info-dir (concat ccs4-elisp-dir "info/"))

;; Add CCS-4 xemacs directories to load path
;; Assumes that ccs4-elisp-dir is valid

(setq load-path (cons ccs4-elisp-dir load-path))
(defvar Info-directory-list my-info-dir)
(setq Info-directory-list (cons ccs4-info-dir Info-directory-list))

(if (not (file-accessible-directory-p my-templates-dir))
    (if (file-accessible-directory-p ccs4-templates-dir)
	(setq my-templates-dir ccs4-templates-dir)))

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
; all of these files live in my-elisp-dir

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
