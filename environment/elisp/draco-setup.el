;;; -*-emacs-lisp-*-
;;; draco-setup.el --- ELisp package for making Draco related stuff easier.
;;
;; Author: Kelly Thompson
;; Created: 30 Nov 2004
;; Version: 0.0.1
;;
;; Copyright (C) 2016-2019 Triad National Security, LLC.
;; All rights reserved.
;;
;;============================================================
;; Commentary:
;;
;; - Optional: Customize the variable draco-env-dirs
;;
;; (setq draco-env-dirs (list
;;			    "/codes/radtran/vendors/environment/beta/"
;;			    "/codes/radtran/vendors/environment/"
;;			    "~/.xemacs" ))
;;
;; - Add the following line to your ~/.xemacs/init.el
;;
;;   (load-library "/fully/qualified/path/to/draco-setup")
;;
;; - Invoke draco-mode with M-x turn-on-draco-mode.  To have
;;   draco-mode invoked automatically when you start XEmacs add the
;;   following line to your init.el somewhere after you have loaded
;;   "draco-setup":
;;
;;   (turn-on-draco-mode)
;;
;; - Customize the mode by giving the command
;;   M-x customize-group <ret> draco-mode <ret>
;;
;; - NOTE: Fundamental mode will not autoload any of the draco-mode
;;   settings!  If you want all new XEmacs buffers to have these
;;   settings applied (including fonts/colors) you can set the
;;   default-mode in XEmacs to be text-mode.  Edit ~/.xemacs/custom.el
;;   and add the following command:
;;
;;   (custom-set-variables
;;      '(default-major-mode ('text-mode) t))
;;
;;============================================================

;; Customizable values (paths, basic setup options, etc.)

(defcustom my-home-dir (concat (getenv "HOME") "/")
"Location of $HOME."
:group 'draco-mode
:type 'string)

(defcustom draco-env-dirs nil "\nList of directories that will be
prepended to the load-path if they exist.  The directories <dir>/elisp
will also be examined and prepended to the the load-path if they
exist.  \nAdd to this list by using the following command in personal
elisp files: \n\t(setq draco-env-dirs (cons \"/path/to/extra/dir/\"))"
:group 'draco-mode
:type 'list)

(if (not draco-env-dirs)
    (setq draco-env-dirs
	  (list (concat my-draco-env-dir "elisp/")
                (concat my-home-dir "draco/environment/")
		(concat my-home-dir ".xemacs/")
		"/usr/projects/draco/environment/" )))

(defcustom draco-env-dir nil
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
                 \"New file ...\" from the XEmacs DRACO menu.
   - tex         currently empty."
:group 'draco-mode
:type 'string)

(defcustom draco-elisp-dir nil
"\nDirectory containing standard CC4 or Draco elisp code."
:group 'draco-mode
:type 'string)

(defcustom draco-templates-dir nil
"\nDirectory containing source code templates that are to be
used with the Draco elisp code (DRACO Menu)."
:group 'draco-mode
:type 'string)

(defcustom draco-colorize-modeline nil
  "*Does the user want us to colorize the modeline
based on the buffer-mode?  When customizing these colors
it may be useful to run the XEmacs command

M-x list-colors-display

to obtain a list of colors known to XEmacs."
:group 'draco-mode
:type '(radio	(const :tag "Yes" t)
		(const :tag "No"  nil)))

(defcustom draco-bibpath
  '("/codes/radtran/vendors/draco-5_7_0/draco/environment/bibfiles:/codes/radtran/vendors/capsaicin-2_5_0/source/doc/bib")
  "String containing a list of directories, each separated by a colon.
Each directory entry should end with a double slash.  Each directory
will be searched for bibfiles associated with \bibliography(file)
command at the end of a LaTeX document.  This path is used by
reftex-mode.;

Consider prepending your local directories."
  :group 'draco-mode
  :type 'list)

(defcustom draco-texpath
  '("/codes/radtran/vendors/draco-5_7_0/draco/environment/latex:/codes/radtran/vendors/draco-5_7_0/draco/environment/bibtex")
  "String containing a list of directories, each separated by a colon.
Each directory entry should end with a double slash.  Each directory
will be searched for LaTeX files (.sty, .bst, .cls, .tex, eps, etc.)
This path is used by reftex-mode.;

Consider prepending your local directories."
  :group 'draco-mode
  :type 'list)

(defcustom draco-code-comment-width 80
  "*Number of characters to use for comments (default 80)"
:group 'draco-mode
:type '(number) )

;; ========================================
;; Use Draco configuration for these modes
;; ========================================

(defcustom draco-want-mppl-mode t
  "*Does the user want to have the Draco minor mode enabled for MPPL
mode?"
:group 'draco-mode
:type '(radio	(const :tag "Yes" t)
		(const :tag "No"  nil)))

(defcustom draco-want-tcl-mode t
  "*Does the user want to have the Draco minor mode enabled for TCL
mode?"
:group 'draco-mode
:type '(radio	(const :tag "Yes" t)
		(const :tag "No"  nil)))

(defcustom draco-want-python-mode t
  "*Does the user want to have  the Draco minor mode enabled for
Python mode?"
:group 'draco-mode
:type '(radio	(const :tag "Yes" t)
		(const :tag "No"  nil)))

(defcustom draco-want-cmake-mode t
  "*Does the user want to have  the Draco minor mode enabled for CMake
mode?"
:group 'draco-mode
:type '(radio	(const :tag "Yes" t)
		(const :tag "No"  nil)))

(defcustom draco-want-autoconf-mode t
  "*Does the user want to have the Draco minor mode enabled for
autoconf mode?"
:group 'draco-mode
:type '(radio	(const :tag "Yes" t)
		(const :tag "No"  nil)))

(defcustom draco-want-makefile-mode t
  "*Does the user want to have the Draco minor mode enabled for
Makefile mode?"
:group 'draco-mode
:type '(radio	(const :tag "Yes" t)
		(const :tag "No"  nil)))

(defcustom draco-want-cc-mode t
  "*Does the user want to have the Draco minor mode enabled for C/C++
mode?"
:group 'draco-mode
:type '(radio	(const :tag "Yes" t)
		(const :tag "No"  nil)))

(defcustom draco-want-auctex-mode t
  "*Does the user want to have the Draco minor mode enabled for AucTeX
mode?"
:group 'draco-mode
:type '(radio	(const :tag "Yes" t)
		(const :tag "No"  nil)))

(defcustom draco-want-f90-mode t
  "*Does the user want to have the Draco minor mode enabled for
f90-mode?"
:group 'draco-mode)
:type '(radio	(const :tag "Yes" t)
		(const :tag "No"  nil))

(defcustom draco-want-fortran-mode t
  "*Does the user want to have the Draco minor mode enabled for
fortran-mode?"
:group 'draco-mode
:type '(radio	(const :tag "Yes" t)
		(const :tag "No"  nil)))

(defcustom draco-want-change-log-mode t
  "*Does the user want to have the Draco minor mode enabled for
change-log-mode?"
:group 'draco-mode
:type '(radio	(const :tag "Yes" t)
		(const :tag "No"  nil)))

(defcustom draco-want-emacs-lisp-mode t
  "*Does the user want to have the Draco minor mode enabled for
emacs-lisp-mode?"
:group 'draco-mode
:type '(radio	(const :tag "Yes" t)
		(const :tag "No"  nil)))

(defcustom draco-want-shell-mode t
  "*Does the user want to have the Draco minor mode enabled for
shell-mode?"
:group 'draco-mode
:type '(radio	(const :tag "Yes" t)
		(const :tag "No"  nil)))

(defcustom draco-want-cvs-mode t
  "*Does the user want to have the Draco minor mode enabled for
cvs-mode?"
:group 'draco-mode
:type '(radio	(const :tag "Yes" t)
		(const :tag "No"  nil)))

(defcustom draco-want-doxymacs-mode nil
  "*Does the user want to have the Draco minor mode enabled for
doxymacs-mode?"
:group 'draco-mode
:type '(radio	(const :tag "Yes" t)
		(const :tag "No"  nil)))

(defcustom draco-want-sh-mode t
  "*Does the user want to have the Draco minor mode enabled for
sh-mode?"
:group 'draco-mode
:type '(radio	(const :tag "Yes" t)
		(const :tag "No"  nil)))

(defcustom draco-want-text-mode t
  "*Does the user want to have the Draco minor mode enabled for
text-mode?"
:group 'draco-mode
:type '(radio	(const :tag "Yes" t)
		(const :tag "No"  nil)))

(defcustom draco-want-sgml-mode t
  "*Does the user want to have the Draco minor mode enabled for
sgml-mode?"
:group 'draco-mode
:type '(radio	(const :tag "Yes" t)
		(const :tag "No"  nil)))

(defcustom draco-want-dired-mode t
  "*Does the user want to have the Draco minor mode enabled for
dired-mode?"
:group 'draco-mode
:type '(radio	(const :tag "Yes" t)
		(const :tag "No"  nil)))

(defcustom draco-want-perl-mode t
  "*Does the user want to have the Draco minor mode enabled for
perl-mode?"
:group 'draco-mode
:type '(radio	(const :tag "Yes" t)
		(const :tag "No"  nil)))

(defcustom draco-want-compilation-mode t
  "*Does the user want to have the Draco minor mode enabled for
compilation-mode?"
:group 'draco-mode
:type '(radio	(const :tag "Yes" t)
		(const :tag "No"  nil)))

(defcustom draco-want-global-keys t
  "*Does the user want to have the Draco global key mappings?"
:group 'draco-mode
:type '(radio	(const :tag "Yes" t)
		(const :tag "No"  nil)))

;; ========================================
;; Update load path using draco-env-dirs.
;; Also, set draco-elisp-dir
;;========================================

(let ((dirlist draco-env-dirs))
  (progn
    ;; if the directory list is empty then guess some possible locations.
    (if (not draco-env-dirs)
	(setq draco-env-dirs
	      (list (concat my-home-dir "draco/environment/")
		    (concat my-home-dir "capsaicin/environment/")
		    (concat my-home-dir ".xemacs/")
		    "/usr/projects/draco/environment/"
		    "/ccs/codes/radtran/vendors/environment/" )))
    ;; process least-imporant directories first.  This ensures that
    ;; the first directory in the list "draco-env-dirs" become the
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
      ;; If we find the draco-mode.el file in this directory then set
      ;; this directory to be the draco-elisp-dir
      (if (file-readable-p (concat ldir "elisp/draco-setup.el"))
	  (progn
	    (setq draco-env-dir ldir)
	    (setq draco-elisp-dir (concat ldir "elisp/"))
	    (setq draco-templates-dir (concat ldir "templates/"))
	    ))
      ;; Set extra Info paths:
      (if (file-accessible-directory-p (concat ldir "info/"))
	  (setq Info-directory-list (cons (concat ldir "info/")
					  Info-directory-list )))

      ;; remove ldir from dirlist and continue the while-loop.
      (setq dirlist (cdr-safe dirlist)))
    ))

(require 'draco-faces)
(require 'draco-new-files)

;; ========================================
;; Setup some defaults
;; ========================================

;; Convert tabs to spaces when appropriate (c++-mode, etc.)

(defconst xemacsp (featurep 'xemacs) "Are we running XEmacs?")
(defconst emacs>=23p (and (not xemacsp) (> emacs-major-version 22))
  "Are we running GNU Emacs 23 or above?")

(setq-default indent-tabs-mode nil)
(defun insert-timestamp ()
  (interactive)
  (insert (format-time-string "%A, %b %d, %Y, %H:%M %P")))

;; Setup that only works for GNU Emacs
(if emacs>=23p
    (progn
      (which-function-mode   1)
      (global-font-lock-mode 1)
      )
  )

;; ========================================
;; Setup mode specific stuff
;; ========================================

; rebuild the auto-mode-alist from scratch using the draco
; information.
(setq auto-mode-alist nil)
(require (quote draco-config-modes))

;; If the bool draco-want-<pkg>-mode is t then setup mode <pkg>
;; using draco specific settings.  This includes turning on
;; draco-mode as a minor mode for each <pkg> mode.
(if draco-want-auctex-mode     (draco-setup-auctex-mode))
(if draco-want-autoconf-mode   (draco-setup-autoconf-mode))
(if draco-want-cc-mode         (draco-setup-cc-mode))
(if draco-want-change-log-mode (draco-setup-change-log-mode))
(if draco-want-compilation-mode (draco-setup-compilation-mode))
(if draco-want-cvs-mode        (draco-setup-cvs-mode))
(if draco-want-dired-mode      (draco-setup-dired-mode))
(if draco-want-perl-mode       (draco-setup-perl-mode))
(if draco-want-doxymacs-mode   (draco-setup-doxymacs-mode))
(if draco-want-emacs-lisp-mode (draco-setup-emacs-lisp-mode))
(if draco-want-f90-mode        (draco-setup-f90-mode))
(if draco-want-fortran-mode    (draco-setup-fortran-mode))
(if draco-want-makefile-mode   (draco-setup-makefile-mode))
;; (if draco-want-mppl-mode       (draco-setup-mppl-mode))
(if draco-want-python-mode     (draco-setup-python-mode))
;; (if draco-want-sgml-mode       (draco-setup-sgml-mode))

;; fontify/indent bash/csh scripts
(if emacs>=23p
    (if draco-want-sh-mode         (draco-setup-sh-mode)))

;(if draco-want-shell-mode      (draco-setup-shell-mode))
(if draco-want-tcl-mode        (draco-setup-tcl-mode))
(if draco-want-text-mode       (draco-setup-text-mode))
;; CMake needs to be set after text-mode so that the file associates
;; for CMakeLists.txt are set to cmake-mode instead of text-mode.
(if draco-want-cmake-mode      (draco-setup-cmake-mode))
(if draco-want-global-keys     (require 'draco-global-keys))

;(kill-buffer "*Compile-Log-Show*")
;(kill-buffer "*Compile-Log*")

;;; draco-setup.el ends here
