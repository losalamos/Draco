;;---------------------------------------------------------------------------;;
;; draco-new-files.el
;; Tom Evans
;; Feb. 17, 1999
;;---------------------------------------------------------------------------;;
;;; $Id$
;;---------------------------------------------------------------------------;;
;; provide macros for the draco development environment, we do not
;; provide key bindings here, that is for the user to decide
;; this file no longer requires tme-hacks.el
;;---------------------------------------------------------------------------;;

;;---------------------------------------------------------------------------;;
;; DRACO ENVIRONMENT FUNCTIONS (INTERACTIVE)
;;---------------------------------------------------------------------------;;
;;  1) setting up a package                              [draco-package]
;;  2) setting up a package test                         [draco-package-test]
;;  3) setting up a package autodoc dir                  [draco-package-doc]
;;  4) setting up a C++ translation unit                 [draco-class]
;;  5) setting up a C++ header                           [draco-cc-head]
;;  6) setting up a C++ header.in                        [draco-cc-headin]
;;  7) setting up a C header                             [draco-c-head]
;;  8) setting up a C header.in                          [draco-c-headin]
;;  9) setting up a C++ implementation file (.cc,.t.hh)  [draco-cc-imp]
;; 10) setting up a C++ instantiation file (_pt.cc)      [draco-cc-pt]
;; 11) setting up a test executable                      [draco-cc-test]
;; 12) setting up a python file                          [draco-python]
;; 13) setting up a specialized makefile                 [draco-make]
;; 14) setting up a Draco Memorandum (LaTeX)             [draco-memo]
;; 15) setting up a Draco Research Note (LaTeX)          [draco-note]
;; 16) setting up a Draco Article (LaTeX)                [draco-article]
;; 17) setting up a Draco Report (LaTeX)                 [draco-report]
;; 18) setting up a Draco Bibliography (LaTeX)           [draco-bib]
;; 19) setting up a Vision & Scope Doc. (LaTeX)          [draco-viscope]
;; 20) setting up a Bug Post-Mortem Memo (LaTeX)         [draco-bug-pm]
;;---------------------------------------------------------------------------;;

(require 'draco-helper-functions)

;;---------------------------------------------------------------------------;;
;; set up a draco package environment

(defun draco-package ()
  "Function to set up a draco package directory with stuff.
The files that we will place into all package directories are:

 1) CMakeLists.txt
 2) config.h.in

These files are based on templates in the draco/templates directory."

  (interactive)

  ;; Checks

;  (if (string= (draco-get-local-dir-name) "src")
;      (error "This command must be run from the <pkg> directory."))
;  (if (string= (draco-get-local-dir-name) "test")
;      (error "This command must be run from the <pkg> directory."))
;  (if (string= (draco-get-local-dir-name) "autodoc")
;      (error "This command must be run from the <pkg> directory."))

  ;; Query for package and namespace names and keep them function
  ;; local. 
  (let 
      ((draco-package-name nil)
       (draco-safe-package-name nil)
       (draco-namespace nil)
       (draco-paper-name nil))

    (draco-guess-names)

    ;; Create new buffers from templates (CMakeLists.txt, config.h.in)
    
    (draco-create-buffer-from-template 
     "CMakeLists.txt"
     (concat draco-templates-dir "/CMakeLists.package.txt"))
    
    (draco-create-buffer-from-template 
     "config.h.in"
     (concat draco-templates-dir "/template.h")
     "config")
    
    ))

;;---------------------------------------------------------------------------;;
;; set up a draco package test environment

(defun draco-package-test ()
  "Function to set up a draco package test directory with stuff.
The files that we will place into all package test directories are:

 1) CMakeLists.txt

These files are based on templates in the draco/templates directory."

  (interactive)
  ;; checks
  
  (if (not (string= (draco-get-local-dir-name) "test"))
      (error "This command must be run from the %s directory." 
	     (concat (draco-get-local-dir-name) "/test")))

  ;; Query for package and namespace names and keep them function
  ;; local. 
  (let 
      ((draco-package-name nil)
       (draco-safe-package-name nil)
       (draco-namespace nil)
       (draco-paper-name nil))

    (draco-guess-names)

  ;; Create new buffers from templates (CMakeLists.txt)

  (draco-create-buffer-from-template
   "CMakeLists.txt"
   (concat draco-templates-dir "/CMakeLists.test.txt"))
))

;;---------------------------------------------------------------------------;;
;; set up a draco package autodoc environment

(defun draco-package-doc ()
  "Function to set up a draco package autodoc directory with stuff.
The files that we will place into all package autodoc directories are:

 1) mainpage.dcc

These files are based on templates in the draco/templates directory."

  (interactive)

  ;; checks

  (if (not (string= (draco-get-local-dir-name) "autodoc"))
      (error "This command must be run from the %s directory." 
	     (concat (draco-get-local-dir-name) "/autodoc")))

  ;; Query for package and namespace names and keep them function
  ;; local. 
  (let 
      ((draco-package-name nil)
       (draco-safe-package-name nil)
       (draco-namespace nil)
       (draco-paper-name nil))

  (draco-guess-names)
  
  ;; Create new buffers from templates (mainpage.dcc):
  
  (draco-create-buffer-from-template 
   (concat draco-package-name ".dcc")
   (concat draco-templates-dir "/mainpage.dcc"))
  ))

;;---------------------------------------------------------------------------;;
;; make a new class (translation unit) in a draco pkg

(defun draco-class (draco-class-name)
  "Function to set a C++ translation unit in draco (.hh,.t.hh,_pt.cc,.cc)"
  
  ;; get the class draco-class-name
  (interactive "sClass Name: ")

  ;; Checks

;  (if (string= (draco-get-local-dir-name) "src")
;      (error "This command must be run from the <pkg> directory."))
;  (if (string= (draco-get-local-dir-name) "test")
;      (error "This command must be run from the <pkg> directory."))
;  (if (string= (draco-get-local-dir-name) "autodoc")
;      (error "This command must be run from the <pkg> directory."))

  ;; Query for package and namespace names and keep them function
  ;; local. 
  (let 
      ((draco-package-name nil)
       (draco-safe-package-name nil)
       (draco-namespace nil)
       (draco-paper-name nil))

  (draco-guess-names)

  ;; Create new buffers from templates (.hh, .cc, .i.hh)

  (draco-create-buffer-from-template 
   (concat draco-class-name ".hh")
   (concat draco-templates-dir "/template.hh")) 

  (draco-create-buffer-from-template
   (concat draco-class-name ".cc")
   (concat draco-templates-dir "/template.cc"))

  (draco-create-buffer-from-template
   (concat draco-class-name ".i.hh")
   (concat draco-templates-dir "/template.i.hh"))

  (draco-create-buffer-from-template
   (concat draco-class-name ".t.hh")
   (concat draco-templates-dir "/template.t.hh"))
))

;;---------------------------------------------------------------------------;;
;; set up a new draco C++ header file

(defun draco-cc-head (draco-class-name)
  "Function to set up a C++ header file.

This function will create two new buffers named draco-class-name.hh
and draco-class-name.i.hh.  The contents of the file template.hh
located at draco-templates-dir will be inserted and special text
markers will be replaced by appropriate strings."

  (interactive "sC++ Header Name: ")

  ;; Checks

;  (if (string= (draco-get-local-dir-name) "src")
;      (error "This command must be run from the <pkg> directory."))
;  (if (string= (draco-get-local-dir-name) "test")
;      (error "This command must be run from the <pkg> directory."))
;  (if (string= (draco-get-local-dir-name) "autodoc")
;      (error "This command must be run from the <pkg> directory."))

  ;; Query for package and namespace names and keep them function
  ;; local. 
  (let 
      ((draco-package-name nil)
       (draco-safe-package-name nil)
       (draco-namespace nil)
       (draco-paper-name nil))

    (draco-guess-names)

  ;; Create the buffer from the template file and replace the text.

  (draco-create-buffer-from-template 
   (concat draco-class-name ".hh")
   (concat draco-templates-dir "/template.hh")) 

  (draco-create-buffer-from-template
   (concat draco-class-name ".i.hh")
   (concat draco-templates-dir "/template.i.hh"))
))

;;---------------------------------------------------------------------------;;
;; set up a new C++ header.in file in draco

(defun draco-cc-headin (draco-header-name)
  "Function to set up a C++ header.in file.

This function will create a new buffer named draco-header-name.hh.in.
The contents of the file template.hh.in located at draco-templates-dir
will be inserted and special text markers will be replaced by
appropriate strings." 

  (interactive "sC++ Header.in Name: ")

  ;; Checks

;  (if (string= (draco-get-local-dir-name) "src")
;      (error "This command must be run from the <pkg> directory."))
;  (if (string= (draco-get-local-dir-name) "test")
;      (error "This command must be run from the <pkg> directory."))
;  (if (string= (draco-get-local-dir-name) "autodoc")
;      (error "This command must be run from the <pkg> directory."))

  ;; Query for package and namespace names and keep them function
  ;; local. 
  (let 
      ((draco-package-name nil)
       (draco-safe-package-name nil)
       (draco-namespace nil)
       (draco-paper-name nil))

    (draco-guess-names)

  ;; Create the buffer from the template file and replace the text.

  (draco-create-buffer-from-template 
   (concat draco-header-name ".hh.in")
   (concat draco-templates-dir "/template.hh")
   draco-header-name)
))

;;---------------------------------------------------------------------------;;
;; set up a new C header file in draco

(defun draco-c-head (draco-header-name)
  "Function to set up a C header file.

This function will create a new buffer named draco-header-name.h.  The
contents of the file template.h located at draco-templates-dir will be
inserted and special text markers will be replaced by appropriate
strings."

  (interactive "sC Header Name: ")

  ;; Checks

;  (if (string= (draco-get-local-dir-name) "src")
;      (error "This command must be run from the <pkg> directory."))
;  (if (string= (draco-get-local-dir-name) "test")
;      (error "This command must be run from the <pkg> directory."))
;  (if (string= (draco-get-local-dir-name) "autodoc")
;      (error "This command must be run from the <pkg> directory."))

  ;; Query for package and namespace names and keep them function
  ;; local. 
  (let 
      ((draco-package-name nil)
       (draco-safe-package-name nil)
       (draco-namespace nil)
       (draco-paper-name nil))

    (draco-guess-names)

  ;; Create the buffer from the template file and replace the text.

  (draco-create-buffer-from-template 
   (concat draco-header-name ".h")
   (concat draco-templates-dir "/template.h")
   draco-header-name)
))


;;---------------------------------------------------------------------------;;
;; set up a new C header.in file in draco

(defun draco-c-headin (draco-header-name)
  "Function to set up a C header.in file.

This function will create a new buffer named draco-header-name.h.in.
The contents of the file template.h.in located at draco-templates-dir
will be inserted and special text markers will be replaced by
appropriate strings."

  (interactive "sC Header.in Name: ")

  ;; Checks

;  (if (string= (draco-get-local-dir-name) "src")
;      (error "This command must be run from the <pkg> directory."))
;  (if (string= (draco-get-local-dir-name) "test")
;      (error "This command must be run from the <pkg> directory."))
;  (if (string= (draco-get-local-dir-name) "autodoc")
;      (error "This command must be run from the <pkg> directory."))

  ;; Query for package and namespace names and keep them function
  ;; local. 
  (let 
      ((draco-package-name nil)
       (draco-safe-package-name nil)
       (draco-namespace nil)
       (draco-paper-name nil))

    (draco-guess-names)

  ;; Create the buffer from the template file and replace the text.

  (draco-create-buffer-from-template
   (concat draco-header-name ".h.in")
   (concat draco-templates-dir "/template.h")
   draco-header-name)
))

;;---------------------------------------------------------------------------;;
;; set up a new C++ implementation file(s) [.cc,.i.hh] in draco

(defun draco-cc-imp (draco-class-name)
  "Function to set up C++ implementation files [.cc,.i.hh] in Draco.

This function will create two new buffers named draco-class-name.cc
and draco-class-name.i.hh.  The contents of the files template.cc and
template.i.hh located at draco-templates-dir will be inserted into
their respective new buffers and special text markers will be replaced
by appropriate strings."

  (interactive "sC++ Base Class Name: ")

  ;; Checks

;  (if (string= (draco-get-local-dir-name) "src")
;      (error "This command must be run from the <pkg> directory."))
;  (if (string= (draco-get-local-dir-name) "test")
;      (error "This command must be run from the <pkg> directory."))
;  (if (string= (draco-get-local-dir-name) "autodoc")
;      (error "This command must be run from the <pkg> directory."))

  ;; Query for package and namespace names and keep them function
  ;; local. 
  (let 
      ((draco-package-name nil)
       (draco-safe-package-name nil)
       (draco-namespace nil)
       (draco-paper-name nil))

    (draco-guess-names)

  ;; Create the buffer from the template file and replace the text.

  (draco-create-buffer-from-template 
   (concat draco-class-name ".cc")
   (concat draco-templates-dir "/template.cc"))

  (draco-create-buffer-from-template
   (concat draco-class-name ".i.hh")
   (concat draco-templates-dir "/template.i.hh"))

  
))

;;---------------------------------------------------------------------------;;
;; set up a new C++ instantiation file(s) [_pt.cc] in draco

(defun draco-cc-pt (draco-class-name)
  "Function to set up C++ instantiation files [_pt.cc] in Draco.

This function will create a new buffer named draco-class-name_pt.cc.
The contents of the file template_pt.cc at draco-templates-dir will be
inserted and special text markers will be replaced by appropriate
strings."

  (interactive "sC++ Instantiation Name: ")

  ;; Checks

;  (if (string= (draco-get-local-dir-name) "src")
;      (error "This command must be run from the <pkg> directory."))
;  (if (string= (draco-get-local-dir-name) "test")
;      (error "This command must be run from the <pkg> directory."))
;  (if (string= (draco-get-local-dir-name) "autodoc")
;      (error "This command must be run from the <pkg> directory."))

  ;; Query for package and namespace names and keep them function
  ;; local. 
  (let 
      ((draco-package-name nil)
       (draco-safe-package-name nil)
       (draco-namespace nil)
       (draco-paper-name nil))

    (draco-guess-names)

  ;; Create the buffer from the template file and replace the text.

  (draco-create-buffer-from-template 
   (concat draco-class-name "_pt.cc")
   (concat draco-templates-dir "/template_pt.cc"))
))

;;---------------------------------------------------------------------------;;
;; set up a new C++ test executable in draco

(defun draco-cc-test (draco-test-exe-name)
  "Function to set up C++ test executable in Draco.

This function will create a new buffer named draco-test-exe-name.cc.
The contents of the file template_c4_test.cc located at
draco-templates-dir will be inserted and special text markers will be
replaced by appropriate strings."

  (interactive "sC++ Test Executable Name: ")

  ;; checks
  
  (if (not (string= (draco-get-local-dir-name) "test"))
      (error "This command must be run from the %s directory." 
	     (concat (draco-get-local-dir-name) "/test")))

  ;; determine if this is parallel or not
  (defvar draco-test-exe-parallel (read-from-minibuffer "Parallel: " "y"))

  ;; Query for package and namespace names and keep them function
  ;; local. 
  (let 
      ((draco-package-name nil)
       (draco-safe-package-name nil)
       (draco-namespace nil)
       (draco-paper-name nil))

    (draco-guess-names)

  ;; Create the buffer from the template file and replace the text.

  (if (or (or (string= draco-test-exe-parallel "y")
	      (string= draco-test-exe-parallel "yes"))
	  (string= draco-test-exe-parallel "Y"))

      (draco-create-buffer-from-template
       (concat draco-test-exe-name ".cc")
       (concat draco-templates-dir "/template_c4_test.cc")
       draco-test-exe-name)
      (draco-create-buffer-from-template 
       (concat draco-test-exe-name ".cc")
       (concat draco-templates-dir "/template_test.cc")
       draco-test-exe-name))
))

;;---------------------------------------------------------------------------;;
;; set up a draco python file

(defun draco-python (draco-python-script-name)
  "Function to spontaneously setup a new Python file in draco"
  
  (interactive "sPython name: ")

  ;; Query for package and namespace names and keep them function
  ;; local. 
  (let 
      ((draco-package-name nil)
       (draco-safe-package-name nil)
       (draco-namespace nil)
       (draco-paper-name nil))

    (draco-guess-names)

  ;; Create the buffer from the template file and replace the text.

  (draco-create-buffer-from-template 
   (concat draco-python-script-name ".py")
   (concat draco-templates-dir "/template.py")
   draco-python-script-name)
))

;;---------------------------------------------------------------------------;;
;; set up a new misc. Makefile.in-type draco

(defun draco-make (draco-makefile-type)
  "Function to set up a Makefile.type file.

Currently \"type\" is one of \"in\" or \"temp.\""

  (interactive "sMakefile.type (type=[test|temp]): ")
;  (defvar draco-makefile-type (read-from-minibuffer "Type (test|temp): " "test")
;    "Makefile template used for creating a new Draco Makefile.
;One of \"test\" or \"temp\".")
  
  ;; If choice is not one of "test" or "temp", then abort.
  (if (not (or (string= draco-makefile-type "test")
	       (string= draco-makefile-type "temp")))
      (error "Makefile type must be one of test or temp."))

  ;; Query for package and namespace names and keep them function
  ;; local. 
  (let 
      ((draco-package-name nil)
       (draco-safe-package-name nil)
       (draco-namespace nil)
       (draco-paper-name nil))

    (draco-guess-names)

  ;; Create the buffer from the template file and replace the text.

  (draco-create-buffer-from-template 
   (concat "Makefile." draco-makefile-type)
   (concat draco-templates-dir "/Makefile." draco-makefile-type)
   (concat "Makefile." draco-makefile-type))
))

;;---------------------------------------------------------------------------;;
;; set up a memo

(defun draco-memo (name)
  "Function to spontaneously setup a new LaTeX LANL style memo"

  (interactive "sMemo Name: ")

  ;; Query for package and namespace names and keep them function
  ;; local. 
  (let 
      ((draco-package-name nil)
       (draco-safe-package-name nil)
       (draco-namespace nil)
       (draco-paper-name nil))

    (draco-guess-names "latex")
    
    ;; Create the buffer from the template file and replace the text.
    
    (draco-create-buffer-from-template 
     (concat draco-paper-name ".tex")
     (concat draco-templates-dir "/draco_memo.tex")
     draco-paper-name)
    ))

;;---------------------------------------------------------------------------;;
;; set up a research note

(defun draco-note (name)
  "Function to spontaneously setup a new LaTeX LANL style research note"

  (interactive "sResearch Note Name: ")

  ;; Query for package and namespace names and keep them function
  ;; local. 
  (let 
      ((draco-package-name nil)
       (draco-safe-package-name nil)
       (draco-namespace nil)
       (draco-paper-name nil))

    (draco-guess-names "latex")
    
    ;; Create the buffer from the template file and replace the text.
    
    (draco-create-buffer-from-template 
     (concat draco-paper-name ".tex")
     (concat draco-templates-dir "/draco_note.tex")
     draco-paper-name)
    ))

;;---------------------------------------------------------------------------;;
;; set up an article based on the LANL style of TME

(defun draco-article (name)
  "Function to spontaneously setup a new LaTeX LANL style article"

  (interactive "sPaper Name: ")

  ;; Query for package and namespace names and keep them function
  ;; local. 
  (let 
      ((draco-package-name nil)
       (draco-safe-package-name nil)
       (draco-namespace nil)
       (draco-paper-name nil))

    (draco-guess-names "latex")
    
    ;; Create the buffer from the template file and replace the text.
    
    (draco-create-buffer-from-template 
     (concat draco-paper-name ".tex")
     (concat draco-templates-dir "/draco_art.tex")
     draco-paper-name)
    ))

;;---------------------------------------------------------------------------;;
;; set up a draco report

(defun draco-report (name)
  "Function to spontaneously setup a new LaTeX LANL style report"

  (interactive "sReport Name: ")

  ;; Query for package and namespace names and keep them function
  ;; local. 
  (let 
      ((draco-package-name nil)
       (draco-safe-package-name nil)
       (draco-namespace nil)
       (draco-paper-name nil))

    (draco-guess-names "latex")
    
    ;; Create the buffer from the template file and replace the text.
    
    (draco-create-buffer-from-template 
     (concat draco-paper-name ".tex")
     (concat draco-templates-dir "/draco_rep.tex")
     draco-paper-name)
    ))

;;---------------------------------------------------------------------------;;
;; set up a new bib file

(defun draco-bib (name)
  "Function to spontaneously setup a new LaTeX BiBTeX file"

  (interactive "sBib Name: ")

  ;; Query for package and namespace names and keep them function
  ;; local. 
  (let 
      ((draco-package-name nil)
       (draco-safe-package-name nil)
       (draco-namespace nil)
       (draco-paper-name nil))

    (draco-guess-names "bibtex")
    
    ;; Create the buffer from the template file and replace the text.
    
    (draco-create-buffer-from-template 
     (concat draco-paper-name ".bib")
     (concat draco-templates-dir "/draco_bib.bib")
     draco-paper-name)
))

;;---------------------------------------------------------------------------;;
;; set up a "project vision/scope statement" memo 

(defun draco-viscope (name)
  "Function to spontaneously setup a new vision/scope memo"

  (interactive "sVision/Scope Memo Name: ")

  ;; Query for package and namespace names and keep them function
  ;; local. 
  (let 
      ((draco-package-name nil)
       (draco-safe-package-name nil)
       (draco-namespace nil)
       (draco-paper-name nil))

    (draco-guess-names "latex")
    
    ;; Create the buffer from the template file and replace the text.
    
    (draco-create-buffer-from-template 
     (concat draco-paper-name ".tex")
     (concat draco-templates-dir "/draco_viscope.tex")
     draco-paper-name)
    ))

;;---------------------------------------------------------------------------;;
;; set up a "Bug Post-Mortem" memo 

(defun draco-bug-pm (name)
  "Function to spontaneously setup a new bug post-mortem memo"

  (interactive "sMemo Name: ")

  ;; Query for package and namespace names and keep them function
  ;; local. 
  (let 
      ((draco-package-name nil)
       (draco-safe-package-name nil)
       (draco-namespace nil)
       (draco-paper-name nil))

    (draco-guess-names "latex")
    
    ;; Create the buffer from the template file and replace the text.
    
    (draco-create-buffer-from-template 
     (concat draco-paper-name ".tex")
     (concat draco-templates-dir "/draco_bug_pm.tex")
     draco-paper-name)
))

;;---------------------------------------------------------------------------;;

(provide 'draco-new-files)
