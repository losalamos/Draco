;;; draco-hacks.el
;;; Tom Evans
;;; Feb. 17, 1999

;;; $Id$

;;---------------------------------------------------------------------------;;
;; provide macros for the draco development environment, we do not
;; provide key bindings here, that is for the user to decide
;; this file no longer requires tme-hacks.el
;;---------------------------------------------------------------------------;;

;;---------------------------------------------------------------------------;;
;; SERVICE FUNCTIONS USED IN DRACO-HACKS.EL
;;---------------------------------------------------------------------------;;

(defvar my-templates-dir nil
  "*directory location of template files (alternates to draco/templates).")

;;---------------------------------------------------------------------------;;
;; Query for the pkg name, with guessed default.

(defun draco-pkg-header ()
  "Function to get a pkg directory name"
  ;;(setq pkg " ")
  ;;(read-from-minibuffer "Package name: " (substring pkg 1)))
  (setq dir (expand-file-name "."))
  (setq parent (expand-file-name ".."))
  
  ;; Figure out the difference in length between dir and parent,
  ;; subtract one (for the "/"), and negate, in order to get that many
  ;; chars off the end of dir.
  
  (setq dlen (- (length dir) (length parent)))
  (setq xlen (* (- dlen 1) -1))
  
  (read-from-minibuffer "Package path: " (substring dir xlen)))

;;---------------------------------------------------------------------------;;
;; Query for the pkg name, with guessed default.

(defun draco-pkg-test-header ()
  "Function to get a pkg/test directory name"
  ;;(setq pkg " ")
  ;;(read-from-minibuffer "Package name: " (substring pkg 1)))
  (setq dir (expand-file-name "."))
  (setq parent (expand-file-name "../.."))
  
  ;; Figure out the difference in length between dir and parent,
  ;; subtract one (for the "/"), and negate, in order to get that many
  ;; chars off the end of dir.
  
  (setq dlen (- (length dir) (length parent)))
  (setq xlen (* (- dlen 1) -1))
  
  (read-from-minibuffer "Package test path: " (substring dir xlen)))

;;---------------------------------------------------------------------------;;
;; Query for the pkg name, with guessed default.

(defun draco-pkg-doc-header ()
  "Function to get a pkg name from pkg/autodoc directory"
  ;;(setq pkg " ")
  ;;(read-from-minibuffer "Package name: " (substring pkg 1)))
  (setq dir (expand-file-name ".."))
  (setq parent (expand-file-name "../.."))
  
  ;; Figure out the difference in length between dir and parent,
  ;; subtract one (for the "/"), and negate, in order to get that many
  ;; chars off the end of dir.
  (setq dlen (- (length dir) (length parent)))
  (setq xlen (* (- dlen 1) -1))
  
  (read-from-minibuffer "Package doc name: " (substring dir xlen)))

;;---------------------------------------------------------------------------;;
;; Query for the pkg name, with guessed default.

(defun draco-pkg-in-test-header ()
  "Function to get a pkg/test directory name"
  ;;(setq pkg " ")
  ;;(read-from-minibuffer "Package name: " (substring pkg 1)))
  (setq dir (expand-file-name ".."))
  (setq parent (expand-file-name "../.."))
  
  ;; Figure out the difference in length between dir and parent,
  ;; subtract one (for the "/"), and negate, in order to get that many
  ;; chars off the end of dir.
  
  (setq dlen (- (length dir) (length parent)))
  (setq xlen (* (- dlen 1) -1))
  
  (read-from-minibuffer "Package for this test: " (substring dir xlen)))
  
;;---------------------------------------------------------------------------;;
;; Query for the safe pkg name from a package name.

(defun draco-safe-pkg (pkg)
  "Function to get a safe pkg from a package name"
  (setq spkg (mapconcat (function (lambda (x)
				    (cond
				     ((eq x ?\+) "")
				     ((eq x ?/) "_")
				     (t (format "%c" x)))))
			pkg "" )))

;;---------------------------------------------------------------------------;;
;; Query for the desired namespace

(defun draco-namespace (guess)
  "Function to query the user for the desired namespace given a default guess"
  
  ;; now get the namespace
  (read-from-minibuffer "Namespace for this package: " guess))

;;---------------------------------------------------------------------------;;
;; find where draco is; we assume we are in a draco sub-directory

(defun draco-dracodir ()
  "Function to determine the draco/ directory location"
  
  ;; see if we are in draco directory, we go a maximum of 20 levels

  ;; definitions of variables
  (setq ldir ".")
  (setq pdir "..")
  (setq found nil)
  (setq goal "draco")
  (setq nlevel 0)
  (setq level nil)

  (while (not found)
    ;; get this level and the parent level
    (setq level (expand-file-name ldir))
    (setq parent (expand-file-name pdir))

    ;; determine lengths of this level and the parent level
    (setq dlen (- (length level) (length parent)))
    (setq xlen (* (- dlen 1) -1))

    ;; get the name of this level and check to see if it is draco
    (setq pkg (substring level xlen))
    (if (string= pkg goal)
	(setq found t))

    ;; go up levels if we haven't found draco
    (setq ldir (concat ldir "/.."))
    (setq pdir (concat pdir "/.."))

    ;; we only go a maximum of 20 levels, so use a counter
    (setq nlevel (+ nlevel 1))
    (cond ((= nlevel 20)
	   (error "Can't find draco/ !!!!\n")
	   (setq found t)))
    )

  ;; return the draco directory
  level)

;;---------------------------------------------------------------------------;;
;; find where draco/templates is; assume we are in a draco sub-directory

(defun rtt-templates ()
  "Find the nearest templates directory.  Look for a templates directory
under the current framework (draco, solon, dante, etc).  If no
templates directory is found then try (my-templates-directory)."
  
  ;; find the draco/environment/templates directory go a maximum of
  ;; twenty levels

  ;; definitions of variables
  (setq ldir ".")
  (setq found nil)
  (setq nlevel 0)

  (while (not found)
    ;; get this level and the parent level
    (setq level (expand-file-name ldir))
    (setq target (concat level "/environment/templates"))

    ;; see if the environment/templates/ directory exists
    (if (file-exists-p target)
	(setq found t))

    ;; go up levels if we haven't found draco
    (setq ldir (concat ldir "/.."))

    ;; we only go a maximum of 20 levels, so use a counter
    (setq nlevel (+ nlevel 1))

    ;; If we reach level 20 and haven't found the templates
    ;; directory then we look for an alternate templates directory
    ;; specified by (my-templates-dir).  If this is a valid
    ;; directory then return its value otherwise return an error.
    (if (= nlevel 20)
	(progn ; directory not found
	  (message "Can't find draco/environment/templates. Looking for alt. template directory.\n")
	  (if my-templates-dir ; has this var been set?
	      (if (file-accessible-directory-p my-templates-dir)
		  (progn
		    (setq target my-templates-dir)
		    (setq found t))))
          ;; if we still don't have a good template directory then
	  ;; return an error.
	  (if (not found) 
	      (progn
		(error "No template directory found!\n")
		(setq found t)))
	 ))) 
  
  ;; return the target
  target)

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
;; 11) setting up a python file                          [draco-python]
;; 12) setting up a specialized makefile                 [draco-make]
;; 13) setting up a test executable                      [draco-cc-test]

;;---------------------------------------------------------------------------;;
;; set up a draco package environment

(defun draco-package ()
  "Function to set up a draco package directory with stuff"

  (interactive)

  ;; the files that we will place into all package directories are:
  ;;
  ;; 1) configure.in
  ;; 2) config.h.in
  ;; 3) Release.hh
  ;; 4) Release.cc
  ;;
  ;; these files are based on templates in the draco/templates
  ;; directory that can be accessed using the rtt-templates() elisp 
  ;; function

  ;; first find templates directory
  (setq tempdir (rtt-templates))

  ;; now get the header for this pkg
  (setq pkg (draco-pkg-header))

  ;; now get the safe pkg name
  (setq spkg (draco-safe-pkg pkg))

  ;; now get the namespace
  (setq nspace (draco-namespace spkg))

  ;; now we get the files that the directory does not currently have
  (setq conffile "configure.ac")
  (setq chinfile "config.h.in")
  (setq vhhfile "Release.hh")
  (setq vccfile "Release.cc")
  
  ;;
  ;; first lets get the configure.in file
  ;;
  (cond ((not (file-exists-p conffile))
	 (setq tempcon (concat tempdir "/configure.package.ac"))
	 (if (file-exists-p conffile)
	     (error "File configure.ac already exists." conffile))
	 
	 (find-file conffile)
	 (insert-file-contents tempcon)
	 
	 ;; insert the goodies into the new configure.in
	 (find-file conffile)
	 (perform-replace "<user>" (user-full-name) nil nil nil )
	 (goto-char (point-min))
	 (perform-replace "<pkg>" pkg nil nil nil )
	 (goto-char (point-min))
	 (perform-replace "<date>"  (current-time-string) nil nil nil)
	 (goto-char (point-min))
	 (perform-replace "<start>" "" nil nil nil )))

  ;;
  ;; now lets get the config.h.in file
  ;;
  (cond ((not (file-exists-p chinfile))
	 (setq tempchin (concat tempdir "/template.h"))
	 (if (file-exists-p chinfile)
	     (error "File config.h.in already exists." chinfile))
	 
	 (find-file chinfile)
	 (insert-file-contents tempchin)
	 
	 ;; insert the goodies into the new config.h.in
	 (find-file chinfile)
	 (perform-replace "<user>" (user-full-name) nil nil nil )
	 (goto-char (point-min))
	 (perform-replace "<class>" "config" nil nil nil )
	 (goto-char (point-min))
	 (perform-replace "<pkg>" pkg nil nil nil )
	 (goto-char (point-min))
	 (perform-replace "<spkg>" spkg nil nil nil )
	 (goto-char (point-min))
	 (perform-replace "<date>"  (current-time-string) nil nil nil)
	 (goto-char (point-min))
	 (perform-replace "<start>" "" nil nil nil )))
  
  ;; 
  ;; now lets make the Release.hh file
  ;;
  (cond ((not (file-exists-p vhhfile))
	 (setq tempvhh (concat tempdir "/Release.hh"))
	 (if (file-exists-p vhhfile)
	     (error "File Release.hh already exists." vhhfile))
	 
	 (find-file vhhfile)
	 (insert-file-contents tempvhh)
	 
	 ;; insert the goodies into the new config.h.in
	 (find-file vhhfile)
	 (perform-replace "<user>" (user-full-name) nil nil nil )
	 (goto-char (point-min))
	 (perform-replace "<pkg>" pkg nil nil nil )
	 (goto-char (point-min))
	 (perform-replace "<spkg>" spkg nil nil nil )
	 (goto-char (point-min))
	 (perform-replace "<namespace>" nspace nil nil nil )
	 (goto-char (point-min))
	 (perform-replace "<date>"  (current-time-string) nil nil nil)
	 (goto-char (point-min))
	 (perform-replace "<start>" "" nil nil nil )))

  ;;
  ;; now lets make the Release.cc file
  ;;
  (cond ((not (file-exists-p vccfile))
	 (setq tempvcc (concat tempdir "/Release.cc"))
	 (if (file-exists-p vccfile)
	     (error "File Release.cc already exists." vccfile))
	 
	 (find-file vccfile)
	 (insert-file-contents tempvcc)
	 
	 ;; insert the goodies into the new config.h.in
	 (find-file vccfile)
	 (perform-replace "<user>" (user-full-name) nil nil nil )
	 (goto-char (point-min))
	 (perform-replace "<pkg>" pkg nil nil nil )
	 (goto-char (point-min))
	 (perform-replace "<spkg>" spkg nil nil nil )
	 (goto-char (point-min))
	 (perform-replace "<namespace>" nspace nil nil nil )
	 (goto-char (point-min))
	 (perform-replace "<date>"  (current-time-string) nil nil nil)
	 (goto-char (point-min))
	 (perform-replace "<start>" "" nil nil nil ))))

;;---------------------------------------------------------------------------;;
;; set up a draco package test environment

(defun draco-package-test ()
  "Function to set up a draco package test directory with stuff"

  (interactive)

  ;; the files that we will place into all package test directories are:
  ;;
  ;; 1) <spkg>.hh
  ;; 2) <spkg>.cc
  ;; 3) Makefile.target
  ;;
  ;; these files are based on templates in the draco/templates
  ;; directory that can be accessed using the rtt-templates() elisp 
  ;; function

  ;; first find templates directory
  (setq tempdir (rtt-templates))

  ;; now get the header for this pkg
  (setq pkg (draco-pkg-test-header))

  ;; now get the safe pkg name
  (setq spkg (draco-safe-pkg pkg))

  ;; now get the namespace
  (setq nspace (draco-namespace spkg))

  ;; now we get the files that the directory does not currently have
  (setq header (concat spkg ".hh"))
  (setq source (concat spkg ".cc"))
  (setq makefile "Makefile.target")

  ;;
  ;; first lets get the header file
  ;;
  (cond ((not (file-exists-p header))
	 (setq temphead (concat tempdir "/pkg_Test.hh"))
	 (if (file-exists-p header);
	     (error "Header file already exists." header))
	 
	 (find-file header)
	 (insert-file-contents temphead)
	 
	 ;; insert the goodies into the new header
	 (find-file header)
	 (perform-replace "<user>" (user-full-name) nil nil nil )
	 (goto-char (point-min))
	 (perform-replace "<pkg>" pkg nil nil nil )
	 (goto-char (point-min))
	 (perform-replace "<date>"  (current-time-string) nil nil nil)
	 (goto-char (point-min))
	 (perform-replace "<spkg>" spkg nil nil nil )
	 (goto-char (point-min))
	 (perform-replace "<namespace>" nspace nil nil nil )
	 (goto-char (point-min))
	 (perform-replace "<start>" "" nil nil nil )))

  ;;
  ;; next lets get the source file
  ;;
  (cond ((not (file-exists-p source))
	 (setq tempsrc (concat tempdir "/pkg_Test.cc"))
	 (if (file-exists-p source);
	     (error "Source file already exists." source))
	 
	 (find-file source)
	 (insert-file-contents tempsrc)
	 
	 ;; insert the goodies into the new source
	 (find-file source)
	 (perform-replace "<user>" (user-full-name) nil nil nil )
	 (goto-char (point-min))
	 (perform-replace "<pkg>" pkg nil nil nil )
	 (goto-char (point-min))
	 (perform-replace "<date>"  (current-time-string) nil nil nil)
	 (goto-char (point-min))
	 (perform-replace "<spkg>" spkg nil nil nil )
	 (goto-char (point-min))
	 (perform-replace "<namespace>" nspace nil nil nil )
	 (goto-char (point-min))
	 (perform-replace "<start>" "" nil nil nil )))

  ;;
  ;; finally lets get the makefile
  ;;
  (cond ((not (file-exists-p makefile))
	 (setq tempmake (concat tempdir "/Makefile.test"))
	 (if (file-exists-p makefile);
	     (error "Makefile.target already exists." makefile))
	 
	 (find-file makefile)
	 (insert-file-contents tempmake)
	 
	 ;; insert the goodies into the new makefile
	 (find-file makefile)
	 (perform-replace "<user>" (user-full-name) nil nil nil )
	 (goto-char (point-min))
	 (perform-replace "<pkg>" pkg nil nil nil )
	 (goto-char (point-min))
	 (perform-replace "<date>"  (current-time-string) nil nil nil)
	 (goto-char (point-min))
	 (perform-replace "<spkg>" spkg nil nil nil )
	 (goto-char (point-min))
	 (perform-replace "<start>" "" nil nil nil )))
  )

;;---------------------------------------------------------------------------;;
;; set up a draco package autodoc environment

(defun draco-package-doc ()
  "Function to set up a draco package autodoc directory with stuff"

  (interactive)

  ;; the files that we will place into all package autodoc directories are:
  ;;
  ;; 1) mainpage.dcc
  ;;
  ;; these files are based on templates in the draco/templates
  ;; directory that can be accessed using the rtt-templates() elisp 
  ;; function

  ;; first find templates directory
  (setq tempdir (rtt-templates))

  ;; now get the header for this pkg
  (setq pkg (draco-pkg-doc-header))

  ;; now get the safe pkg name
  (setq spkg (draco-safe-pkg pkg))

  ;; now we get the files that the directory does not currently have
  (setq mainpage (concat pkg ".dcc"))

  ;;
  ;; get the mainpage
  ;;
  (cond ((not (file-exists-p mainpage))
	 (setq tempmain (concat tempdir "/mainpage.dcc"))
	 (if (file-exists-p mainpage);
	     (error "mainpage already exists." mainpage))
	 
	 (find-file mainpage)
	 (insert-file-contents tempmain)
	 
	 ;; insert the goodies into the new makefile
	 (find-file mainpage)
	 (perform-replace "<user>" (user-full-name) nil nil nil )
	 (goto-char (point-min))
	 (perform-replace "<pkg>" pkg nil nil nil )
	 (goto-char (point-min))
	 (perform-replace "<date>"  (current-time-string) nil nil nil)
	 (goto-char (point-min))
	 (perform-replace "<spkg>" spkg nil nil nil )
	 (goto-char (point-min))
	 (perform-replace "<start>" "" nil nil nil )))
  )

;;---------------------------------------------------------------------------;;
;; make a new class (translation unit) in a draco pkg

(defun draco-class (name)
  "Function to set a C++ translation unit in draco (.hh,.t.hh,.cc)"
  
  ;; get the class name
  (interactive "sClass Name: ")

  ;; get the package directory
  (setq pkg (draco-pkg-header))

  ;; now get the safe pkg name
  (setq spkg (draco-safe-pkg pkg))

  ;; now get the namespace
  (setq nspace (draco-namespace spkg))

  ;; find templates directory
  (setq tempdir (rtt-templates))

  ;;
  ;; get the proper template for the header file
  ;;

  (setq nfile (concat name ".hh"))
  (setq tfile (concat tempdir "/template.hh"))
  (if (file-exists-p nfile)
      (error "File %s already exists." nfile))

  (find-file nfile)
  (insert-file-contents tfile)

  ;; insert the goodies
  
  (find-file nfile)
  (perform-replace "<pkg>" pkg nil nil nil )
  (goto-char (point-min))
  (perform-replace "<spkg>" spkg nil nil nil )
  (goto-char (point-min))
  (perform-replace "<namespace>" nspace nil nil nil )
  (goto-char (point-min))
  (perform-replace "<user>" (user-full-name) nil nil nil )
  (goto-char (point-min))
  (perform-replace "<class>" name nil nil nil )
  (goto-char (point-min))
  (perform-replace "<date>"  (current-time-string) nil nil nil )
  (goto-char (point-min))
  (perform-replace "<start>" "" nil nil nil )

  ;;
  ;; get the proper template for the implementation file (.cc)
  ;;

  (setq nfile (concat name ".cc"))
  (setq tfile (concat tempdir "/template.cc"))
  (if (file-exists-p nfile)
      (error "File %s already exists." nfile))

  (find-file nfile)
  (insert-file-contents tfile)

  ;; insert the goodies
  
  (find-file nfile)
  (perform-replace "<pkg>" pkg nil nil nil )
  (goto-char (point-min))
  (perform-replace "<spkg>" spkg nil nil nil )
  (goto-char (point-min))
  (perform-replace "<namespace>" nspace nil nil nil )
  (goto-char (point-min))
  (perform-replace "<user>" (user-full-name) nil nil nil )
  (goto-char (point-min))
  (perform-replace "<class>" name nil nil nil )
  (goto-char (point-min))
  (perform-replace "<basename>" name nil nil nil )
  (goto-char (point-min))
  (perform-replace "<date>"  (current-time-string) nil nil nil )
  (goto-char (point-min))
  (perform-replace "<start>" "" nil nil nil )

  ;;
  ;; get the proper template for the implementation file (.t.hh)
  ;;

  (setq nfile (concat name ".t.hh"))
  (setq tfile (concat tempdir "/template.t.hh"))
  (if (file-exists-p nfile)
      (error "File %s already exists." nfile))

  (find-file nfile)
  (insert-file-contents tfile)

  ;; insert the goodies
  
  (find-file nfile)
  (perform-replace "<pkg>" pkg nil nil nil )
  (goto-char (point-min))
  (perform-replace "<spkg>" spkg nil nil nil )
  (goto-char (point-min))
  (perform-replace "<namespace>" nspace nil nil nil )
  (goto-char (point-min))
  (perform-replace "<user>" (user-full-name) nil nil nil )
  (goto-char (point-min))
  (perform-replace "<class>" name nil nil nil )
  (goto-char (point-min))
  (perform-replace "<basename>" name nil nil nil )
  (goto-char (point-min))
  (perform-replace "<date>"  (current-time-string) nil nil nil )
  (goto-char (point-min))
  (perform-replace "<start>" "" nil nil nil )

  ;;
  ;; get the proper template for the implementation file (.t.hh)
  ;;

  (setq nfile (concat name ".i.hh"))
  (setq tfile (concat tempdir "/template.i.hh"))
  (if (file-exists-p nfile)
      (error "File %s already exists." nfile))

  (find-file nfile)
  (insert-file-contents tfile)

  ;; insert the goodies
  
  (find-file nfile)
  (perform-replace "<pkg>" pkg nil nil nil )
  (goto-char (point-min))
  (perform-replace "<spkg>" spkg nil nil nil )
  (goto-char (point-min))
  (perform-replace "<namespace>" nspace nil nil nil )
  (goto-char (point-min))
  (perform-replace "<user>" (user-full-name) nil nil nil )
  (goto-char (point-min))
  (perform-replace "<class>" name nil nil nil )
  (goto-char (point-min))
  (perform-replace "<date>"  (current-time-string) nil nil nil )
  (goto-char (point-min))
  (perform-replace "<start>" "" nil nil nil ))

;;---------------------------------------------------------------------------;;
;; set up a new draco C++ header file

(defun draco-cc-head (name)
  "Function to set up a C++ header file"

  (interactive "sC++ Header Name: ")

  ;; get package directory
  (setq pkg (draco-pkg-header))
  
  ;; now get the safe pkg name
  (setq spkg (draco-safe-pkg pkg))

  ;; now get the namespace
  (setq nspace (draco-namespace spkg))

  (setq tempdir (rtt-templates))
  (setq nfile (concat name ".hh"))
  (setq tfile (concat tempdir "/template.hh"))
  (if (file-exists-p nfile)
      (error "File %s already exists." nfile))
  
  (find-file nfile)
  (insert-file-contents tfile)

  ;; do Furnish style customizations

  (find-file nfile)
  (perform-replace "<pkg>" pkg nil nil nil )
  (goto-char (point-min))
  (perform-replace "<spkg>" spkg nil nil nil )
  (goto-char (point-min))
  (perform-replace "<namespace>" nspace nil nil nil )
  (goto-char (point-min))
  (perform-replace "<user>" (user-full-name) nil nil nil )
  (goto-char (point-min))
  (perform-replace "<class>" name nil nil nil )
  (goto-char (point-min))
  (perform-replace "<date>"  (current-time-string) nil nil nil )
  (goto-char (point-min))
  (perform-replace "<start>" "" nil nil nil )

  ;; Add a Class.i.hh file
  (setq nfile (concat name ".i.hh"))
  (setq tfile (concat tempdir "/template.i.hh"))
  (if (file-exists-p nfile)
      (error "File %s already exists." nfile))
  
  (find-file nfile)
  (insert-file-contents tfile)

  ;; do Furnish style customizations

  (find-file nfile)
  (perform-replace "<pkg>" pkg nil nil nil )
  (goto-char (point-min))
  (perform-replace "<spkg>" spkg nil nil nil )
  (goto-char (point-min))
  (perform-replace "<namespace>" nspace nil nil nil )
  (goto-char (point-min))
  (perform-replace "<user>" (user-full-name) nil nil nil )
  (goto-char (point-min))
  (perform-replace "<class>" name nil nil nil )
  (goto-char (point-min))
  (perform-replace "<date>"  (current-time-string) nil nil nil )
  (goto-char (point-min))
  (perform-replace "<start>" "" nil nil nil ))

;;---------------------------------------------------------------------------;;
;; set up a new C++ header.in file in draco

(defun draco-cc-headin (name)
  "Function to set up a C++ header.in file"

  (interactive "sC++ Header.in Name: ")

  ;; get package directory
  (setq pkg (draco-pkg-header))
  
  ;; now get the safe pkg name
  (setq spkg (draco-safe-pkg pkg))

  ;; now get the namespace
  (setq nspace (draco-namespace spkg))

  (setq tempdir (rtt-templates))
  (setq nfile (concat name ".hh.in"))
  (setq tfile (concat tempdir "/template.hh"))
  (if (file-exists-p nfile)
      (error "File %s already exists." nfile))
  
  (find-file nfile)
  (insert-file-contents tfile)

  ;; do Furnish style customizations

  (find-file nfile)
  (perform-replace "<pkg>" pkg nil nil nil )
  (goto-char (point-min))
  (perform-replace "<spkg>" spkg nil nil nil )
  (goto-char (point-min))
  (perform-replace "<namespace>" nspace nil nil nil )
  (goto-char (point-min))
  (perform-replace "<user>" (user-full-name) nil nil nil )
  (goto-char (point-min))
  (perform-replace "<class>" name nil nil nil )
  (goto-char (point-min))
  (perform-replace "<date>"  (current-time-string) nil nil nil )
  (goto-char (point-min))
  (perform-replace "<start>" "" nil nil nil ))

;;---------------------------------------------------------------------------;;
;; set up a new C header file in draco

(defun draco-c-head (name)
  "Function to set up a C header file"

  (interactive "sC Header Name: ")

  ;; get package directory
  (setq pkg (draco-pkg-header))
  
  ;; now get the safe pkg name
  (setq spkg (draco-safe-pkg pkg))

  (setq tempdir (rtt-templates))
  (setq nfile (concat name ".h"))
  (setq tfile (concat tempdir "/template.h"))
  (if (file-exists-p nfile)
      (error "File %s already exists." nfile))
  
  (find-file nfile)
  (insert-file-contents tfile)

  ;; do Furnish style customizations

  (find-file nfile)
  (perform-replace "<pkg>" pkg nil nil nil )
  (goto-char (point-min))
  (perform-replace "<spkg>" spkg nil nil nil )
  (goto-char (point-min))
  (perform-replace "<user>" (user-full-name) nil nil nil )
  (goto-char (point-min))
  (perform-replace "<class>" name nil nil nil )
  (goto-char (point-min))
  (perform-replace "<date>"  (current-time-string) nil nil nil )
  (goto-char (point-min))
  (perform-replace "<start>" "" nil nil nil ))

;;---------------------------------------------------------------------------;;
;; set up a new C header.in file in draco

(defun draco-c-headin (name)
  "Function to set up a C header.in file"

  (interactive "sC Header.in Name: ")

  ;; get package directory
  (setq pkg (draco-pkg-header))
  
  ;; now get the safe pkg name
  (setq spkg (draco-safe-pkg pkg))

  (setq tempdir (rtt-templates))
  (setq nfile (concat name ".h.in"))
  (setq tfile (concat tempdir "/template.h"))
  (if (file-exists-p nfile)
      (error "File %s already exists." nfile))
  
  (find-file nfile)
  (insert-file-contents tfile)

  ;; do Furnish style customizations

  (find-file nfile)
  (perform-replace "<pkg>" pkg nil nil nil )
  (goto-char (point-min))
  (perform-replace "<spkg>" spkg nil nil nil )
  (goto-char (point-min))
  (perform-replace "<user>" (user-full-name) nil nil nil )
  (goto-char (point-min))
  (perform-replace "<class>" name nil nil nil )
  (goto-char (point-min))
  (perform-replace "<date>"  (current-time-string) nil nil nil )
  (goto-char (point-min))
  (perform-replace "<start>" "" nil nil nil ))

;;---------------------------------------------------------------------------;;
;; set up a new C++ implementation file(s) [.cc,.t.hh] in draco

(defun draco-cc-imp (name)
  "Function to set up C++ implementation files [.cc,.t.hh] in draco"

  (interactive "sC++ Base Class Name: ")

  ;; get the templates directory
  (setq tempdir (rtt-templates))

  ;; get the package directory
  (setq pkg (draco-pkg-header))

  ;; now get the safe pkg name
  (setq spkg (draco-safe-pkg pkg))

  ;; now get the namespace
  (setq nspace (draco-namespace spkg))

  ;;
  ;; get the proper template for the implementation file (.cc)
  ;;

  (setq nfile (concat name ".cc"))
  (setq tfile (concat tempdir "/template.cc"))
  (if (file-exists-p nfile)
      (error "File %s already exists." nfile))

  (find-file nfile)
  (insert-file-contents tfile)

  ;; insert the goodies
  
  (find-file nfile)
  (perform-replace "<user>" (user-full-name) nil nil nil )
  (goto-char (point-min))
  (perform-replace "<pkg>" pkg nil nil nil )
  (goto-char (point-min))
  (perform-replace "<spkg>" spkg nil nil nil )
  (goto-char (point-min))
  (perform-replace "<namespace>" nspace nil nil nil )
  (goto-char (point-min))
  (perform-replace "<class>" name nil nil nil )
  (goto-char (point-min))
  (perform-replace "<basename>" name nil nil nil )
  (goto-char (point-min))
  (perform-replace "<date>"  (current-time-string) nil nil nil )
  (goto-char (point-min))
  (perform-replace "<start>" "" nil nil nil )

  ;;
  ;; get the proper template for the implementation file (.t.hh)
  ;;

  (setq nfile (concat name ".t.hh"))
  (setq tfile (concat tempdir "/template.t.hh"))
  (if (file-exists-p nfile)
      (error "File %s already exists." nfile))

  (find-file nfile)
  (insert-file-contents tfile)

  ;; insert the goodies
  
  (find-file nfile)
  (perform-replace "<user>" (user-full-name) nil nil nil )
  (goto-char (point-min))
  (perform-replace "<class>" name nil nil nil )
  (goto-char (point-min))
  (perform-replace "<pkg>" pkg nil nil nil )
  (goto-char (point-min))
  (perform-replace "<spkg>" spkg nil nil nil )
  (goto-char (point-min))
  (perform-replace "<namespace>" nspace nil nil nil )
  (goto-char (point-min))
  (perform-replace "<basename>" name nil nil nil )
  (goto-char (point-min))
  (perform-replace "<date>"  (current-time-string) nil nil nil )
  (goto-char (point-min))
  (perform-replace "<start>" "" nil nil nil ))

;;---------------------------------------------------------------------------;;
;; set up a new C++ instantiation file(s) [_pt.cc] in draco

(defun draco-cc-pt (name)
  "Function to set up C++ instantiation files [_pt.cc] in draco"

  (interactive "sC++ Instantiation Name: ")

  ;; get the templates directory
  (setq tempdir (rtt-templates))

  ;; get the package directory
  (setq pkg (draco-pkg-header))

  ;; now get the safe pkg name
  (setq spkg (draco-safe-pkg pkg))

  ;; now get the namespace
  (setq nspace (draco-namespace spkg))

  ;;
  ;; get the proper template for the instantiation file (.cc)
  ;;

  (setq name (concat name "_pt"))
  (setq nfile (concat name ".cc"))
  (setq tfile (concat tempdir "/template.cc"))
  (if (file-exists-p nfile)
      (error "File %s already exists." nfile))

  (find-file nfile)
  (insert-file-contents tfile)

  ;; insert the goodies
  
  (find-file nfile)
  (perform-replace "<user>" (user-full-name) nil nil nil )
  (goto-char (point-min))
  (perform-replace "<class>" name nil nil nil )
  (goto-char (point-min))
  (perform-replace "<basename>" name nil nil nil )
  (goto-char (point-min))
  (perform-replace "<pkg>" pkg nil nil nil )
  (goto-char (point-min))
  (perform-replace "<spkg>" spkg nil nil nil )
  (goto-char (point-min))
  (perform-replace "<namespace>" nspace nil nil nil )
  (goto-char (point-min))
  (perform-replace "<date>"  (current-time-string) nil nil nil )
  (goto-char (point-min))
  (perform-replace "<start>" "" nil nil nil ))

;;---------------------------------------------------------------------------;;
;; set up a draco python file

(defun draco-python (name)
  "Function to spontaneously setup a new Python file in draco"
  
  (interactive "sPython name: ")
  
  (setq tempdir (rtt-templates))
  (setq nfile (concat name ".py"))
  (setq tfile (concat tempdir "/template.py"))
  (if (file-exists-p nfile)
      (error "File already exists." nfile))
  
  (find-file nfile)
  (insert-file-contents tfile)

  ;; do Furnish style customizations

  (find-file nfile)
  (perform-replace "<user>" (user-full-name) nil nil nil )
  (goto-char (point-min))
  (perform-replace "<basename>" name nil nil nil )
  (goto-char (point-min))
  (perform-replace "<date>"  (current-time-string) nil nil nil)
  (goto-char (point-min))
  (perform-replace "<start>" "" nil nil nil ))

;;---------------------------------------------------------------------------;;
;; set up a new misc. Makefile.in-type draco

(defun draco-make (name)
  "Function to set up a Makefile.type file"

  (interactive "sMakefile: ")

  (setq type (read-from-minibuffer "Type: " "in"))

  (setq tempdir (rtt-templates))
  (setq nfile (concat  "Makefile." type))
  (setq tfile (concat tempdir "/Makefile.temp"))
  (if (file-exists-p nfile)
      (error "File %s already exists." nfile))
  
  (find-file nfile)
  (insert-file-contents tfile)

  ;; do Furnish style customizations

  (find-file nfile)
  (perform-replace "<name>" name nil nil nil )
  (goto-char (point-min))
  (perform-replace "<type>" type nil nil nil )
  (goto-char (point-min))
  (perform-replace "<user>" (user-full-name) nil nil nil )
  (goto-char (point-min))
  (perform-replace "<date>"  (current-time-string) nil nil nil )
  (goto-char (point-min))
  (perform-replace "<start>" "" nil nil nil ))

;;---------------------------------------------------------------------------;;
;; set up a new C++ test executable in draco

(defun draco-cc-test (name)
  "Function to set up C++ test executable in draco"

  (interactive "sC++ Test Executable Name: ")

  ;; get the templates directory
  (setq tempdir (rtt-templates))

  ;; get the package directory
  (setq pkg (draco-pkg-in-test-header))

  ;; get test package
  (setq tpkg (draco-pkg-test-header))

  ;; now get the safe pkg name
  (setq spkg (draco-safe-pkg tpkg))

  ;; now get the namespace
  (setq nspace (draco-namespace spkg))

  ;; determine if this is parallel or not
  (setq parallel (read-from-minibuffer "Parallel: " "y"))

  (setq tfile (concat tempdir "/template_test.cc"))
  (if (string= parallel "y")
      (setq tfile (concat tempdir "/template_c4_test.cc")))

  ;;
  ;; get the proper template for the test executable file (.cc)
  ;;

  (setq nfile (concat name ".cc"))
  (if (file-exists-p nfile)
      (error "File %s already exists." nfile))
  
  (find-file nfile)
  (insert-file-contents tfile)

  ;; insert the goodies
 
  (find-file nfile)
  (perform-replace "<user>" (user-full-name) nil nil nil )
  (goto-char (point-min))
  (perform-replace "<pkg>" pkg nil nil nil )
  (goto-char (point-min))
  (perform-replace "<tpkg>" tpkg nil nil nil )
  (goto-char (point-min))
  (perform-replace "<spkg>" spkg nil nil nil )
  (goto-char (point-min))
  (perform-replace "<namespace>" nspace nil nil nil )
  (goto-char (point-min))
  (perform-replace "<class>" name nil nil nil )
  (goto-char (point-min))
  (perform-replace "<date>"  (current-time-string) nil nil nil )
  (goto-char (point-min))
  (perform-replace "<start>" "" nil nil nil ))

;;---------------------------------------------------------------------------;;
;; LATEX TEMPLATES USEFUL IN DRACO
;;---------------------------------------------------------------------------;;
;; 1) Memo                               [draco-memo]
;; 2) Research Note                      [draco-note]
;; 3) article or LA-UR of small size     [draco-article]
;; 4) LA Report or LA=UR of large size   [draco-report]
;; 5) BiBTeX files                       [draco-bib]
;; 6) Project Vision/Scope Statement     [draco-viscope] 

;;---------------------------------------------------------------------------;;
;; set up a memo

(defun draco-memo (name)
  "Function to spontaneously setup a new LaTeX LANL style memo"

  (interactive "sMemo Name: ")

  (setq tempdir (rtt-templates))
  (setq nfile (concat name ".tex"))
  (setq tfile (concat tempdir "/draco_memo.tex"))
  (if (file-exists-p nfile)
      (error "File %s already exists." nfile))

  (find-file nfile)
  (insert-file-contents tfile)

  ;; insert the goodies

  (find-file nfile)
  (perform-replace "<user>" (user-full-name) nil nil nil )
  (goto-char (point-min))
  (perform-replace "<papername>" name nil nil nil )
  (goto-char (point-min))
  (perform-replace "<start>" "" nil nil nil ))

;;---------------------------------------------------------------------------;;
;; set up a research note

(defun draco-note (name)
  "Function to spontaneously setup a new LaTeX LANL style research note"

  (interactive "sResearch Note Name: ")

  (setq tempdir (rtt-templates))
  (setq nfile (concat name ".tex"))
  (setq tfile (concat tempdir "/draco_note.tex"))
  (if (file-exists-p nfile)
      (error "File %s already exists." nfile))

  (find-file nfile)
  (insert-file-contents tfile)

  ;; do Furnish style customizations

  (find-file nfile)
  (perform-replace "<user>" (user-full-name) nil nil nil )
  (goto-char (point-min))
  (perform-replace "<papername>" name nil nil nil )
  (goto-char (point-min))
  (perform-replace "<start>" "" nil nil nil ))

;;---------------------------------------------------------------------------;;
;; set up an article based on the LANL style of TME

(defun draco-article (name)
  "Function to spontaneously setup a new LaTeX LANL style article"

  (interactive "sPaper Name: ")

  (setq tempdir (rtt-templates))
  (setq nfile (concat name ".tex"))
  (setq tfile (concat tempdir "/draco_art.tex"))
  (if (file-exists-p nfile)
      (error "File %s already exists." nfile))

  (find-file nfile)
  (insert-file-contents tfile)

  ;; do Furnish style customizations

  (find-file nfile)
  (perform-replace "<user>" (user-full-name) nil nil nil )
  (goto-char (point-min))
  (perform-replace "<papername>" name nil nil nil )
  (goto-char (point-min))
  (perform-replace "<start>" "" nil nil nil ))

;;---------------------------------------------------------------------------;;
;; set up a draco report

(defun draco-report (name)
  "Function to spontaneously setup a new LaTeX LANL style report"

  (interactive "sReport Name: ")

  (setq tempdir (rtt-templates))
  (setq nfile (concat name ".tex"))
  (setq tfile (concat tempdir "/draco_rep.tex"))
  (if (file-exists-p nfile)
      (error "File %s already exists." nfile))

  (find-file nfile)
  (insert-file-contents tfile)

  ;; do Furnish style customizations

  (find-file nfile)
  (perform-replace "<user>" (user-full-name) nil nil nil )
  (goto-char (point-min))
  (perform-replace "<papername>" name nil nil nil )
  (goto-char (point-min))
  (perform-replace "<start>" "" nil nil nil ))

;;---------------------------------------------------------------------------;;
;; set up a new bib file

(defun draco-bib (name)
  "Function to spontaneously setup a new LaTeX BiBTeX file"

  (interactive "sBib Name: ")

  (setq tempdir (rtt-templates))
  (setq nfile (concat name ".bib"))
  (setq tfile (concat tempdir "draco_bib.bib"))
  (if (file-exists-p nfile)
      (error "File %s already exists." nfile))

  (find-file nfile)
  (insert-file-contents tfile)

  ;; do Furnish style customizations

  (find-file nfile)
  (perform-replace "<user>" (user-full-name) nil nil nil )
  (goto-char (point-min))
  (perform-replace "<papername>" name nil nil nil )
  (goto-char (point-min))
  (perform-replace "<start>" "" nil nil nil ))

;;---------------------------------------------------------------------------;;
;; set up a "project vision/scope statement" memo 

(defun draco-viscope (name)
  "Function to spontaneously setup a new vision/scope memo"

  (interactive "sVision/Scope Memo Name: ")

  (setq tempdir (rtt-templates))
  (setq nfile (concat name ".tex"))
  (setq tfile (concat tempdir "/draco_viscope.tex"))
  (if (file-exists-p nfile)
      (error "File %s already exists." nfile))

  (find-file nfile)
  (insert-file-contents tfile)

  ;; insert the goodies

  (find-file nfile)
  (perform-replace "<user>" (user-full-name) nil nil nil )
  (goto-char (point-min))
  (perform-replace "<papername>" name nil nil nil )
  (goto-char (point-min))
  (perform-replace "<start>" "" nil nil nil ))

;;---------------------------------------------------------------------------;;
;; set up draco-hacks
;;---------------------------------------------------------------------------;;

(defun draco-setup ()
  (message "Configuring DRACO support."))

(provide 'draco-hacks)

;;---------------------------------------------------------------------------;;
;; end of draco-hacks.el
;;---------------------------------------------------------------------------;;


