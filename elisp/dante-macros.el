;;; dante-macros.el
;;; Kelly Thompson
;;; June 4, 2001

;;; $Id$

;;---------------------------------------------------------------------------;;
;; provide macros for the Dante development environment, we do not
;; provide key bindings here.
;;---------------------------------------------------------------------------;;

;; We borrow some elisp code developed by Tom Evans that is found in
;; draco-hack.el. 
(require 'draco-hacks)

;;---------------------------------------------------------------------------;;
;; SERVICE FUNCTIONS USED IN DANTE-MACROS.EL
;;---------------------------------------------------------------------------;;

(defun create-buffer-with-file (newfilename templatename packagename
					    safepackagename classname)
  "
*Create a new buffer named newfilename that contains a copy of the
contents found in the file templatename (a fully qualified path).  The
template file may contain any of the following strings which will be
expanded when the new buffer is created:

  <user> - replaced by (user-full-name).
  <date> - replaced by (current-time-string).
  <pkg>  - replaced by the function argument packagename.
  <spkg> - replaced by the function argument safepackagename.
"
  (interactive)
  (progn
    (cond ((not (file-exists-p newfilename))
	   (if (file-exists-p newfilename)
	       (error "File" newfilename " already exists." newfilename))
	   
	   ; create a new buffer and insert the template
	   (find-file newfilename)
	   (insert-file-contents templatename)
	   
	   ;; insert the goodies into the new configure.in
	   (find-file newfilename)
	   (perform-replace "<user>" (user-full-name) nil nil nil )
	   (goto-char (point-min))
	   (perform-replace "<pkg>" packagename nil nil nil )
	   (goto-char (point-min))
	   (perform-replace "<class>" classname nil nil nil )
	   (goto-char (point-min))
	   (perform-replace "<spkg>" safepackagename nil nil nil )
	   (goto-char (point-min))
	   (perform-replace "<date>"  (current-time-string) nil nil nil)
	   (goto-char (point-min))
	   (perform-replace "<start>" "" nil nil nil )))
    ))

;;---------------------------------------------------------------------------;;
;; Query for the class name, with guessed default = package name.

(defun dante-class-name (pkg)
  "Function to get a class name.  The default value will be the
package (i.e. directory) name."
  (read-from-minibuffer "Class: " pkg)
)


;;---------------------------------------------------------------------------;;
;; DANTE ENVIRONMENT FUNCTIONS (INTERACTIVE)
;;---------------------------------------------------------------------------;;
;;  1) setting up a package                          [dante-shadow-package]


;;---------------------------------------------------------------------------;;

(defun dante-shadow-package ()
  "
*Function to set up a Dante Shadow Object package directory with stuff.

The files that we will place into all package directories are:
  
  Build: configure.in, Sources.mk, config.h.in,
       and static_depend.in
  C++: package.cc, package.hh and package_flat.hh
  F90: package_flat.F, package_assert.F, package_opq.F, package_class.F,
       libpackage.F, package_realkind.F and package_class.F
  test: Create test directory and put an empty static_depend.in in it.
  doc: Create an empty doc directory.
"
  (interactive)

  ;; these files are based on templates in the
  ;; danteV2/DanteCXXShadows/templates directory that can be accessed
  ;; using the rtt-templates() elisp function

  ;; first find templates directory
  (setq tempdir (rtt-templates))

  ;; now get the header for this pkg
  (setq pkg (draco-pkg-header))

  ;; now get the safe pkg name
  (setq spkg (draco-safe-pkg pkg))

  ;; Request the new class name (default to the package name)
  (setq class-name (dante-class-name pkg))

  ;;
  ;; first lets get the configure.in file
  ;;
  (create-buffer-with-file "configure.in"
			   (concat tempdir "/configure.package.in")
			   pkg spkg class-name)
  ;; 
  ;; now lets get the sources.mk file
  ;;
  (create-buffer-with-file "Sources.mk"
			   (concat tempdir "/sources.package.mk")
			   pkg spkg class-name)
  ;;
  ;; now lets get the config.h.in file
  ;;
  (create-buffer-with-file "config.h.in"
			   (concat tempdir "/config.h.in")
			   pkg spkg class-name)
  ;;
  ;; now lets get the site.h file
  ;;
  ; (create-buffer-with-file "site.h"
  ;			   (concat tempdir "/site.package.h")
  ;			   pkg spkg class-name)
  ;;
  ;; now lets get the siteF.h file
  ;;
  ; (create-buffer-with-file "siteF.h"
  ;			   (concat tempdir "/siteF.package.h")
  ;			   pkg spkg class-name)
  ;;
  ;; now lets get the static_depend.in file
  ;;
  (create-buffer-with-file "static_depend.in"
			   (concat tempdir "/static_depend.package.in")
			   pkg spkg class-name)
  ;;
  ;; first lets get the class.cc file
  ;;
  (create-buffer-with-file (concat class-name ".cc")
                           (concat tempdir "/template.cc")
                           pkg spkg class-name)
  ;; 
  ;; now lets get the class.hh file
  ;;
  (create-buffer-with-file (concat class-name ".hh")
			   (concat tempdir "/template.hh")
			   pkg spkg class-name)
  ;; 
  ;; now lets get the class_flat.hh file
  ;;
  (create-buffer-with-file (concat class-name "_flat.hh")
			   (concat tempdir "/template_flat.hh")
			   pkg spkg class-name)
  ;; 
  ;; now lets get the class_flat.F file
  ;;
  (create-buffer-with-file (concat class-name "_flat.F")
			   (concat tempdir "/template_flat.F")
			   pkg spkg class-name)
  ;; 
  ;; now lets get the pkg_assert.F file
  ;;
  (create-buffer-with-file (concat spkg "_assert.F")
			   (concat tempdir "/template_assert.F")
			   pkg spkg class-name)
  ;; 
  ;; now lets get the class_opq.F file
  ;;
  (create-buffer-with-file (concat class-name "_opq.F")
			   (concat tempdir "/template_opq.F")
			   pkg spkg class-name)
  ;; 
  ;; now lets get the pkg_class.F file
  ;;
  (create-buffer-with-file (concat spkg "_class.F")
			   (concat tempdir "/template_class.F")
			   pkg spkg class-name)
  ;; 
  ;; now lets get the libpkg.F file
  ;;
  (create-buffer-with-file (concat "lib" spkg ".F")
			   (concat tempdir "/libtemplate.F")
			   pkg spkg class-name)
  ;; 
  ;; now lets get the libpkg.F file
  ;;
  (create-buffer-with-file (concat spkg "_realkind.F")
			   (concat tempdir "/template_realkind.F")
			   pkg spkg class-name)
  ;; 
  ;; now create the test directory and put an empty static_depend file
  ;; in it.
  ;;
  (if (not (file-accessible-directory-p "test"))
      (make-directory "test"))
  (if (not (file-exists-p "test/static_depend.in"))
      (shell-command "touch test/static_depend.in"))
  
  ;;
  ;; Create an empty doc directory if one does not already exist.
  ;;
  (if (not (file-accessible-directory-p "doc"))
      (make-directory "doc"))
  
  ) ; end dante-shadow-package

;;---------------------------------------------------------------------------;;
;; set up dante-macros
;;---------------------------------------------------------------------------;;

(defun dante-setup ()
  (message "Configuring Dante support."))

(provide 'dante-macros)

;;---------------------------------------------------------------------------;;
;; end of dante-macros.el
;;---------------------------------------------------------------------------;;


