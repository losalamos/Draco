;;; rtt-menu.el
;;; Kelly (K.T.) Thompson
;;; May 25, 2001

;;; $Id$

;;---------------------------------------------------------------------------;;
;; Provide menu access to many emacs macros written specifically for RTT
;;---------------------------------------------------------------------------;;

; We need the rtt-hacks package to be loaded so that the "New files"
; section works.
(require 'draco-hacks)

; We need the Config-pkg package to be loaded so that the "Modes"
; section works.
(require 'Config-pkg)

; The function dante-shadow-package is defined in dante-macros.el
(require 'dante-macros)

; Provide a way for individuals to indicate if they want to use the
; RTT menu.  Presumably, people who spend most of their time working
; on one of Draco like codes would have the following in their .emacs:
;     (setq want-rtt-menu t)
; The default value will be nill
(defvar want-rtt-menu nil
  "Deos the user wnat to use the RTT pull down menu?")

;;---------------------------------------------------------------------------;;
;; RTT MENU
;;---------------------------------------------------------------------------;;

(defun rtt-full-menu-setup ()
  (message "Configuring RTT Menu (rtt-full-menu-setup)")

      (let* ((current-menubar (default-value 'current-menubar))       
	     (edit-menu (assoc "Edit" current-menubar)))
	
	(unless (assoc "Transpose" edit-menu) ; Idempotence
	  
	  ;; Rtt Submenu
	  (add-submenu
	   nil
	   '("RTT"
	     ("New Files"
	      ["New C++ class"     draco-class t]
	      ["New C++ header"    draco-cc-head t]
	      ["New C++ header.in" draco-cc-headin t]
	      ["New C header"      draco-c-head t]
	      ["New C header.in"   draco-c-headin t]
	      ["New C++ template impl. file (.t.hh)" draco-cc-imp t]
	      ["New C++ instantiation file (_pt.cc)" draco-cc-pt t]
	      ["New Python file"   draco-python t]
	      ["New specialized makefile" draco-make t]
;	      ("Draco"
	       ["New C++ package" draco-package t]
	       ["New C++ package/test" draco-package-test t]
               ["New C++ package/autodoc" draco-package-doc t]
	       ["New C++ unit test executable" draco-cc-test t]
;	       )
	      ("Dante"
	       ["New Dante shadow package" dante-shadow-package t]
	       )
	      )
;	     "-----"
	     ("Editing RTT files"
	      ["Find companion file"        rtt-find-companion-file t]
	      ["Insert C++ comment block"   rtt-insert-class-doc t]
	      ["Insert C++ comment divider" rtt-insert-comment-divider t]
	      ["Insert C++ function comment divider" rtt-insert-function-doc t]
	      "---"
	      ["Insert a Makefile divider"  tme-makefile-comment-divider t]
	      "---"
	      ["Insert ChangeLog Record"  add-change-log-entry t]
	      )
;	     "-----"
	     ("Documents"
	      ["Create a memo"                 draco-memo t]
	      ["Create a note"                 draco-note t]
	      ["Create an article"             draco-article t]
	      ["Create a report"               draco-report t]
	      ["Create a BiB-TeX File"         draco-bib t]
	      ["Create a Vision/Scope Memo"    draco-viscope t]
	      ["Create a Bug Post-Mortem Memo" draco-bug-pm t]
	      )
;	     "---"
	; No longer using GNATS bug reporting stuff
 
	;       ("Problem Reports"
	;	["Report a problem"        send-pr t]
	;	["Edit a problem report"   edit-pr t]
	;	["Query a problem report"  query-pr t]
	;	)
	     ("Modes"
	      
	      ["MPPL"     (progn
			    (if (not want-mppl-mode) 
				(progn
				  (setq want-mppl-mode t)
				  (draco-init-mppl-mode))
			      (message "MPPL mode already configured")) 
			    (mppl-mode)
			    (font-lock-fontify-buffer))
	       :style toggle
	       :selected want-mppl-mode]
	      
	      ["TCL"       (progn
			     (if (not want-tcl-mode)
				 (progn
				   (setq want-tcl-mode t)
				   (draco-init-tcl-mode))
			       (message "TCL mode already configured")) 
			     (tcl-mode)
			     (font-lock-fontify-buffer))
	       :style toggle
	       :selected want-tcl-mode]
	      
	      ["Python"    (progn
			     (if (not want-python-mode)
				 (progn
				   (setq want-python-mode t)
				   (draco-init-python-mode))
			       (message "Python mode already configured"))
			     (python-mode)
			     (font-lock-fontify-buffer))
	       :style toggle
	       :selected want-python-mode]
	      
	      ["Makefile"  (progn
			     (if (not want-makefile-mode)
				 (progn
				   (setq want-makefile-mode t)
				   (draco-init-makefile-mode))
			       (message "Makefile mode already configured"))
			     (makefile-mode)
			     (font-lock-fontify-buffer))
	       :style toggle
	       :selected want-makefile-mode]
	      
	      ["C/C++"  (progn
			  (if (not want-cc-mode)
			      (progn
				(setq want-cc-mode t)
				(draco-init-cc-mode))
			    (message "C/C++ mode already configured"))
			  (cc++-mode)
			  (font-lock-fontify-buffer))
	       :style toggle
	       :selected want-cc-mode]
	      
	      ["AucTeX" (if (not want-auctex-pkg)
			    (progn
			      (setq want-auctex-pkg t)
			      (draco-init-auctex-mode))
			  (message "AucTeX mode already configured"))
	       :style toggle
	       :selected want-auctex-pkg]
	      
	      ["F90"    (progn
			  (if (not want-f90-mode)
			      (progn
				(setq want-f90-mode t)
				(draco-init-f90-mode))
			    (message "F90 mode already configured"))
			  (f90-mode)
			  (font-lock-fontify-buffer))
	       :style toggle
	       :selected want-f90-mode]
	      
	      ["FORTAN" (progn
			  (if (not want-fortran-mode)
			      (progn
				(setq want-fortran-mode t)
				(draco-init-fortran-mode))
			    (message "FORTRAN mode already configured"))
			  (fortran-mode)
			  (font-lock-fontify-buffer))
	       :style toggle
	       :selected want-fortran-mode]
	      
	      ["ChangeLog" (progn
			  (if (not want-change-log-mode)
			      (progn
				(setq want-change-log-mode t)
				(draco-init-change-log-mode))
			    (message "CHANGE-LOG mode already configured"))
			  (change-log-mode)
			  (font-lock-fontify-buffer))
	       :style toggle
	       :selected want-change-log-mode]
	      
	      )
; Menu to get into pcl-cvs-mode (see
; http://www.lns.cornell.edu/public/COMP/info/xemacs/pcl-cvs/pcl-cvs_toc.html) 
	   ("XEmacs Features"
	    ["Speedbar"            speedbar-frame-mode               t]
            ["CVS Add ChangeLog"   add-change-log-entry-other-window t]
	    ["CVS Examine"         cvs-examine                       t]
	    ["CVS Status"          cvs-status                        t]
	    ["CVS Checkout"        cvs-checkout                      t]
	    )
	   )

	   "Draco menubar")))
    ) ; end rtt-menu-setup()  

;;---------------------------------------------------------------------------;;
;; set up rtt-menu
;;---------------------------------------------------------------------------;;

(defun rtt-menu-setup ()
  (rtt-full-menu-setup)
  (message "Configuring RTT menu."))

(provide 'rtt-menu)

;;---------------------------------------------------------------------------;;
;; end of rtt-menu.el
;;---------------------------------------------------------------------------;;


;; Notes

; See
; usr/lib/xemacs/xemacs-packages/lisp/view-process/view-process-xemacs.el
; for more examples on using radio buttons, etc in menus.
