;; Config-pkg.el
;; Geoffrey Furnish, Kelly Thompson
;; 23 June 1997
;;
;; $Id$
;;
;; Configure a variety of packages, upon request of user.

;; Put the following in your .emacs or use ccs4defaults.el

; (setq config-pkg-verbose t
;      config-pkg-colorize-modeline t
;      want-mppl-mode     nil
;      want-tcl-mode      t
;      want-python-mode   t
;      want-makefile-mode t
;      want-cc-mode       t
;      want-auctex-pkg    t
;      want-f90-mode      t
;      want-emacs-lisp-mode t
;      want-shell-mode    t
;      want-cvs-mode      t
;      want-doxymacs-mode t
;      want-autoconf-mode t
;      want-sgml-mode     t
;      want-fortran-mode  t)
; (load-library "Config-pkg")

(defvar config-pkg-verbose nil
  "*Does the user want us to be chatty?")

(defvar config-pkg-colorize-modeline nil
  "*Does the user want us to colorize the modeline 
based on the buffer-mode?  When customizing these colors 
it may be useful to run the XEmacs command 

M-x list-colors-display 

to obtain a list of colors known to XEmacs.")

;; ========================================
;; MPPL
;; ========================================

(defvar want-mppl-mode nil
  "*Does the user want to have MPPL mode?")

(defun draco-init-mppl-mode ()
  "Autoload mppl-mode, append appropriate suffixes to auto-mode-alist
and add turn-on-auto-fill to the mppl-mode-hook."
  (interactive)
  (progn
    (if config-pkg-verbose
	(message "configuring MPPL mode"))
    (autoload 'mppl-mode   "mppl-mode" nil t)
    (setq auto-mode-alist
	  (append '(("\\.i$" . mppl-mode)
		    ("\\.p$" . mppl-mode)
		    ("\\.m$" . mppl-mode)
		    ("\\.pm4$" . mppl-mode)
		    ) auto-mode-alist))
    (add-hook 'mppl-mode-hook 'turn-on-auto-fill)))

(if want-mppl-mode (draco-init-mppl-mode))

;; ========================================
;; TCL
;; ========================================

(defvar want-tcl-mode nil
  "*Does the user want to have TCL mode?")

(defun draco-init-tcl-mode ()
  "Autoload tcl-mode and append the appropriate suffixes to 
auto-mode-alist."
  (interactive)
  (progn
    (if config-pkg-verbose
	(message "Configuring TCL mode"))
    (autoload 'tcl-mode   "tcl-mode" nil t)
    (setq auto-mode-alist
	  (append '(("\\.tcl$" . tcl-mode)
		    ("\\.itk$" . tcl-mode)
		    ("\\.ith$" . tcl-mode)
		    ("\\.itm$" . tcl-mode)
		    ) auto-mode-alist))
    (add-hook 'tcl-mode-hook 'turn-on-auto-fill)))

(if want-tcl-mode (draco-init-tcl-mode))

;; ========================================
;; Python
;; ========================================

(defvar want-python-mode nil
  "*Does the user want to have Python mode?")

(defun draco-init-python-mode ()
  "Autoload python-mode and append the appropriate suffixes to
auto-mode-alist."
  (interactive)
      (progn
	(if config-pkg-verbose
	    (message "Configuring Python mode"))
      (autoload 'python-mode "python-mode" "Python editing mode." t)
      (setq auto-mode-alist
	    (cons '("\\.py$" . python-mode) auto-mode-alist))
      (defun rtt-python-mode-hook ()
	"RTT hooks added to Python mode."
	(defvar py-indent-offset 4)
	(local-set-key [(f5)] 'tme-makefile-divider)
	(local-set-key [(f6)] 'tme-makefile-comment-divider))
      (add-hook 'python-mode-hook 'rtt-python-mode-hook)
      (add-hook 'python-mode-hook 'turn-on-auto-fill)))

(if want-python-mode (draco-init-python-mode))

;; ========================================
;; Autoconf
;; ========================================

(defvar want-autoconf-mode nil
  "*Does the user want to have autoconf mode?")

(defun draco-init-autoconf-mode ()
  "Autoload autoconf-mode and append the appropriate suffixes to
auto-mode-alist."
  (interactive)
      (progn
	(if config-pkg-verbose
	    (message "Configuring Autoconf mode"))
      (setq auto-mode-alist
	    (append '(("\\.m4$" . autoconf-mode)
		      ("\\.ac$" . autoconf-mode)
		      ("\\.in$" . autoconf-mode))
		    auto-mode-alist))
      (defun rtt-autoconf-mode-hook ()
	"RTT hooks added to autoconf mode"
	(setq tab-stop-list '(3 7 11 15 19 23 27 31 35 39 43 47 51 55 59 63 67 71 75 79 83))
	(local-set-key [(f5)] 'tme-makefile-divider)
	(local-set-key [(f6)] 'tme-makefile-comment-divider)
	(local-set-key [(f9)] 'tme-m4-divider)
	(local-set-key [(tab)] 'tab-to-tab-stop)
	(local-set-key [(f10)] 'tme-m4-comment-divider))
      (add-hook 'autoconf-mode-hook 'rtt-autoconf-mode-hook)
      (add-hook 'autoconf-mode-hook 'turn-on-auto-fill)))

(if want-autoconf-mode (draco-init-autoconf-mode))

;; ========================================
;; Makefile
;; ========================================

(defvar want-makefile-mode nil
  "*Does the user want to have Makefile mode?")

(defun draco-init-makefile-mode ()
  "Autoload makefile-mode and append the appropriate suffixes to
auto-mode-alist."
  (interactive)
    (progn
      (if config-pkg-verbose
	  (message "Configuring Makefile mode"))
      (autoload 'makefile-mode   "make-mode" nil t)
      (setq auto-mode-alist
	    (append '(("makefile" . makefile-mode)
		      ("Makefile" . makefile-mode)
		      ) auto-mode-alist))
      (defun rtt-makefile-mode-hook ()
	"RTT hooks added to Makefile mode."
	(local-set-key [(f5)] 'tme-makefile-divider)
	(local-set-key [(f6)] 'tme-makefile-comment-divider))
      (add-hook 'makefile-mode-hook 'rtt-makefile-mode-hook t)
      (add-hook 'makefile-mode-hook 'turn-on-font-lock)
      (add-hook 'makefile-mode-hook 'turn-on-auto-fill t)))

(if want-makefile-mode (draco-init-makefile-mode))


;; ========================================
;; C++
;; ========================================

(defvar want-cc-mode nil
  "*Does the user want to have C/C++ mode?")

(defun draco-init-cc-mode ()
  "Autoload c++-mode, c-mode and append the appropriate suffixes to 
auto-mode-alist."
  (interactive)
    (progn
      (if config-pkg-verbose
	  (message "Configuring C/C++ mode"))
      (autoload 'c++-mode "cc-mode" "C++ Editing Mode" t)
      (autoload 'c-mode   "cc-mode" "C Editing Mode" t)
      (if config-pkg-colorize-modeline 
	  (add-hook 'c++-mode-hook        
		    '(lambda () 
		       (set-face-background 'modeline 
					    "skyblue" (current-buffer))
		       (set-face-foreground 'modeline 
					    "black"   (current-buffer)))))
      (if config-pkg-colorize-modeline 
	  (add-hook 'c-mode-hook        
		    '(lambda () 
		       (set-face-background 'modeline 
					    "pink" (current-buffer))
		       (set-face-foreground 'modeline 
					    "black"   (current-buffer)))))
      (setq auto-mode-alist
	    (append '(("\\.C$"   . c++-mode)
		      ("\\.cc$"  . c++-mode)
                      ("\\.pt$"  . c++-mode)
		      ("\\.hh$"  . c++-mode)
		      ("\\.hpp$"  . c++-mode)
		      ("\\.cpp$"  . c++-mode)
                      ("\\.hh.in$" . c++-mode)
		      ("\\.h.in$"  . c-mode)
		      ("\\.c$"   . c-mode)   ; to edit C code
		      ("\\.h$"   . c-mode)   ; to edit C code
 		      ("\\.dot$" . c-mode)  ; for dot files
		      ) auto-mode-alist))
	(require 'fl-keywords)
	(add-hook 'c-mode-common-hook 'add-draco-dbc-font-lock-keywords)
      (add-hook 'c++-mode-hook 'turn-on-font-lock)
      (add-hook 'c++-mode-hook 'turn-on-auto-fill)
      (add-hook 'c-mode-hook 'turn-on-font-lock)
      (add-hook 'c-mode-hook 'turn-on-auto-fill)))

(if want-cc-mode (draco-init-cc-mode))

;; ========================================
;; AUCTEX
;; ========================================

(defvar want-auctex-pkg nil
  "*Does the user want to have AucTeX mode?")

(defun draco-init-auctex-mode ()
  "Loads the tex-site package."
  (interactive)
    (progn
      (if config-pkg-verbose
	  (message "Configuring AucTeX mode."))
      (require 'tex-site)
      (setq auto-mode-alist
	    (append '(("\\.tex$" . tex-mode)
		      ("\\.bib$" . bibtex-mode)
                      ("\\.bst$" . tex-mode)
		      ("\\.bbl$" . tex-mode)
                      ("\\.blg$" . tex-mode)
		      ("\\.idx$" . tex-mode)
		      ("\\.ilg$" . tex-mode)   
		      ("\\.ind$" . tex-mode)  
 		      ("\\.toc$" . tex-mode)
		      ("\\.lof$" . tex-mode)
		      ("\\.lot$" . tex-mode)
		      ("\\.cls$" . tex-mode)
		      ("\\.sty$" . tex-mode)            
		      ) auto-mode-alist))
      (add-hook 'TeX-mode-hook 'turn-on-auto-fill)))

(if want-auctex-pkg (draco-init-auctex-mode))

;; ========================================
;; FORTRAN-90
;; ========================================

(defvar want-f90-mode nil
  "*Does the user want to have f90-mode?")

(defun draco-init-f90-mode ()
  "Autoload f90-mode and append the approriate suffixes to
auto-mode-alist."
  (interactive)
  (progn
    (if config-pkg-verbose
	(message "Configuring F90 mode."))
    (if config-pkg-colorize-modeline 
	(add-hook 'f90-mode-hook        
		  '(lambda () 
		     (set-face-background 'modeline 
					  "orange" (current-buffer))
		     (set-face-foreground 'modeline 
					  "black"   (current-buffer)))))
    (setq auto-mode-alist
	  (append
	   '(("\\.f90$"  . f90-mode)
	     ("\\.F$"    . f90-mode)
	     ("\\.FH$"   . f90-mode)
	     ("\\.fm4$"  . f90-mode)
	     ) auto-mode-alist))
    ;; let .F denone Fortran and not freeze files
    (defvar crypt-freeze-vs-fortran nil)
    (add-hook 'f90-mode-hook 'turn-on-auto-fill)))

(if want-f90-mode (draco-init-f90-mode))

;; ========================================
;; FORTRAN
;; ========================================

(defvar want-fortran-mode nil
  "*Does the user want to have fortran-mode?")

(defun draco-init-fortran-mode ()
  "Autoload fortran-mode and append the approriate suffixes to
auto-mode-alist."
  (interactive)
  (progn
    (if config-pkg-verbose
	(message "Configuring FORTRAN mode."))
    (setq auto-mode-alist
	  (append
	   '(("\\.for$"  . fortran-mode)
	     ("\\.f$"    . fortran-mode) 
	     ("\\.id$"   . fortran-mode)
	     ("\\.fh$"   . fortran-mode)
	     ) auto-mode-alist))
    ;; let .F denone Fortran and not freeze files
    (defvar crypt-freeze-vs-fortran nil)
    (add-hook 'fortran-mode-hook 'turn-on-auto-fill)))
  
(if want-fortran-mode (draco-init-fortran-mode))


;; ========================================
;; ChangeLog
;; ========================================

(defvar want-change-log-mode nil
  "*Does the user want to have change-log-mode?")

(defun draco-init-change-log-mode ()
  "Autoload change-log-mode and append the approriate suffixes to
auto-mode-alist."
  (interactive)
  (progn
    (if config-pkg-verbose
	(message "Configuring ChangeLog mode."))
    (add-hook 'change-log-mode-hook 'turn-on-font-lock)
    (if config-pkg-colorize-modeline 
	(add-hook 'change-log-mode-hook        
		  '(lambda () 
		       (set-face-background 'modeline 
					    "bisque3" (current-buffer))
		       (set-face-foreground 'modeline 
					    "black"   (current-buffer)))))
    (setq auto-mode-alist
	  (append
	   '(("ChangeLog"  . change-log-mode)
	     ) auto-mode-alist))
 (add-hook 'change-log-mode-hook 'turn-on-auto-fill)))
  
(if want-change-log-mode (draco-init-change-log-mode))

;; ========================================
;; Emacs Lisp
;; ========================================

(defvar want-emacs-lisp-mode nil
  "*Does the user want to have emacs-lisp-mode?")

(defun draco-init-emacs-lisp-mode ()
  "Autoload emacs-lisp-mode, append the approriate suffixes to
auto-mode-alist and set up some customizations for RTT."
  (interactive)
  (progn
    (if config-pkg-verbose
	(message "Configuring emacs-lisp-mode."))
    (autoload 'emacs-lisp-mode "emacs-lisp-mode" "Emacs Lisp Editing Mode" t)
    (add-hook 'emacs-lisp-mode-hook 'turn-on-font-lock)
    (if config-pkg-colorize-modeline 
	(add-hook 'emacs-lisp-mode-hook        
		  '(lambda () 
		     (set-face-background 'modeline 
					  "tan" (current-buffer))
		     (set-face-foreground 'modeline 
					  "black"   (current-buffer)))))
    (setq auto-mode-alist
	  (append
	   '(("\\.el$"  . emacs-lisp-mode)
	     ) auto-mode-alist))
    (defun rtt-elisp-mode-hook ()
      "Hooks added to Elisp mode"
      (local-set-key [(f5)] 'tme-elisp-divider)
      (local-set-key [(f6)] 'tme-elisp-comment-divider))
    (add-hook 'emacs-lisp-mode-hook 'rtt-elisp-mode-hook)
    (add-hook 'emacs-lisp-mode-hook 'turn-on-font-lock)
    (add-hook 'emacs-lisp-mode-hook 'turn-on-auto-fill)))
  
(if want-emacs-lisp-mode (draco-init-emacs-lisp-mode))

;; ========================================
;; Interactive Shell
;; ========================================

(defvar want-shell-mode nil
  "*Does the user want to have shell-mode?")

(defun draco-init-shell-mode ()
  "Autoload shell-mode and set up some customizations for RTT."
  (interactive)
  (progn
    (if config-pkg-verbose
	(message "Configuring shell-mode."))
    (autoload 'shell-mode "shell-mode" "Interactive Shell Mode" t)
    (if config-pkg-colorize-modeline 
	(add-hook 'shell-mode-hook        
		  '(lambda () ;; M-x list-colors-display
		     (set-face-background 'modeline 
					  "thistle" (current-buffer))
		     (set-face-foreground 'modeline 
					  "black"   (current-buffer)))))
    (add-hook 'shell-mode-hook 'turn-on-font-lock)))
  
(if want-shell-mode (draco-init-shell-mode))

;; ========================================
;; CVS Mode
;; ========================================

(defvar want-cvs-mode nil
  "*Does the user want to have cvs-mode?")

(defun draco-init-cvs-mode ()
  "Autoload cvs-mode and set up some customizations for RTT."
  (interactive)
  (progn
    (if config-pkg-verbose
	(message "Configuring cvs-mode."))
    (autoload 'cvs-mode "cvs-mode" "Interactive Shell Mode" t)
    (defun rtt-cvs-edit-mode-hook ()
      "Setup the PCL-CVS cvs-edit-mode with rtt prefs."
      (auto-fill-mode t)
      (setq fill-prefix "  "))
    (add-hook 'cvs-edit-mode-hook 'rtt-cvs-edit-mode-hook)
    (if config-pkg-colorize-modeline 
	(add-hook 'cvs-mode-hook        
		  '(lambda () ;; M-x list-colors-display
		     (set-face-background 'modeline 
					  "honeydew" (current-buffer))
		     (set-face-foreground 'modeline 
					  "black"   (current-buffer)))))
    (defvar cvs-erase-input-buffer        nil)
    (defvar cvs-inhibit-copyright-message t  )
    ; If this variable is set to any non-nil value
    ; `cvs-mode-remove-handled' will be called every time you check in
    ; files, after the check-in is ready. See section 5.11 Removing handled
    ; entries.
    (defvar cvs-auto-remove-handled t)

    ; If this variable is set to any non-nil value, directories that do not
    ; contain any files to be checked in will not be listed in the `*cvs*'
    ; buffer. 
    (defvar cvs-auto-remove-handled-directories t)
    (require 'pcl-cvs)
    (define-key global-map [(meta shift f5)]
      'cvs-mode-add-change-log-entry-other-window)
    )
  )

(if want-cvs-mode (draco-init-cvs-mode))

;; ========================================
;; Doxymacs Mode
;; ========================================

(defvar want-doxymacs-mode t
  "*Does the user want to have doxymacs-mode?")

(defun draco-init-doxymacs-mode ()
  "Autoload doxymacs-mode and set up some customizations for RTT."
  (interactive)
  (progn
    (if config-pkg-verbose
	(message "Configuring doxymacs-mode."))
    (autoload 'doxymacs-mode "doxymacs-mode" "Interactive Shell Mode" t)
    (require 'doxymacs)
    (defvar doxymacs-doxygen-style "Qt")
    (add-hook 'c-mode-common-hook 'doxymacs-mode)
    (defun rtt-doxymacs-font-lock-hook ()
      (if (or (eq major-mode 'c-mode) (eq major-mode 'c++-mode))
	  (doxymacs-font-lock)))
     (add-hook 'font-lock-mode-hook 'rtt-doxymacs-font-lock-hook)
    )
  )

(if want-doxymacs-mode (draco-init-doxymacs-mode))

;; ========================================
;; Shell mode
;; ========================================

(defvar want-sh-mode nil
  "*Does the user want to have sh-mode?")

(defun draco-init-sh-mode ()
  "Autoload sh-mode, append the approriate suffixes to
auto-mode-alist and set up some customizations for RTT."
  (interactive)
  (progn
    (if config-pkg-verbose
	(message "Configuring sh-mode."))
    (autoload 'sh-mode "sh-mode" "Bourne Shell Editing Mode" t)
    (add-hook 'sh-mode-hook 'turn-on-font-lock)
    (if config-pkg-colorize-modeline 
	(add-hook 'sh-mode-hook        
		  '(lambda () 
		     (set-face-background 'modeline 
					  "palegoldenrode" (current-buffer))
		     (set-face-foreground 'modeline 
					  "black"   (current-buffer)))))
    (setq auto-mode-alist
	  (append
	   '(("\\.bash." . sh-mode)
	     ) auto-mode-alist))
    (require 'sh-script)
    (sh-set-shell "bash")
    (defun rtt-sh-mode-hook ()
      "Hooks added to shell mode"
      (local-set-key [(f5)] 'tme-makefile-divider)
      (local-set-key [(f6)] 'tme-makefile-comment-divider))
    (add-hook 'sh-mode-hook 'rtt-sh-mode-hook)
    (add-hook 'sh-mode-hook 'turn-on-font-lock)
    (add-hook 'sh-mode-hook 'turn-on-auto-fill)))

(if want-sh-mode (draco-init-sh-mode))

;; ========================================
;; SGML mode
;; ========================================

(defvar want-sgml-mode nil
  "*Does the user want to have sgnk-mode?")

(defun draco-init-sgml-mode ()
  "Autoload sgml-mode, append the approriate suffixes to
auto-mode-alist and set up some customizations for RTT."
  (interactive)
  (progn
    (if config-pkg-verbose
	(message "Configuring sgml-mode."))
    (autoload 'sgml-mode "sgml-mode" "SGML Shell Editing Mode" t)
    (add-hook 'sgml-mode-hook 'turn-on-font-lock)
    (if config-pkg-colorize-modeline 
	(add-hook 'sgml-mode-hook        
		  '(lambda () 
		     (set-face-background 'modeline 
					  "thistle" (current-buffer))
		     (set-face-foreground 'modeline
					  "black"   (current-buffer)))))
    (setq auto-mode-alist
	  (append
	   '(("\\.sgml$" . sh-mode)
	     ) auto-mode-alist))
;    (defun rtt-sgml-mode-hook ()
;      "Hooks added to shell mode"
;      (local-set-key [(f5)] 'tme-makefile-divider)
;      (local-set-key [(f6)] 'tme-makefile-comment-divider))
;    (add-hook 'sgml-mode-hook 'rtt-sgml-mode-hook)
    (add-hook 'sgml-mode-hook 'turn-on-font-lock)
    (add-hook 'sgml-mode-hook 'turn-on-auto-fill)))

(if want-sgml-mode (draco-init-sgml-mode))

;;---------------------------------------------------------------------------;;
;; set up Config-pkg
;;---------------------------------------------------------------------------;;

(defun Config-pkg ()
  (message "Configuring DRACO menu."))

(provide 'Config-pkg)
