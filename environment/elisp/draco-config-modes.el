;; ======================================================================
;; draco-config-modes.el
;; Kelly Thompson
;; 8 Dec 2004
;;
;; $Id$
;;
;; Configure a variety of packages, upon request of user.
;;
;; Must load the draco-setup.el package from XEmacs startup files
;; (init.el or .emacs).   Then use the XEmacs menubar to select 
;; Options --> Advanced --> Group.  At the prompt enter "draco".  Set 
;; and save options for Draco elisp setup.
;; ======================================================================

(require (quote draco-mode))

;; ========================================
;; Compilation-mode
;; ========================================

(defun draco-setup-compilation-mode ()
  "Autoload compilation-mode and add mode local keybindings to the compilation-mode-hook."
  (interactive)
  (progn
;    (autoload 'compilation-mode "compilation-mode" nil t)
    (defun draco-compilation-mode-hook ()
      "draco-mode hooks added to MPPL mode."
      (turn-on-draco-mode)
      (turn-on-auto-fill))
    (add-hook 'compilation-mode-hook 'draco-compilation-mode-hook)))

;; ========================================
;; MPPL
;; ========================================

(defun draco-setup-mppl-mode ()
  "Autoload mppl-mode, append appropriate suffixes to auto-mode-alist
and add turn-on-auto-fill to the mppl-mode-hook."
  (interactive)
  (progn
    (autoload 'mppl-mode   "mppl-mode" nil t)
    (setq auto-mode-alist
	  (append '(("\\.i$" . mppl-mode)
		    ("\\.p$" . mppl-mode)
		    ("\\.m$" . mppl-mode)
		    ("\\.pm4$" . mppl-mode)
		    ) auto-mode-alist))
    (defun draco-mppl-mode-hook ()
      "draco-mode hooks added to MPPL mode."
      (turn-on-draco-mode)
      (turn-on-auto-fill))
    (add-hook 'mppl-mode-hook 'draco-mppl-mode-hook)))

;; ========================================
;; TCL
;; ========================================

(defun draco-setup-tcl-mode ()
  "Autoload tcl-mode and append the appropriate suffixes to 
auto-mode-alist."
  (interactive)
  (progn
    (autoload 'tcl-mode   "tcl-mode" nil t)
    (setq auto-mode-alist
	  (append '(("\\.tcl$" . tcl-mode)
		    ("\\.itk$" . tcl-mode)
		    ("\\.ith$" . tcl-mode)
		    ("\\.itm$" . tcl-mode)
		    ) auto-mode-alist))
    (defun draco-tcl-mode-hook ()
      "draco-mode hooks added to TCL mode."
      (turn-on-draco-mode)
      (turn-on-auto-fill))
    (add-hook 'tcl-mode-hook 'draco-tcl-mode-hook)))

;; ========================================
;; Python
;; ========================================

(defun draco-setup-python-mode ()
  "Autoload python-mode and append the appropriate suffixes to
auto-mode-alist."
  (interactive)
      (progn
      (autoload 'python-mode "python-mode" "Python editing mode." t)
      (setq auto-mode-alist
	    (cons '("\\.py$" . python-mode) auto-mode-alist))
      (defun draco-python-mode-hook ()
	"DRACO hooks added to Python mode."
	(defvar py-indent-offset 4)
	(defvar indent-tabs-mode nil)
	(local-set-key [(control c)(control c)] 'comment-region)
	(local-set-key [(f5)] 'draco-makefile-divider)
	(local-set-key [(f6)] 'draco-makefile-comment-divider)
	(turn-on-draco-mode)
	(turn-on-auto-fill))
      (add-hook 'python-mode-hook 'draco-python-mode-hook)))

;; ========================================
;; CMake
;; ========================================

(defun draco-setup-cmake-mode ()
  "Autoload cmake-mode and append the appropriate suffixes to
auto-mode-alist."
  (interactive)
  (progn
    (autoload 'cmake-mode "cmake-mode" "CMake editing mode." t)
    (setq auto-mode-alist
          (append '(("\\.cmake"         . cmake-mode)
                    ("CMakeLists\\.txt" . cmake-mode)
                    ("CMakeCache\\.txt" . cmake-mode)
                    ("\\.cmake\\.in"    . cmake-mode))
                  auto-mode-alist))
    (defun draco-cmake-mode-hook ()
      "DRACO hooks added to CMake mode."
      (defvar cmake-tab-width 2)
      (local-set-key [(control c)(control c)] 'comment-region)
      (local-set-key [(f5)] 'draco-makefile-divider)
      (local-set-key [(f6)] 'draco-makefile-comment-divider)
      (turn-on-draco-mode)
      (turn-on-auto-fill))
    (add-hook 'cmake-mode-hook 'draco-cmake-mode-hook)))

;; ========================================
;; Autoconf
;; ========================================

(defun draco-setup-autoconf-mode ()
  "Autoload autoconf-mode and append the appropriate suffixes to
auto-mode-alist."
  (interactive)
      (progn
	(setq auto-mode-alist
	      (append '(("\\.m4$" . autoconf-mode)
			("\\.ac$" . autoconf-mode)
			("\\.in$" . autoconf-mode))
		      auto-mode-alist))
	
	(defun draco-menu-insert-comments-m4 () 
	  "Submenu for inserting comments (context sensitive)."
	  (list "Insert comment..."
		["Insert m4 divider" draco-m4-divider t]
		["Insert m4 comment divider" draco-m4-comment-divider t]
		["Insert Makefile divider" draco-makefile-divider t]
		["Insert Makefile comment divider" draco-makefile-comment-divider t]))
	
	(defun draco-autoconf-mode-hook ()
	  "DRACO hooks added to autoconf mode"
	  (setq tab-stop-list '(3 7 11 15 19 23 27 31 35 39 43 47 51 55 59 63 67 71 75 79 83))
	  (local-set-key [(f5)] 'draco-m4-divider)
	  (local-set-key [(f6)] 'draco-m4-comment-divider)
	  (local-set-key [(shift f5)] 'draco-makefile-divider)
	  (local-set-key [(shift f6)] 'draco-makefile-comment-divider)
	  (local-set-key [(tab)] 'tab-to-tab-stop)
	  (draco-mode-update-menu (draco-menu-insert-comments-m4))
	  (turn-on-auto-fill)
	  (turn-on-draco-mode))
	(add-hook 'autoconf-mode-hook 'draco-autoconf-mode-hook)))

;; ========================================
;; Makefile
;; ========================================

(defun draco-setup-makefile-mode ()
  "Autoload makefile-mode and append the appropriate suffixes to
auto-mode-alist.

- Autoload make-mode.
- Register files named Makefile.* and makefile.* with this mode.
- Create menu items for inserting Makefile dividers.
- Set keybindings:
  [ f5 ] - Insert Makefile divider
  [ f6 ] - Insert Makefile comment divider
- Turn on draco-mode
- Turn on font-lock
- Turn on auto-fill"
  (interactive)
  (progn
    (autoload 'makefile-mode   "make-mode" nil t)
    (setq auto-mode-alist
	  (append '(("makefile.*" . makefile-mode)
		    ("Makefile.*" . makefile-mode)
		    ) auto-mode-alist))
    (defun draco-menu-insert-comments-makefile () 
      "Submenu for inserting comments (context sensitive)."
      (list "Insert comment..."
	    ["Insert Makefile divider" draco-makefile-divider t]
	    ["Insert Makefile comment divider"
	     draco-makefile-comment-divider t]))
    (defun draco-makefile-mode-hook ()
      "DRACO hooks added to Makefile mode."
      (draco-mode-update-menu (draco-menu-insert-comments-makefile))
      (local-set-key [(f5)] 'draco-makefile-divider)
      (local-set-key [(f6)] 'draco-makefile-comment-divider)
      (turn-on-font-lock)
      (turn-on-auto-fill)
      (turn-on-draco-mode))
    (add-hook 'makefile-mode-hook 'draco-makefile-mode-hook t)))

;; ========================================
;; C++
;; ========================================

(defun draco-setup-cc-mode ()
  "Autoload c++-mode, c-mode and append the appropriate suffixes to 
auto-mode-alist.

- Autoload c-mode, c++-mode and doxymacs-mode.
- Associate files *.C, *.cc, *.pt, *.hh, *.hpp, *.cpp, *.hh.in,
  *.h.in, *.c, *.h, *.dcc, *.dcc.in and *.dot with this mode.
- Create and install menu items for inserting C++/C/Doxygen comment
  blocks.
- Set the C++ indentation style to match Draco source code.
- Setup colorized modeline (if requested).
- Set fill-column to 78
- Set hotkeys:
  [ f5 ]   - Insert C++ divider
  [ f6 ]   - Insert C++ comment divider
  [ M-f5 ] - Insert C++ function comment
  [ M-f6 ] - Insert C++ class comment
  [ S-f5 ] - Insert doxygen function comment
  [ S-f6 ] - Insert doxygen file comment
  [ C-f5 ] - Insert doxygen multiline comment
  [ C-f6 ] - Insert doxygen singleline comment
- Turn on font lock.
- Turn on auto fill.
- Turn on draco-mode."
  (interactive)
    (progn
      (autoload 'c++-mode "cc-mode" "C++ Editing Mode" t)
      (autoload 'c-mode   "cc-mode" "C Editing Mode" t)
      (autoload 'doxymacs-mode "doxymacs-mode" "Doxygen Editing Mode" t)
      (setq auto-mode-alist
	    (append '(("\\.C$"      . c++-mode)
		      ("\\.cc$"     . c++-mode)
                      ("\\.pt$"     . c++-mode)
		      ("\\.hh$"     . c++-mode)
		      ("\\.hpp$"    . c++-mode)
		      ("\\.cpp$"    . c++-mode)
                      ("\\.hh.in$"  . c++-mode)
		      ("\\.h.in$"   . c-mode)
		      ("\\.c$"      . c-mode)   ; to edit C code
		      ("\\.h$"      . c-mode)   ; to edit C code
		      ("\\.dcc$"    . c-mode)   ; to edit C code
		      ("\\.dcc.in$" . c-mode)   ; to edit C code
 		      ("\\.dot$"    . c-mode)  ; for dot files
		      ) auto-mode-alist))
      (defun draco-menu-insert-comments-cc () 
	"Submenu for inserting comments (context sensitive)."
	(list "Insert comment..."
	      ["Insert C++ divider" draco-insert-comment-divider t]
	      ["Insert C++ comment block"   draco-cc-divider t]
	      ["Insert C++ function divider" draco-insert-function-doc t]
	      ["Insert C++ class comment block"   draco-insert-class-doc       t]
	      ["Insert Doxygen singleline comment" doxymacs-insert-blank-singleline-comment t ]
	      ["Insert Doxygen multiline comment" doxymacs-insert-blank-multiline-comment t ]
	      ["Insert Doxygen file comment" doxymacs-insert-file-comment t ]
	      ["Insert Doxygen function comment" doxymacs-insert-function-comment t ]
	      ["Insert Doxygen member comment" doxymacs-insert-member-comment t ]
	      ["Insert Doxygen grouping comment" doxymacs-insert-grouping-comment t ]
	      ["Insert C divider" draco-c-comment-divider t]
	      ["Insert C comment block"   draco-c-divider t]
	      ))

            ;; Borrowed from
            ;; http://www.esperi.demon.co.uk/nix/xemacs/personal/init-prog-modes.html

            ;; Also see help for XEmacs variable c-offsets-alist 
            ;; \C-h v c-offset-alist
      (defun draco-setup-c-mode ()
	"Setup C, C++ mode for Draco Developers.

This is run in the C-mode-common-hook to set up indentation and other
parameters on creation of buffers managed by cc-mode.el for Nix's personal coding style."
	(c-add-style 
	 "draco" '
	 (
          ; Tab indent == 4 spaces
	  (c-basic-offset . 4)       
          ; K&R? Blugh. Not usin' *that*.
	  ;(c-recognize-knr-p . nil)
          ; Do nil for lone comments
	  ;(c-comment-only-line-offset . (0 . 0)) 
          ; We don't use *-prefixed comments
  	  ;(c-block-comment-prefix . "") 
          ; Even with no code before them
	  ;(c-indent-comments-syntactically-p . t) 
          ; Make function calls look nice
	  ;(c-cleanup-list . (space-before-funcall compact-empty-funcall)) 
          ; Snap #s to the first column
	  (c-electric-pound-behavior . 'alignleft) 
          ; Regexp to find the starting brace of a block
	  ;(defun-prompt-regexp . " ") 
	  (c-offsets-alist . (
			      (access-label . -2 )
			      (block-close . 0)
			      (block-open  . 0)
			      (case-label  . +)
			      (class-close . 0)
			      (class-open  . 0)
			      (defun-close . 0)
			      (defun-open  . 0)
			      (do-while-closure  . 0)
			      (else-clause       . 0)
			      ;;(extern-lang-close . +)
			      ;;(extern-lang-open  . +)
			      (extern-lang-close . 0)
			      (extern-lang-open  . 0)
			      (inline-close      . 0)
			      (inline-open       . 0)
			      (innamespace       . 0)
			      (statement-case-intro . +)
			      (statement-cont    . c-lineup-math)
			      (substatement-open . 0)
			      )))))
;      (defun draco-clean-up-common-hook ()
;	"Clean up the common hook.

;Removes the other style setup functions from the hook, so that they
;only get run once, rather than repeatedly."
;	(remove-hook 'c-mode-common-hook 'draco-setup-c-mode)
;	(remove-hook 'c-mode-common-hook 'draco-clean-up-common-hook))

;      (defun draco-setup-this-c-mode-buffer ()
;	"Set up this buffer for C mode.
;Things (like auto-hungry-state setting and style setting) that should not
;be removed from the `c-mode-common-hook' after the first call by
;`draco-clean-up-c-common-hook', but which should rather take effect
;separately for each buffer."
;	(if (gawd-personal-code-p)
;	    (progn
;	      (c-make-styles-buffer-local t) ; Don't let this interfere with user styles
;	      (c-toggle-hungry-state 1)
;	      (c-set-style "draco"))))



      (if draco-colorize-modeline 
	  (add-hook 'c++-mode-hook        
		    '(lambda () 
		       (set-face-background 'modeline 
					    "skyblue" (current-buffer))
		       (set-face-foreground 'modeline 
					    "black"   (current-buffer)))))
      (if draco-colorize-modeline 
	  (add-hook 'c-mode-hook        
		    '(lambda () 
		       (set-face-background 'modeline 
					    "pink" (current-buffer))
		       (set-face-foreground 'modeline 
					    "black"   (current-buffer)))))
      
      (defun draco-c-mode-hook ()
	"DRACO hooks added to C/C++ mode.

- Sets c-style to \"draco\"
- Sets fill-column to 78
- Sets f5/f6 as hot keys to insert dividers.
- Turns on auto-fill"
	(draco-setup-c-mode)
	(c-set-style "draco")
	(local-set-key "\C-m" 'newline-and-indent)
	(set-fill-column 78)
	(local-set-key [(f5)] 'draco-cc-divider)
	(local-set-key [(f6)] 'draco-insert-comment-divider)
	(local-set-key [(meta f5)] 'draco-insert-function-doc)
	(local-set-key [(meta f6)] 'draco-insert-class-doc)
	(local-set-key [(shift f5)] 'doxymacs-insert-function-comment)
	(local-set-key [(shift f6)] 'doxymacs-insert-file-comment)
	(local-set-key [(control f5)] 'doxymacs-insert-blank-multiline-comment)
	(local-set-key [(control f6)] 'doxymacs-insert-blank-singleline-comment)
	(draco-mode-update-menu (draco-menu-insert-comments-cc))
	(turn-on-font-lock)
	(turn-on-auto-fill))
      (add-hook 'c-mode-common-hook 'draco-c-mode-hook)
      (add-hook 'c-mode-common-hook 'imenu-add-menubar-index)
;      (add-hook 'c-mode-common-hook 'draco-add-style)
      (add-hook 'c-mode-common-hook 'turn-on-draco-mode)
      (add-hook 'font-lock-mode-hook
		'(lambda ()
		   (if (or (eq major-mode 'c-mode) (eq major-mode 'c++-mode))
		       (draco-font-lock))))
      ))


;; ========================================
;; AUCTEX
;; ========================================

(defun draco-setup-auctex-mode ()
  "Loads the tex-site package."
  (interactive)
    (progn
      ; (require 'tex-site)
      ; (require 'reftex)
      ; (require 'bib-cite)

    (defun draco-menu-insert-comments-tex () 
      "Submenu for inserting comments (context sensitive)."
      (list "Insert comment..."
	    ["Insert LaTeX divider" draco-latex-divider t]
	    ["Insert LaTeX comment divider" draco-latex-comment-divider t]))

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

      (setq reftex-enable-partial-scans           t
	    reftex-save-parse-info                t
	    reftex-use-multiple-selection-buffers t
	    reftex-plug-into-AUCTeX               t)

      (setq reftex-texpath-environment-variables draco-texpath)
      (setq reftex-bibpath-environment-variables draco-bibpath)
      
      (defun draco-latex-mode-hook ()
	"DRACO hooks added to LaTeX and BibTex modes."
	(local-set-key [(f5)] 'draco-latex-divider)
	(local-set-key [(f6)] 'draco-latex-comment-divider)
	(local-set-key "\C-c %" 'comment-region)
	(draco-mode-update-menu (draco-menu-insert-comments-tex))
	(turn-on-bib-cite)
	(turn-on-reftex)
	(turn-on-auto-fill)
	(turn-on-draco-mode)
	)
      (add-hook 'TeX-mode-hook  'draco-latex-mode-hook)
      (add-hook 'bibtex-mode-hook 'draco-latex-mode-hook)
      (add-hook 'tex-mode-hook  'turn-on-draco-mode)
      (add-hook 'tex-mode-hook  'draco-latex-mode-hook)
;      (add-hook 'LaTeX-mode-hook  'turn-on-draco-mode)
;      (add-hook 'LaTeX-mode-hook  'turn-on-auto-fill)))
      ))

;; ========================================
;; FORTRAN-90
;; ========================================

(defun draco-setup-f90-mode ()
  "Autoload f90-mode and append the approriate suffixes to
auto-mode-alist."
  (interactive)
  (progn

    (defun draco-menu-insert-comments-f90 () 
      "Submenu for inserting comments (context sensitive)."
      (list "Insert comment..."
	    ["Insert f90 subroutine divider" draco-f90-subroutine-divider       t]
	    ["Insert f90 comment divider" draco-f90-comment-divider t]))
    
    (if draco-colorize-modeline 
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

    (defun draco-f90-mode-hook ()
      "Hooks added to F90 mode"
      (local-set-key [(f5)]         'draco-f90-subroutine-divider)
      (local-set-key [(control f6)] 'draco-f90-insert-document)
      (local-set-key [(f6)]         'draco-f90-comment-divider)
      (draco-mode-update-menu (draco-menu-insert-comments-f90))
      (set-fill-column 80))
     ;; let .F denone Fortran and not freeze files
    (defvar crypt-freeze-vs-fortran nil)
    (add-hook 'f90-mode-hook 'draco-f90-mode-hook)
    (add-hook 'f90-mode-hook 'turn-on-draco-mode)
    (add-hook 'f90-mode-hook 'turn-on-auto-fill)))

;; ========================================
;; FORTRAN
;; ========================================

(defun draco-setup-fortran-mode ()
  "Autoload fortran-mode and append the approriate suffixes to
auto-mode-alist."
  (interactive)
  (progn

    (defun draco-menu-insert-comments-f77 () 
      "Submenu for inserting comments (context sensitive)."
      (list "Insert comment..."
	    ["Insert f77 subroutine divider" draco-f77-subroutine-divider t]
	    ["Insert f77 comment divider" draco-f77-comment-divider t]))

    (if draco-colorize-modeline 
	(add-hook 'fortran-mode-hook        
		  '(lambda () 
		     (set-face-background 'modeline 
					  "orange" (current-buffer))
		     (set-face-foreground 'modeline 
					  "black"
					  (current-buffer)))))
    ;; let .F denone Fortran and not freeze files
    (defvar crypt-freeze-vs-fortran nil)
    (setq auto-mode-alist
	  (append
	   '(("\\.for$"  . fortran-mode)
	     ("\\.f$"    . fortran-mode) 
	     ("\\.g$"    . fortran-mode) 
	     ("\\.id$"   . fortran-mode)
	     ("\\.fh$"   . fortran-mode)
	     ) auto-mode-alist))
    (defun draco-fortran-mode-hook ()
      "Hooks added to F77 mode"
      (local-set-key [(f5)]         'draco-f77-subroutine-divider)
      (local-set-key [(control f6)] 'draco-f77-insert-document)
      (local-set-key [(f6)]         'draco-f77-comment-divider)
      (draco-mode-update-menu (draco-menu-insert-comments-f77)))
    (add-hook 'fortran-mode-hook 'draco-fortran-mode-hook)
    (add-hook 'fortran-mode-hook 'turn-on-draco-mode)
    (add-hook 'fortran-mode-hook 'turn-on-auto-fill)))
  
;; ========================================
;; ChangeLog
;; ========================================

(defun change-log-font-lock ()
  "Turn on font-lock for ChangeLog keywords."
  (interactive)
  (let 
      ((old 
	(if (eq (car-safe font-lock-keywords) t)
		(cdr font-lock-keywords)
		    font-lock-keywords)))
    (setq font-lock-keywords (append old change-log-font-lock-keywords))
    ))

(defun draco-setup-change-log-mode ()
  "Autoload change-log-mode and append the approriate suffixes to
auto-mode-alist."
  (interactive)
  (progn
    (autoload 'change-log-mode "change-log-mode" 
      "ChangeLog Editing Mode" t)
    (require 'add-log)
    (if draco-colorize-modeline 
	(add-hook 'change-log-mode-hook        
		  '(lambda () 
		     (set-face-background 'modeline 
					  "bisque3" (current-buffer))
		     (set-face-foreground 'modeline 
					  "black"
					  (current-buffer)))))
    (setq auto-mode-alist
	  (append
	   '(("ChangeLog"  . change-log-mode)
	     ) auto-mode-alist))
    (add-hook 'change-log-mode-hook 'turn-on-font-lock)
    (add-hook 'font-lock-mode-hook
	      '(lambda ()
		 (if (eq major-mode 'change-log-mode)
		     (change-log-font-lock))))
;    (add-hook 'change-log-mode-hook 'turn-on-draco-mode)
; See http://www.emacswiki.org/cgi-bin/wiki?OutlineMode
    (add-hook 'change-log-mode-hook
              '(lambda ()
                 (set (make-local-variable 'outline-regexp) 
                      "[[:digit:]]+")))
    (add-hook 'change-log-mode-hook 'turn-on-auto-fill)))
  
;; ========================================
;; Emacs Lisp
;; ========================================

(defun draco-setup-emacs-lisp-mode ()
  "Autoload emacs-lisp-mode, append the approriate suffixes to
auto-mode-alist and set up some customizations for DRACO."
  (interactive)
  (progn
    (autoload 'emacs-lisp-mode "emacs-lisp-mode" "Emacs Lisp Editing Mode" t)

    (defun draco-menu-insert-comments-elisp () 
      "Submenu for inserting comments (context sensitive)."
      (list "Insert comment..."
	    ["Insert elisp comment block"   draco-elisp-divider       t]
	    ["Insert elisp comment divider" draco-elisp-comment-divider t]))
    
    (setq auto-mode-alist
	  (append
	   '(("\\.el$"  . emacs-lisp-mode)
	     (".emacs$"  . emacs-lisp-mode)
	     ) auto-mode-alist))
    
    (if draco-colorize-modeline 
	(add-hook 'emacs-lisp-mode-hook        
		  '(lambda () 
		     (set-face-background 'modeline 
					  "tan" (current-buffer))
		     (set-face-foreground 'modeline 
					  "black"   (current-buffer)))))
    (defun draco-elisp-mode-hook ()
      "Hooks added to Elisp mode"
      (local-set-key [(f5)] 'draco-elisp-divider)
      (local-set-key [(f6)] 'draco-elisp-comment-divider)
      (draco-mode-update-menu (draco-menu-insert-comments-elisp))
      (local-set-key [(control c)(control c)] 'comment-region))
    (add-hook 'emacs-lisp-mode-hook 'turn-on-draco-mode)
    (add-hook 'emacs-lisp-mode-hook 'draco-elisp-mode-hook)
    (add-hook 'emacs-lisp-mode-hook 'turn-on-font-lock)
    (add-hook 'emacs-lisp-mode-hook 'turn-on-auto-fill)))

;; ========================================
;; Interactive Shell
;; ========================================

(defun draco-setup-shell-mode ()
  "Autoload shell-mode and set up some customizations for DRACO."
  (interactive)
  (progn
    (autoload 'shell-mode "shell-mode" "Interactive Shell Mode" t)
    (if draco-colorize-modeline 
	(add-hook 'shell-mode-hook        
		  '(lambda () ;; M-x list-colors-display
		     (set-face-background 'modeline 
					  "thistle" (current-buffer))
		     (set-face-foreground 'modeline 
					  "black"   (current-buffer)))))
    (add-hook 'shell-mode-hook 'turn-on-draco-mode)
    (add-hook 'shell-mode-hook 'turn-on-font-lock)))
  
;; ========================================
;; CVS Mode
;; http://www.xemacs.org/Documentation/packages/html/pcl-cvs_5.html#SEC13
;; ========================================

(defun draco-setup-cvs-mode ()
  "Autoload cvs-mode and set up some customizations for DRACO."
  (interactive)
  (progn
;    (autoload 'cvs-examine "pcl-cvs" "CVS mode" t)
;    (require 'pcl-cvs)

    (defun draco-menu-extras-cvs () 
      "Submenu for inserting comments (context sensitive)."
      (list "CVS extras..."
	    ["CVS help" cvs-help t]))

    (defun draco-cvs-edit-mode-hook ()
      "Setup the PCL-CVS cvs-edit-mode with draco prefs."
      (auto-fill-mode t)
      (setq fill-prefix "  ")
      (draco-mode-update-menu (draco-menu-extras-cvs)))
    (add-hook 'cvs-mode-hook 'draco-cvs-edit-mode-hook)
    (add-hook 'cvs-mode-hook 'turn-on-draco-mode)
    (if draco-colorize-modeline 
	(add-hook 'cvs-mode-hook        
		  '(lambda () ;; M-x list-colors-display
		     (set-face-background 'modeline 
					  "honeydew" (current-buffer))
		     (set-face-foreground 'modeline 
					  "black"   (current-buffer)))))
    (setq cvs-erase-input-buffer        nil
	  cvs-inhibit-copyright-message t  
	  cvs-status-flags "-q"
	  vc-dired-terse-display nil )
    ; If this variable is set to any non-nil value
    ; `cvs-mode-remove-handled' will be called every time you check in
    ; files, after the check-in is ready. See section 5.11 Removing handled
    ; entries.
    (setq cvs-auto-remove-handled t)

    ; If this variable is set to any non-nil value, directories that do not
    ; contain any files to be checked in will not be listed in the `*cvs*'
    ; buffer. 
    (setq cvs-auto-remove-handled-directories t)
    (setq cvs-commit-buffer-require-final-newline t)
    ))

;; ========================================
;; Doxymacs Mode
;; ========================================

(defun draco-setup-doxymacs-mode ()
  "Autoload doxymacs-mode and set up some customizations for DRACO."
  (interactive)
  (progn
    (autoload 'doxymacs-mode "doxymacs-mode" "Doxygen Editing Mode" t)
    (require 'doxymacs)
    (defvar doxymacs-doxygen-style "Qt")
    (add-hook 'c-mode-common-hook 'doxymacs-mode)
    (defun draco-doxymacs-font-lock-hook ()
      (if (or (eq major-mode 'c-mode) (eq major-mode 'c++-mode))
	  (doxymacs-font-lock)))
     (add-hook 'font-lock-mode-hook 'draco-doxymacs-font-lock-hook)
    )
  )

;; ========================================
;; Shell mode
;; ========================================

(defun draco-setup-sh-mode ()
  "Autoload sh-mode, append the approriate suffixes to
auto-mode-alist and set up some customizations for DRACO."
  (interactive)
  (progn
    (autoload 'sh-mode "sh-mode" "Bourne Shell Editing Mode" t)

    (defun draco-menu-insert-comments-makefile () 
      "Submenu for inserting comments (context sensitive)."
      (list "Insert comment..."
	    ["Insert Makefile divider"         draco-makefile-divider         t]
	    ["Insert Makefile comment divider" draco-makefile-comment-divider t]))

    (if draco-colorize-modeline 
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
    (defun draco-sh-mode-hook ()
      "Hooks added to shell mode"
      (local-set-key [(f5)] 'draco-makefile-divider)
      (local-set-key [(f6)] 'draco-makefile-comment-divider)
      (draco-mode-update-menu (draco-menu-insert-comments-makefile)))
    (add-hook 'sh-mode-hook 'draco-sh-mode-hook)
    (add-hook 'sh-mode-hook 'turn-on-draco-mode)
    (add-hook 'sh-mode-hook 'turn-on-font-lock)
    (add-hook 'sh-mode-hook 'turn-on-auto-fill)))

;; ========================================
;; SGML mode
;; ========================================

(defun draco-setup-sgml-mode ()
  "Autoload sgml-mode, append the approriate suffixes to
auto-mode-alist and set up some customizations for DRACO."
  (interactive)
  (progn
    (autoload 'sgml-mode "sgml-mode" "SGML Shell Editing Mode" t)
    (if draco-colorize-modeline 
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
    (add-hook 'sgml-mode-hook 'turn-on-draco-mode)
    (add-hook 'sgml-mode-hook 'turn-on-font-lock)
    (add-hook 'sgml-mode-hook 'turn-on-auto-fill)))

;; ========================================
;; text mode
;; ========================================

(defun draco-setup-text-mode ()
  "Autoload text-mode, append the approriate suffixes to
auto-mode-alist and set up some customizations for DRACO."
  (interactive)
  (progn
    (autoload 'text-mode "text-mode" "Text Editing Mode" t)
    (if draco-colorize-modeline 
	(add-hook 'text-mode-hook        
		  '(lambda () 
		     (set-face-background 'modeline 
					  "wheat" (current-buffer))
		     (set-face-foreground 'modeline
					  "black"   (current-buffer)))))
    (setq auto-mode-alist
	  (append
	   '(("\\.text$" . text-mode)
	     ("\\.txt$"  . text-mode)
	     ("\\.log$"  . text-mode)
	     ("^README*" . text-mode)
	     ) auto-mode-alist))
    (add-hook 'text-mode-hook 'turn-on-draco-mode)
    (add-hook 'text-mode-hook 'turn-on-font-lock)
    (add-hook 'text-mode-hook 'turn-on-auto-fill)))

;; ========================================
;; Dired mode
;; ========================================

(defun draco-setup-dired-mode ()
  "Autoload dired-mode, append the approriate suffixes to
auto-mode-alist and set up some customizations for DRACO."
  (interactive)
  (progn
    (autoload 'dired-mode "dired-mode" "Dired File Services Mode" t)
    (add-hook 'dired-mode-hook 'turn-on-font-lock)
    (defun draco-dired-mode-hook ()
      "Hooks added to dired-mode"
      (local-set-key [(f5)] 'dired-redisplay-subdir)
      (if draco-colorize-modeline 
	  (add-hook 'dired-mode-hook        
		    '(lambda () 
		       (set-face-background 'modeline 
					    "thistle" (current-buffer))
		       (set-face-foreground 'modeline
					    "black"   (current-buffer))))))
    (add-hook 'dired-setup-keys-hook 'draco-dired-mode-hook)
    (add-hook 'dired-mode-hook 'turn-on-draco-mode)))

;; ========================================
;; Perl mode
;; ========================================

(defun draco-setup-perl-mode ()
  "Autoload dired-mode, append the approriate suffixes to
auto-mode-alist and set up some customizations for DRACO."
  (interactive)
  (progn
    (autoload 'perl-mode "perl-mode" "Perl Mode" t)
    (add-hook 'perl-mode-hook 'turn-on-font-lock)
    (defun draco-perl-mode-hook ()
      "Hooks added to perl-mode"
      (local-set-key [(f5)] 'dired-redisplay-subdir)
      (if draco-colorize-modeline 
	  (add-hook 'perl-mode-hook        
		    '(lambda () 
		       (set-face-background 'modeline 
					    "thistle" (current-buffer))
		       (set-face-foreground 'modeline
					    "black"   (current-buffer))))))
    (setq auto-mode-alist
          (append
           '(("\\.perl$" . perl-mode)
             ) auto-mode-alist))
    (add-hook 'perl-mode-hook 'turn-on-draco-mode)
    )
  )

;;---------------------------------------------------------------------------;;
;; Provide these functions from draco-config-modes
;;---------------------------------------------------------------------------;;

(provide 'draco-config-modes)
