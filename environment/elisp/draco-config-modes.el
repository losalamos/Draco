;; ======================================================================
;; draco-config-modes.el
;;
;; Copyright (C) 2016-2018 Los Alamos National Security, LLC
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

(defvar tcl-indent-level 2)
(defun draco-setup-tcl-mode ()
  "Autoload tcl-mode and append the appropriate suffixes to
auto-mode-alist."
  (interactive)
  (progn
    (autoload 'tcl "tcl" nil t)
    (setq auto-mode-alist
	  (append '(("\\.tcl$" . tcl-mode)
		    ("\\.itk$" . tcl-mode)
		    ("\\.ith$" . tcl-mode)
		    ("\\.itm$" . tcl-mode)
		    ) auto-mode-alist))
    (defun draco-tcl-mode-hook ()
      "draco-mode hooks added to TCL mode."
      (turn-on-draco-mode)
      (turn-on-auto-fill)
      (setq tcl-indent-level 2))
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
        (set-fill-column draco-code-comment-width)
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
      (turn-on-auto-fill)
      (set-fill-column draco-code-comment-width)
      (require 'fill-column-indicator)
      (fci-mode))
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
- Set fill-column to (draco-code-comment-width).
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
		      ("\\.cxx$"    . c++-mode)
                      ("\\.pt$"     . c++-mode)
		      ("\\.hh$"     . c++-mode)
		      ("\\.hpp$"    . c++-mode)
		      ("\\.cpp$"    . c++-mode)
                      ("\\.hh.in$"  . c++-mode)
		      ("\\.h.in$"   . c-mode)
		      ("\\.c$"      . c-mode)   ; to edit C code
		      ("\\.cu$"     . c-mode)   ; to edit CUDA kernels
		      ("\\.h$"      . c-mode)   ; to edit C code
		      ("\\.dcc$"    . c++-mode)   ; to edit C code
		      ("\\.dcc.in$" . c++-mode)   ; to edit C code
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
          ; Tab indent == 2 spaces
	  (c-basic-offset . 2)
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
			      ;(extern-lang-close . 0)
			      ;(extern-lang-open  . 0)
                              (inextern-lang     . 0)
			      (inline-close      . 0)
			      (inline-open       . 0)
			      (innamespace       . 0)
			      (statement-case-intro . +)
			      (statement-cont    . c-lineup-math)
			      (substatement-open . 0)
			      )
                           )
          )
         )
        )

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
- Sets fill-column to (draco-code-comment-width)
- Sets f5/f6 as hot keys to insert dividers.
- Turns on auto-fill"
	(draco-setup-c-mode)
	(c-set-style "draco")
	(local-set-key "\C-m" 'newline-and-indent)
	(set-fill-column draco-code-comment-width)
        (require 'fill-column-indicator)
        (fci-mode)
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

      ;(setq reftex-texpath-environment-variables draco-texpath)
      ;(setq reftex-bibpath-environment-variables draco-bibpath)

      (setq reftex-bibpath-environment-variables
            '("BIBINPUTS" "TEXBIB" "!kpsewhich -show-path=.bib"))

      (defun draco-latex-mode-hook ()
	"DRACO hooks added to LaTeX and BibTex modes."
	(local-set-key [(f5)] 'draco-latex-divider)
	(local-set-key [(f6)] 'draco-latex-comment-divider)
	(local-set-key "\C-c %" 'comment-region)
	(draco-mode-update-menu (draco-menu-insert-comments-tex))
	;(turn-on-bib-cite)
	(turn-on-reftex)
        (turn-on-reftex)
        (setq reftex-plug-into-AUCTeX t)
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
	     ("\\.F90$"  . f90-mode)
	     ("\\.FH$"   . f90-mode)
	     ("\\.fm4$"  . f90-mode)
	     ) auto-mode-alist))

    (defun draco-f90-mode-hook ()
      "Hooks added to F90 mode. See https://jblevins.org/log/f90-mode"
      (local-set-key [(f5)]         'draco-f90-subroutine-divider)
      (local-set-key [(control f6)] 'draco-f90-insert-document)
      (local-set-key [(f6)]         'draco-f90-comment-divider)
      (draco-mode-update-menu (draco-menu-insert-comments-f90))
      (set-fill-column draco-code-comment-width)
      (abbrev-mode 1)
      (setq f90-font-lock-keywords f90-font-lock-keywords-3)
      (setq f90-beginning-ampersand nil)
      (setq f90-associate-indent 0)
      (require 'fill-column-indicator)
      (fci-mode))
     ;; let .F denone Fortran and not freeze files
    (defvar crypt-freeze-vs-fortran nil)
    (add-hook 'f90-mode-hook 'draco-f90-mode-hook)
    (add-hook 'f90-mode-hook 'turn-on-draco-mode)
    (add-hook 'f90-mode-hook 'turn-on-auto-fill)
    ; should add this sometime
    ; (add-hook 'font-lock-mode-hook
    ;          '(lambda ()
    ;             (if (major-mode 'f90-mode)
    ;                 (draco-f90-font-lock)))) ; create this function
                                        ; based on draco-font-lock but
                                        ; for f90
    ))

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
    (add-hook 'fortran-mode-hook 'turn-on-auto-fill)
    ))

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
    (set-fill-column draco-code-comment-width)
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
    (require 'psvn)
    (setq svn-status-verbose nil)

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
    ))

;; ========================================
;; Doxymacs Mode
;; ========================================

(require 'doxymacs)
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
      (draco-mode-update-menu (draco-menu-insert-comments-makefile))
      (set-fill-column draco-code-comment-width)
      (setq sh-basic-offset 2 sh-indentation 2)
      )
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
      (setq dired-listing-switches "-alh")
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


;; ========================================
;; ECB & CEDET
;;
;; http://alexott.net/en/writings/emacs-devenv/EmacsCedet.html
;;
;; When installing cedet from trunk:
;; 1. run make from top level
;; 2. run make from the contrib directory.
;;
;; Options to add to ~/.emacs
;;
;; I set these from within emacs.  The most notable is probably
;; ecb-prescan-directories-for-emptyness, which was triggering
;; "Permission denied" errors when trying to walk directories on HPC
;; machines.
;; (custom-set-variables
;;  ;; custom-set-variables was added by Custom.
;;  ;; If you edit it by hand, you could mess it up, so be careful.
;;  ;; Your init file should contain only one such instance.
;;  ;; If there is more than one, they won't work right.
;;  '(ecb-layout-window-sizes (quote (("left8" (ecb-directories-buffer-name 0.3548387096774194 . 0.28888888888888886) (ecb-sources-buffer-name 0.3548387096774194 . 0.24444444444444444) (ecb-methods-buffer-name 0.3548387096774194 . 0.28888888888888886) (ecb-history-buffer-name 0.3548387096774194 . 0.16666666666666666)))))
;;  '(ecb-options-version "2.40")
;;  '(ecb-prescan-directories-for-emptyness nil)
;;  '(ecb-primary-secondary-mouse-buttons (quote mouse-1--mouse-2)))
;; ========================================
(defun draco-start-ecb ()
"Start Emacs Code Browser"
  (interactive)
  (if (> emacs-major-version 23) ;; need emacs-24 or newer.
      (progn
        (defvar cedetver "latest" "CEDET version.")
        (defvar draco-vendor-dir "/ccs/codes/radtran/vendors" "Draco vendor dir")
        ;; Find and set the draco vendor directory.
        (if (file-accessible-directory-p "/ccs/codes/radtran/vendors")
            (setq draco-vendor-dir "/ccs/codes/radtran/vendors"))
        ;; if a vendor directory exists on the local disk, use it
        ;; instead of the NFS mounted location.
        (if (file-accessible-directory-p "/var/tmp/vendors")
            (setq draco-vendor-dir "/var/tmp/vendors"))
        ;; HPC
        (if (file-accessible-directory-p "/usr/projects/draco/vendors")
            (setq draco-vendor-dir "/usr/projects/draco/vendors"))
        ;; Darwin
        (if (file-accessible-directory-p "/projects/opt/draco/vendors")
            (setq draco-vendor-dir "/projects/opt/draco/vendors"))

        ;;
        ;; CEDET
        ;;
        (if (file-accessible-directory-p
             (concat draco-vendor-dir
                     "/elisp/cedet-" cedetver "/lisp"))
            (progn
              (load-file
               (concat draco-vendor-dir "/elisp/cedet-" cedetver "/cedet-devel-load.el"))
              (load-file
               (concat draco-vendor-dir "/elisp/cedet-" cedetver "/contrib/cedet-contrib-load.el"))

              ;; Add further minor-modes to be enabled by semantic-mode.
              ;; See doc-string of `semantic-default-submodes' for other things
              ;; you can use here.  Set these modes before enabling Semantic.
              (add-to-list 'semantic-default-submodes
                           'global-cedet-m3-minor-mode t)
              (add-to-list 'semantic-default-submodes
                           'global-semantic-mru-bookmark-mode t)

              ;; Enable Semantic
              (semantic-mode 1)
              (require 'semantic/ia)

              ;; Enable GCC-specific support
              (require 'semantic/bovine/gcc)

              ;; Parse headers when idle
              (setq semantic-idle-work-update-headers-flag t)

              ;; Decrease idle threshold
              (setq semantic-idle-scheduler-idle-time 0.3)

              ;; Decrease idle work threshold
              ;; (setq semantic-idle-scheduler-work-idle-time 30)

              ;; Allow creation of ebrowse databases
              ;; (setq semanticdb-default-system-save-directory
              ;;       "~/.emacs.d/semanticdb")

              ;; Customize CEDET key bindings
              (defun my-cedet-bindings-hook ()
                (local-set-key [(meta return)] 'semantic-ia-complete-symbol-menu)
                (local-set-key [(control return)] 'semantic-ia-complete-symbol-menu)
                (local-set-key "\C-c?" 'semantic-ia-complete-symbol)

                (local-set-key "\C-c>" 'semantic-complete-analyze-inline)
                (local-set-key "\C-c=" 'semantic-decoration-include-visit)

                (local-set-key "\C-cj" 'semantic-ia-fast-jump)
                (local-set-key "\C-cq" 'semantic-ia-show-doc)
                (local-set-key "\C-ca" 'semantic-ia-show-summary)
                (local-set-key "\C-cp" 'semantic-analyze-proto-impl-toggle))
              (add-hook 'semantic-init-hooks 'my-cedet-bindings-hook)

              ;; Enable tag folding
              (defun my-semantic-folding-hook ()
                (semantic-tag-folding-mode 1)
                (local-set-key "\C-c-" 'semantic-tag-folding-fold-block)
                (local-set-key "\C-c+" 'semantic-tag-folding-show-block))
              (add-hook 'semantic-init-hooks 'my-semantic-folding-hook)

              ))

        ;;
        ;; Global gtags
        ;;
        (defvar gtagsglobalver "6.2.9" "Global gtags version.")
        (if (file-accessible-directory-p (concat
                                          draco-vendor-dir
                                          "/elisp/global-" gtagsglobalver "/share/gtags"))
            (progn
              (load-file (concat draco-vendor-dir "/elisp/global-"
                                 gtagsglobalver "/share/gtags/gtags.el"))
              (autoload 'gtags-mode "gtags" "" t)

              ;; Enable GNU GLOBAL support
              (when (cedet-gnu-global-version-check t)
                (semanticdb-enable-gnu-global-databases 'c-mode t)
                (semanticdb-enable-gnu-global-databases 'c++-mode t)
                (semanticdb-enable-gnu-global-databases 'f90-mode t))

              ;; Enable gtags in C mode and F90 mode
              (add-hook 'c-mode-common-hook '(lambda () (gtags-mode 1)))
              (add-hook 'f90-mode-hook '(lambda () (gtags-mode 1)))

              ;; Switch some bindings for C and F90 mode
              (defun my-gtags-bindings-hook ()
                (local-set-key "\M-."  'gtags-find-tag))
              (add-hook 'c-mode-common-hook 'my-gtags-bindings-hook)
              (add-hook 'f90-mode-hook 'my-gtags-bindings-hook)

              ))

        ;;
        ;; Eassist
        ;;
        (require 'eassist)

        ;; Teach eassist about .hh files
        (setq eassist-header-switches '(("h" . ("cc" "c"))
                                        ;; ("i.hh" . ("t.hh" "cc" "hh"))
                                        ;; ("t.hh" . ("cc" "hh" "i.hh"))
                                        ("hh" . ("i.hh" "t.hh" "cc"))
                                        ("cc" . ("hh" "i.hh" "t.hh"))
                                        ("H" . ("C" "CC"))
                                        ("c" . ("h"))
                                        ("C" . ("H"))))

        ;; Customize eassist key bindings
        (defun my-eassist-bindings-hook ()
          (local-set-key "\C-ct" 'eassist-switch-h-cpp)
          (local-set-key [(shift f8)] 'eassist-switch-h-cpp)
          (local-set-key "\C-ce" 'eassist-list-methods)
          (local-set-key "\C-c\C-r" 'semantic-symref))
        (add-hook 'c-mode-common-hook 'my-eassist-bindings-hook)

        ;;
        ;; EDE (Project Management)
        ;;
        (global-ede-mode 1)

        ;; This trips on "unsafe" operations in some projects
        ;; (ede-enable-generic-projects)

        ;; Tell EDE to use Global to locate files
        (setq ede-locate-setup-options
              '(ede-locate-global ede-locate-base))

        ;; Set up an EDE project for Draco+Jayenne
        ;; If the user has not set
        ;; file-at-root-level-draco, then set it to
        ;; ~/.draco_ede (this assumes that draco and jayenne are checked
        ;; out at $HOME).
        (if (not (boundp 'file-at-root-level-draco))
            (defvar file-at-root-level-draco "~/.draco_ede" ))
        ;; If the file does not exist, create it as an empty file.
        (if (not (file-exists-p file-at-root-level-draco))
            (write-region "" nil file-at-root-level-draco))

        ;; Link Draco, Jayenne and Capsaicin sources so semantic can find
        ;; all the sources.
        (ede-cpp-root-project "draco"
                              :name "Draco Project"
                              :file file-at-root-level-draco
                              :include-path '("/draco/src"
                                              "/jayenne/clubimc/src"
                                              "/jayenne/wedgehog/src"
                                              "/jayenne/milagro/src"
                                              "/capsaicin/src")
                              ;; :system-include-path '("~/exp/include")
                              ;; :spp-table '(("isUnix" . "")
                              ;;              ("BOOST_TEST_DYN_LINK" . "" )
                              )

        ;; Set up an EDE project for EAP

        ;; If the user has not set file-at-root-level-eap, then set it
        ;; to ~/cassio/.eap_ede
        (if (not (boundp 'file-at-root-level-eap ))
            (defvar file-at-root-level-eap "~/cassio/.eap_ede" ))
        ;; If the file does not exist, create it as an empty file.
        (if (not (file-exists-p file-at-root-level-eap))
            (if (file-accessible-directory-p "~/cassio" )
                (write-region "" nil file-at-root-level-eap)))
        ;; If the directory does not exist the above will fail.
        (if (file-exists-p file-at-root-level-eap)
            (ede-cpp-root-project "EAP"
                              :name "EA Project"
                              :file file-at-root-level-eap
                              :include-path
                              '("/Source.othello"
                                "/Source.othello/Draco"
                                "/Source.othello/IMC"
                                "/Source.othello/Sn"
                                "/Source.othello/cassio"
                                "/Source.rh"
                                "/Source.rh/Analysis"
                                "/Source.rh/Bibliography"
                                "/Source.rh/Comm"
                                "/Source.rh/Comm/Dummy_mpi"
                                "/Source.rh/Common"
                                "/Source.rh/Common_Util"
                                "/Source.rh/EOS"
                                "/Source.rh/ExternalSRCs"
                                "/Source.rh/Graphics"
                                "/Source.rh/Graphics/pv"
                                "/Source.rh/Graphics/pv/Adaptor"
                                "/Source.rh/Graphics/pv/Examples"
                                "/Source.rh/Gravity"
                                "/Source.rh/HEBurn"
                                "/Source.rh/HEBurn/Contours"
                                "/Source.rh/HEBurn/Contours/bin"
                                "/Source.rh/Hydro"
                                "/Source.rh/Hydro/OneFileHydros"
                                "/Source.rh/Hydro/OpenCLHydro"
                                "/Source.rh/Hydro/OpenCLHydro/Unit"
                                "/Source.rh/IO"
                                "/Source.rh/Iso"
                                "/Source.rh/Laser"
                                "/Source.rh/MMS"
                                "/Source.rh/MaterialInterface"
                                "/Source.rh/MaterialInterface/VOFQ"
                                "/Source.rh/Mesh"
                                "/Source.rh/Mesh/gem_old_test"
                                "/Source.rh/OpenCLUtil"
                                "/Source.rh/OpenCLUtil/Unit"
                                "/Source.rh/Parser"
                                "/Source.rh/Plasma"
                                "/Source.rh/Plasma/tests"
                                "/Source.rh/Radiation"
                                "/Source.rh/Rage"
                                "/Source.rh/Roxane"
                                "/Source.rh/Roxane/OpenCL"
                                "/Source.rh/Roxane/OpenCL/Unit"
                                "/Source.rh/Roxane_Util"
                                "/Source.rh/Setup"
                                "/Source.rh/Setup/SpicaCSG"
                                "/Source.rh/Solvers"
                                "/Source.rh/Solvers/dierckx"
                                "/Source.rh/Strength"
                                "/Source.rh/Strength/TEPLA"
                                "/Source.rh/TNBurn"
                                "/Source.rh/Turbulence"
                                "/Source.rh/Util"
                                "/Source.rh/Util_basic"
                                "/Source.rh/build"
                                "/Source.rh/dump_reader"
                                "/Source.rh/l7"
                                "/Source.rh/l7/config"
                                "/Source.rh/l7/libsrc"
                                "/Source.rh/l7/libsrc/l7"
                                "/Source.rh/l7/tests"
                                "/Source.rh/l7/tests/F7test"
                                "/Source.rh/l7/tests/L7test"
                                "/Source.rh/leLinkLib"
                                "/Source.rh/leLinkLib/leLinkLib"
                                "/Source.rh/leLinkLib/sageFiles"
                                "/Source.rh/m4"
                                "/Source.rh/xRage")
                              :header-match-regexp
                              "\\.\\(h\\(h\\|xx\\|pp\\|\\+\\+\\)?\\|H\\|f90\\)$\\|\\<\\w+$")
        ) ;; endif

        ;;
        ;; ECB
        ;;
        (defvar ecbver "latest" "Version of Emacs Code Browser.")
        (if (file-accessible-directory-p (concat draco-vendor-dir "/elisp/ecb-" ecbver))
            (progn
              ; prevents ecb failing on start
              (setq ecb-version-check nil)
              (defadvice ecb-check-requirements (around no-version-check activate compile)
                "ECB version checking code is very old so that it thinks that the latest
cedet/emacs is not new enough when in fact it is years newer than the latest
version that it is aware of. So simply bypass the version check."
                (if (or (< emacs-major-version 23)
                        (and (= emacs-major-version 23)
                             (< emacs-minor-version 3)))
                    ad-do-it))

              (add-to-list 'load-path (concat draco-vendor-dir
                                              "/elisp/ecb-" ecbver))
              ; (require 'ecb-autoloads) ;; only for old version of ecb
              (setq ecb-tip-of-the-day nil)
              (setq ecb-source-path '(my-home-dir))

              (require 'ecb)
              ;; M-x ecb-activate
              ;; M-x ecb-byte-compile
              ;; M-x ecb-show-help
              ))
        )
    (message "CEDT/ECB setup required emacs version 24 or newer.")
    )
  )

;;---------------------------------------------------------------------------;;
;; Provide these functions from draco-config-modes
;;---------------------------------------------------------------------------;;

(provide 'draco-config-modes)
