;;; tme-hacks.el
;;; Tom Evans
;;; Feb 26, 1998

;;; $Id$
;;; $Name$

;;---------------------------------------------------------------------------;;
;; Provide additional functionality to the rtt-hacks.el package
;;---------------------------------------------------------------------------;;

;;---------------------------------------------------------------------------;;

;; setup some nifty comment stuff

;;---------------------------------------------------------------------------;;

(defun tme-f90-subroutine-divider ()
  (interactive)
  (beginning-of-line)
  (insert "!-----------------------------------------------------------------------------!\n")
  (insert "! \n")
  (insert "!-----------------------------------------------------------------------------!\n\n")
  (previous-line 3)
  (end-of-line)
)

(defun tme-f90-comment-divider ()
  (interactive)
  (beginning-of-line)
  (insert "!-----------------------------------------------------------------------------!\n")
  (end-of-line)
)

(defun tme-f77-subroutine-divider ()
  (interactive)
  (beginning-of-line)
  (insert "c-----------------------------------------------------------------------------c\n")
  (insert "c \n")
  (insert "c-----------------------------------------------------------------------------c\n\n")
  (previous-line 3)
  (end-of-line)
)

(defun tme-f77-comment-divider ()
  (interactive)
  (beginning-of-line)
  (insert "c-----------------------------------------------------------------------------c\n")
  (end-of-line)
)    

(defun tme-memo-dist ()
  (interactive)
  (beginning-of-line)
  (insert "CCS-4  MS D409:\\\\ \n")
  (insert "CCS-2  MS D413:\\\\ \n")
  (insert "CCS-DO MS B297:\\\\ \n")
  (insert "Feiereisen, Bill,  CCS-DO MS B297\\\\ \n")
  (insert "Lee, Stephen,      CCS-DO MS B297\\\\ \n")
  (insert "Baker, Randal,     CCS-4  MS D409\\\\ \n")
  (insert "Budge, Kent,       CCS-4  MS D409\\\\ \n")
  (insert "Buksas, Michael,   CCS-4  MS D409\\\\ \n")
  (insert "Carrington, David, CCS-4  MS D409\\\\ \n")
  (insert "Dahl, Jon,         CCS-4  MS D409\\\\ \n")
  (insert "Densmore, Jeffery, CCS-4  MS D409\\\\ \n")
  (insert "Evans, Thomas,     CCS-4  MS D409\\\\ \n")
  (insert "Hungerford, Aimee, CCS-4  MS D409\\\\ \n")
  (insert "Olson, Gordon,     CCS-4  MS D409\\\\ \n")
  (insert "Thompson, Kelly,   CCS-4  MS D409\\\\ \n")
  (insert "Turner, Scott,     CCS-4  MS D409\\\\ \n")
  (insert "Urbatsch, Todd,    CCS-4  MS D409\\\\ \n")
  (insert "Ward, Robert,      CCS-4  MS B296\\\\ \n")
  (insert "Warsa, James,      CCS-4  MS D409\\\\ \n")
  (insert "Dendy, Edward,     CCS-2  MS D413\\\\ \n")
  (insert "Dilts, Gary,       CCS-2  MS D413\\\\ \n")
  (insert "Hall, Michael,     CCS-2  MS D413\\\\ \n")
  (insert "Kothe, Douglas,    CCS-2  MS D413\\\\ \n")
  (insert "Lowrie, Robert,    CCS-2  MS D413\\\\ \n")
  (insert "Morel, Jim,        CCS-2  MS D413\\\\ \n")
  (insert "Rider, William,    CCS-2  MS D413\\\\ \n")
  (insert "Turner, John,      CCS-2  MS D413\\\\ \n")
  (insert "Wingate, Beth,     CCS-2  MS D413\\\\ \n")
  (end-of-line)
)

(defun tme-latex-divider ()
  (interactive)
  (beginning-of-line)
  (insert "%%---------------------------------------------------------------------------%%\n")
  (insert "%% \n")
  (insert "%%---------------------------------------------------------------------------%%\n\n")
  (previous-line 3)
  (end-of-line)
)

(defun tme-latex-comment-divider ()
  (interactive)
  (beginning-of-line)
  (insert "%%---------------------------------------------------------------------------%%\n")
  (end-of-line)
)

(defun tme-makefile-divider ()
  (interactive)
  (beginning-of-line)
  (insert "##---------------------------------------------------------------------------##\n")
  (insert "## \n")
  (insert "##---------------------------------------------------------------------------##\n\n")
  (previous-line 3)
  (end-of-line)
)

(defun tme-m4-divider ()
  (interactive)
  (beginning-of-line)
  (insert "dnl-------------------------------------------------------------------------dnl\n")
  (insert "dnl \n")
  (insert "dnl-------------------------------------------------------------------------dnl\n\n")
  (previous-line 3)
  (end-of-line)
)

(defun tme-m4-comment-divider ()
  (interactive)
  (beginning-of-line)
  (insert "dnl-------------------------------------------------------------------------dnl\n")
  (end-of-line)
)

(defun tme-makefile-comment-divider ()
  (interactive)
  (beginning-of-line)
  (insert "##---------------------------------------------------------------------------##\n")
  (end-of-line)
)

(defun tme-c-comment-divider ()
  (interactive)
  (beginning-of-line)
  (insert "/*---------------------------------------------------------------------------*/\n")
  (end-of-line)
)

(defun tme-html-comment-divider ()
  (interactive)
  (beginning-of-line)
  (insert "<!---------------------------------------------------------------------------->\n")
  (end-of-line)
)

(defun tme-c-divider ()
  (interactive)
  (beginning-of-line)
  (insert "/*---------------------------------------------------------------------------*/\n")
  (insert "/* \n")
  (insert "/*---------------------------------------------------------------------------*/\n\n")
  (previous-line 3)
  (end-of-line)
)

(defun tme-elisp-comment-divider ()
  (interactive)
  (beginning-of-line)
  (insert ";;---------------------------------------------------------------------------;;\n")
  (end-of-line)
)

(defun tme-elisp-divider ()
  (interactive)
  (beginning-of-line)
  (insert ";;---------------------------------------------------------------------------;;\n")
  (insert ";; \n")
  (insert ";;---------------------------------------------------------------------------;;\n\n")
  (previous-line 3)
  (end-of-line)
)

;;---------------------------------------------------------------------------;;

;; Define some mode-hooks.

;; use f5 and f6 for comment dividers and such in different modes
;; Also configure tabs and font-lock for some modes.

;; when adding hooks the syntax is:
;;  (add-hook 'mode-hook 'tme-function) if adding a function
;;  (add-hook 'mode-hook tme-variable) if adding a variable

(defun tme-latex-mode-hook ()
  "Hooks added to LaTeX-mode"
  (local-set-key [(f5)] 'tme-latex-divider)
  (local-set-key [(f6)] 'tme-latex-comment-divider))

(add-hook 'LaTeX-mode-hook 'tme-latex-mode-hook t)

(defun tme-bibtex-mode-hook ()
  "Hooks added to BiBTeX mode"
  (local-set-key [(f5)] 'tme-latex-divider)
  (local-set-key [(f6)] 'tme-latex-comment-divider)
  (auto-fill-mode))

(add-hook 'bibtex-mode-hook 'tme-bibtex-mode-hook t)

(defun tme-makefile-mode-hook ()
  "Hooks added to Makefile mode"
  (local-set-key [(f5)] 'tme-makefile-divider)
  (local-set-key [(f6)] 'tme-makefile-comment-divider))
 
(defun tme-c-mode-hook ()
  "Hooks added to c-mode"
  (local-set-key [(f5)] 'tme-c-divider)
  (local-set-key [(f6)] 'tme-c-comment-divider))

(defun tme-autoconf-mode-hook ()
  "Hooks added to autoconf mode"
  (setq tab-stop-list '(3 7 11 15 19 23 27 31 35 39 43 47 51 55 59 63 67 71 75 79 83))
  (local-set-key [(f9)] 'tme-m4-divider)
  (local-set-key [(tab)] 'tab-to-tab-stop)
  (local-set-key [(f10)] 'tme-m4-comment-divider))

(add-hook 'c-mode-hook 'tme-c-mode-hook t)
(add-hook 'sh-mode-hook 'tme-makefile-mode-hook t)
(add-hook 'autoconf-mode-hook 'tme-makefile-mode-hook t)
(add-hook 'autoconf-mode-hook 'tme-autoconf-mode-hook t)

(defun tme-f90-mode-hook ()
  "Hooks added to F90 mode"
  (local-set-key [(f5)] 'tme-f90-subroutine-divider)
  (local-set-key [(control f6)] 'tme-f90-insert-document)
  (local-set-key [(f6)] 'tme-f90-comment-divider)
  (set-fill-column 80))

(add-hook 'f90-mode-hook 'tme-f90-mode-hook t)

(defun tme-f77-mode-hook ()
  "Hooks added to F77 mode"
  (local-set-key [(f5)] 'tme-f77-subroutine-divider)
  (local-set-key [(control f6)] 'tme-f77-insert-document)
  (local-set-key [(f6)] 'tme-f77-comment-divider))

(add-hook 'fortran-mode-hook 'tme-f77-mode-hook t)

(defun tme-elisp-mode-hook ()
  "Hooks added to Elisp mode"
  (local-set-key [(f5)] 'tme-elisp-divider)
  (local-set-key [(f6)] 'tme-elisp-comment-divider))

(add-hook 'emacs-lisp-mode-hook 'tme-elisp-mode-hook t)

(defun tme-changelog-mode-hook ()
  "Hooks added to ChangeLog mode"
  (turn-on-auto-fill))

(add-hook 'change-log-mode-hook 'tme-changelog-mode-hook t)

;;---------------------------------------------------------------------------;;

;;; setup tme-hacks

(defun tme-setup ()
  (message "Configuring TME support."))

(provide 'tme-hacks)

;;---------------------------------------------------------------------------;;
;; end of tme-hacks.el
;;---------------------------------------------------------------------------;;
