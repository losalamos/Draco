;; lucid-setup.el

(message "working on xemacs-setup")

;;---------------------------------------------------------------------------;;
; Make sure that the info package is loaded before we append to
; Info-directory-list. 

(require 'info)

(defvar ccs4-info-dirs
  (list 
   (concat my-home-dir ".xemacs/info/")
   "/codes/radtran/vendors/gsl/Linux/info/"
   "/usr/local/autotools_2.13/info/"
   "/usr/local/gcc-3.4.2/info/"
   "/usr/local/cvs1.12/info/"
   "/usr/info") 
  "\nA list of directories that may contain GNU Info files.
\nAppend to list by adding the following command in personal elisp script:
\n\t (setq ccs4-info-dirs (cons \"/full/path/to/info/dir\") ")

;; Update Info-directory-list using directories provided by ccs4-info-dirs.
(defun ccs4-append-info-dirs (ldirs)
  (interactive)
  (let ((dirlist ldirs))
    (progn
      ;; process least-imporant directories first.  This ensures that
      ;; the first directory in the list "ccs4-env-dirs" become the
      ;; first directory in the Info-directory-list.
      (setq dirlist (reverse dirlist))
      ;; Loop over all directories in the list provided.  Add valid
      ;; elisp directories to the load path (prepend).
      (while dirlist
	(setq ldir (car dirlist))
	(if (file-accessible-directory-p ldir)
	    (setq Info-directory-list (cons ldir Info-directory-list)))
	;; remove ldir from dirlist and continue the while-loop.
	(setq dirlist (cdr-safe dirlist)))
      ))
  )
(ccs4-append-info-dirs ccs4-info-dirs)

;;---------------------------------------------------------------------------;;

;; Stuff for automatic creation of C++ files and classes.

;(defconst gmf-template-dir
;  (expand-file-name "~furnish/lisp/templates/")
;  "The directory where you keep your template files." )

;(defun gmf-new-cc-file (name)
;  "Function to spontaneously setup a new C++ file."
;  (interactive "sBase Name: ")

;  (setq nfile (concat name ".cc"))
;  (setq tfile (concat gmf-template-dir "template.cc"))
;  (if (file-exists-p nfile)
;      (error "File %s already exists." nfile))

;  (find-file nfile)
;  (insert-file-contents tfile)

;; Now perform customizations

;  (goto-char (point-min))
;  (perform-replace "<basename>" name nil nil nil )
;  (goto-char (point-min))
;  (perform-replace "<date>" (current-time-string) nil nil nil )

;  ; Now, hack out the <start> token, since we don't know what to do
;  ; for a genearal file.
;  (perform-replace "<start>" "" nil nil nil )
;)

;;---------------------------------------------------------------------------;;
;(defun gmf-new-class (name)
;  (interactive "sClass name: ")

;; Make sure we can create the new files.

;  (setq  hfile (concat name ".h"))
;  (setq ccfile (concat name ".cc"))
;  (if (or (file-exists-p  hfile)
;	  (file-exists-p ccfile))
;      (error "Cannot create class %s, files already exist." name))

;; Now locate the template files.

;  (setq  thfile (concat gmf-template-dir "template.h"))
;  (setq tccfile (concat gmf-template-dir "template.cc"))
;  (if (not (and (file-exists-p  thfile)
;		(file-exists-p tccfile)))
;      (error "Cannot access class template files."))

;; Now try to copy the template files over.

;  (copy-file  thfile  hfile)
;  (copy-file tccfile ccfile)

;; Now load the header and customize.

;  (find-file hfile)
;  (perform-replace "<class>" name nil nil nil )
;  (goto-char (point-min))
;  (perform-replace "<date>" (current-time-string) nil nil nil )

;; Now load the cc file and customize.

;  (find-file ccfile)
;  (perform-replace "<basename>" name nil nil nil )
;  (goto-char (point-min))
;  (perform-replace "<date>" (current-time-string) nil nil nil )
;  (goto-char (point-min))

;  ; Now, change the <start> token to include the header.
;  (perform-replace "<start>" "" nil nil nil )
;  (insert "#include \"" name ".h\"\n\n")
;)
;;;---------------------------------------------------------------------------;;

(defun gmf-new-itcl-method (name)
  "Construct a new [incr Tcl] method."
  (interactive "sMethod name: ")
  (insert "#-----------------------------------------------------------------------------#\n")
  (insert "# " name " - \n")
  (insert "#-----------------------------------------------------------------------------#\n")
  (insert "\n")
  (insert "    method " name " {} {\n")
  (insert "    }\n")
  (insert "\n")
  (previous-line 6)
  (end-of-line)
  )

(defun rtt-previous-buffer ()
  "Switch to previous buffer."
  (interactive)
  (switch-to-buffer (other-buffer (buffer-name))))

; Some functions for positioning the cursor in the other window.  I
; mainly wanted gmf-bot-other-window for use with compile mode, and
; did the top variant for orthogonality.

(defun gmf-bot-other-window ()
  (interactive)
  (other-window 1)
  (end-of-buffer))

(defun gmf-top-other-window ()
  (interactive)
  (other-window 1)
  (goto-char (point-min)))

(defun gmf-pan-screen-up ()
  (interactive)
  (scroll-up 1))

(defun gmf-pan-screen-down ()
  (interactive)
  (scroll-down 1))

(defun gmf-scroll-down ()
  (interactive)
  (setq wline (count-lines (window-start)
			   (save-excursion (beginning-of-line) (point))))
  (scroll-down)
  (move-to-window-line wline))

(defun rtt-save-and-kill ()
  "Save buffer and then kill it."
  (interactive)
  (if (buffer-file-name (current-buffer))
      (save-buffer))
  (kill-buffer (buffer-name)))

(require 'scroll-in-place)

;; Set up automatic buffer variables

(defun my-c-mode-common-hook ()
  (define-key c++-mode-map 'button3 'kill-region)
  )

(add-hook 'c-mode-hook		'turn-on-font-lock)
(add-hook 'c-mode-common-hook   'my-c-mode-common-hook)
(add-hook 'c++-mode-hook	'turn-on-font-lock)
(add-hook 'dired-mode-hook	'turn-on-font-lock)
(add-hook 'sh-mode-hook	        'turn-on-font-lock)

(add-hook 'mppl-mode-hook         'turn-on-font-lock)
(add-hook 'tcl-mode-hook          'my-tcl-mode-hook)
(add-hook 'latex-mode-hook        'turn-on-font-lock)
(add-hook 'fortran-mode-hook      'turn-on-font-lock)
(add-hook 'm4-mode-hook           'turn-on-font-lock)
(add-hook 'perl-mode-hook         'turn-on-font-lock)
(add-hook 'Newsgroup-mode-hook    'turn-on-font-lock)
(add-hook 'LaTeX-mode-hook        'turn-on-font-lock)
(add-hook 'bibtex-mode-hook       'turn-on-font-lock)
(add-hook 'f90-mode-hook          'turn-on-font-lock)
(add-hook 'autoconf-mode-hook     'turn-on-font-lock)

(defun my-tcl-mode-hook ()
  (turn-on-font-lock)
  (tcl-auto-fill-mode)
  (define-key tcl-mode-map 'button3 'kill-region)
  (define-key tcl-mode-map [(f5)] 'gmf-new-itcl-method)
)

;;; fast-lock is a package which speeds up the highlighting of files
;;; by saving information about a font-locked buffer to a file and
;;; loading that information when the file is loaded again.  This
;;; requires a little extra disk space be used.
;;;
;;; Normally fast-lock puts the cache file (the filename appended with
;;; .flc) in the same directory as the file it caches.  You can
;;; specify an alternate directory to use by setting the variable
;;; fast-lock-cache-directories.

(add-hook 'font-lock-mode-hook 'turn-on-fast-lock)

(require 'func-menu)
(define-key global-map 'f12 'function-menu)
(add-hook 'find-file-hooks 'fume-add-menubar-entry)

(require 'infer-mode)

(defun mjl-insert-function-doc ()
  (interactive)
  (beginning-of-line)
  (insert "//---------------------------------------------------------------------------//\n")
  (insert "// \n")
  (insert "//---------------------------------------------------------------------------//\n\n")
  (previous-line 3)
  (end-of-line)
)

(defun mjl-insert-C-function-doc ()
  (interactive)
  (insert "/*--------------------------------------------------------------------------*\\\n")
  (insert " * \n")
  (insert " *\n")
  (insert " * \n")
  (insert "\\*--------------------------------------------------------------------------*/\n")
  (previous-line 4)
  (end-of-line)
)

(defun mjl-insert-class-doc ()
  (interactive)
  (insert "//===========================================================================//\n")
  (insert "// class \n")
  (insert "//===========================================================================//\n")
  (previous-line 2)
  (end-of-line)
)

(defun mjl-insert-comment-divider ()
  (interactive)
  (beginning-of-line)
  (insert "//---------------------------------------------------------------------------//\n")
)

; ps-print customization stuff.

(require 'ps-print)

(global-set-key 'print 'ps-print-buffer-with-faces)
(global-set-key '(shift print) 'ps-print-region-with-faces)

; Go pick up my font-lock-keyword definitions.

(require 'fl-keywords)


