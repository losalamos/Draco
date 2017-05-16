;; -*-emacs-lisp-*-
;;======================================================================
;; file   : draco-mode.el
;; Author : Kelly Thompson
;; Created: 30 Nov 2004
;; Brief  : ELisp package for making Draco related stuff easier.
;; Note   : Copyright (C) 2016-2017 Los Alamos National Security, LLC.
;;======================================================================
;; Commentary:
;;
;; - Customize the variable ccs4-env-dirs
;;
;; - Put (require 'draco-mode) in your ~/.xemacs/init.el
;;
;; - Invoke draco-mode with M-x draco-mode.  To have draco-mode
;;   invoked automatically when in C/C++ mode, put
;;
;;   (add-hook 'c-mode-common-hook 'draco-mode)
;;
;;   in your .xemacs/init.el.
;;
;; - If you want Draco keywords fontified use M-x draco-font-lock.
;;   To do it automatically, add the following to your .xemacs/init.el:
;;
;;   (defun my-draco-font-lock-hook ()
;;     (if (or (eq major-mode 'c-mode) (eq major-mode 'c++-mode))
;;         (draco-font-lock)))
;;   (add-hook 'font-lock-mode-hook 'my-draco-font-lock-hook)
;;
;;   This will add the Draco keywords to c-mode and c++-mode only.
;;============================================================
;; Front matter and variables

(defconst draco-mode-version "1.0.0"
  "Draco-mode version number")

(defun draco-mode-version ()
  "Report the current version of draco-mode in the minibuffer."
  (interactive)
  (message "Using draco-mode version %s" draco-mode-version))

(defgroup draco-mode nil
  "Find documentation created by draco-mode, and create draco-mode comments."
  :group 'tools)

(defcustom darco-hooks nil
  "Draco-mode hooks."
  :group 'draco-mode)

(defcustom draco-companion-file-alist nil
  "A list of possible file pairs.  The file type of the current buffer
is located within the list (first regexp).  The filename extension
is loaded from the 2nd part of the pair and a filename with this
new extension is loaded (if it exists)."
  :type 'list
  :group 'draco-mode)

(defcustom want-pooma-style-by-default  nil
"Do you want to use alternate POOMA-style indentation?"
:group 'draco-mode
:type '(radio (const :tag "No" nil)
	      (const :tag "Yes" t)))

(defcustom want-draco-menu                t
"\nIf non-nil, an RTT menu will be placed in the XEmacs menubar.
This menu includes selections to create new files from RTT templates
and selections for frequently used XEmacs features (CVS, compile,
etc.)"
:group 'draco-mode
:type '(radio	(const :tag "Yes" t)
		(const :tag "No"  nil)))

(defcustom draco-mode-style
  "Qt"
  "The style of code indentation and keywords.

Must be one of, \"Qt\" or \"C++\". Setting this variable
to anything else will generate errors."
  :type '(radio	(const :tag "Qt" "Qt")
		(const :tag "C++" "C++"))
  :group 'draco-mode)

(defcustom draco-code-comment-width
  92
  "The number of colums wide to use for inserted comment blocks."
  :group 'draco-mode)

;;============================================================
;; End of customizable variables
;;============================================================


;;============================================================
;; Keymap
;;============================================================

(defvar draco-mode-map (make-sparse-keymap)
  "Keymap for draco minor mode.")

;; XEmacs key binding annoyances:
(define-key draco-mode-map "\C-x\C-k" 'kill-buffer)
(define-key draco-mode-map "\e?"      'help-for-help)
(define-key draco-mode-map "\C-x?"    'describe-key-briefly)

;; XEmacs Accessories
(define-key draco-mode-map [(f11)]    'grep)
;; support for F11 on sunray 2F
(define-key draco-mode-map [(SunF36)] 'grep)
(define-key draco-mode-map [(control f11)] 'speedbar)
(define-key draco-mode-map "\C-csb"        'speedbar)

;; CVS
(define-key draco-mode-map [(f9)]               'svn-examine)
(define-key draco-mode-map [(meta f9)]          'svn-status)
(define-key draco-mode-map [(control f9)]       'svn-checkout)

;; Mouse commands
(define-key draco-mode-map [(button3)]          'kill-region)
(define-key draco-mode-map [(meta button3)]     'delete-rectangle)
(define-key draco-mode-map [(control button1)]  'popup-buffer-menu)
(define-key draco-mode-map [(control button2)]  'function-menu)
(define-key draco-mode-map [(control button3)]  'popup-mode-menu)

;; Keypad commands
(define-key draco-mode-map [(kp-multiply)]       'start-kbd-macro)
(define-key draco-mode-map [(kp-subtract)]       'end-kbd-macro)
(define-key draco-mode-map [(kp-add)] 	         'call-last-kbd-macro)
;(define-key draco-mode-map [(kp-divide)] 	 'byte-compile-file)
;(define-key draco-mode-map [(kp-enter)] 	 'other-window)
;(define-key draco-mode-map [(control kp-enter)] 'gmf-top-other-window)
;(define-key draco-mode-map [(shift kp-enter)]   'gmf-bot-other-window)
;(define-key draco-mode-map [(kp-7)] 	         'font-lock-fontify-buffer)
;(define-key draco-mode-map [(delete)] 	         'delete-char)
;(define-key draco-mode-map [(hpDeleteChar)]     'delete-char)
;(define-key draco-mode-map [(shift   delete)]   'delete-char)
;(define-key draco-mode-map [(control delete)]   'delete-char)
;(define-key draco-mode-map [(delete)]           'delete-backward-char)
;(define-key draco-mode-map [(prior)]            'scroll-down)
;(define-key draco-mode-map [(next)]             'scroll-up)

;; Compiling

(define-key draco-mode-map [(f1)]              'compile)
(define-key draco-mode-map [(control meta g)]  'what-line)
(define-key draco-mode-map [(f3)]              'previous-error)
(define-key draco-mode-map [(f4)]              'next-error)

;; Insert ChangeLog comment.
(define-key draco-mode-map [(shift control f5)] 'add-change-log-entry)

;; Buffer management

(define-key draco-mode-map [(f7)]               'draco-save-and-kill-current-buffer)
(define-key draco-mode-map [(shift f7)]         'delete-window)
(define-key draco-mode-map [(control f7)]       'kill-this-buffer)
(define-key draco-mode-map [(f8)]               'draco-toggle-previous-buffer)
(define-key draco-mode-map [(control f8)]       'draco-find-companion-file)

; Some bindings to drop Brief style bookmarks.
; Shift Fx sets bookmakr x, Control Fx returns to bookmark x.

;(define-key draco-mode-map [(meta f1)]     [(control x) r (space) a])
;(define-key draco-mode-map [(control f1)]  [(control x) r j a])
;(define-key draco-mode-map [(meta f2)]     [(control x) r (space) b])
;(define-key draco-mode-map [(control f2)]  [(control x) r j b])
;(define-key draco-mode-map [(meta f3)]     [(control x) r (space) c])
;(define-key draco-mode-map [(control f3)]  [(control x) r j c])
;(define-key draco-mode-map [(meta f4)]     [(control x) r (space) d])
;(define-key draco-mode-map [(control f4)]  [(control x) r j d])

(define-key draco-mode-map [(control f2)] 'font-lock-fontify-buffer)

;;============================================================
;; autoload
;;============================================================

;;;###autoload
(defun turn-on-draco-mode (&optional ignored)
  "Turn on draco minor mode unconditionally."
  (interactive)
  (draco-mode 1))

;;;###autoload
(defun turn-off-draco-mode (&optional ignored)
  "Turn off draco minor mode unconditionally."
  (interactive)
  (draco-mode 0))

;;;###autoload
(or (assoc 'draco-mode minor-mode-alist)
    (setq minor-mode-alist
	  (cons '(draco-mode " draco") minor-mode-alist)))

;;;###autoload
(or (assoc 'draco-mode minor-mode-map-alist)
    (setq minor-mode-map-alist
	  (cons (cons 'draco-mode draco-mode-map)
		minor-mode-map-alist)))

;;============================================================
;; This stuff has to do with fontification
;;
;; For RegEx help see http://dp.iit.bme.hu/mosml/doc/telepites-emacs-sml.txt
;;============================================================

(require 'font-lock)
(make-face 'font-lock-draco-dbc-face)

(defconst draco-mode-keywords
  (list
   (list
    ;; DOXYGEN - One shot keywords that take no arguments
    (concat "\\([@\\\\]\\(brief\\|li\\|\\(end\\)?code\\|sa"
	    "\\|note\\|\\(end\\)?verbatim\\|return\\|arg\\|fn"
	    "\\|hideinitializer\\|showinitializer\\|interface"
	    "\\|internal\\|nosubgrouping\\|author\\|date\\|endif"
	    "\\|invariant\\|post\\|pre\\|remarks\\|since\\|test\\|version"
	    "\\|\\(end\\)?htmlonly\\|\\(end\\)?latexonly\\|f\\$\\|file"
	    "\\|mainpage\\|name\\|overload\\|typedef\\|deprecated\\|par"
	    "\\|addindex\\|line\\|skip\\|skipline\\|until\\|see"
	    "\\|endlink\\)\\)\\>")
    '(0 font-lock-keyword-face prepend))
   ;; attention, warning, etc. given a different font
   (list
    "\\([@\\\\]\\(attention\\|warning\\|todo\\|bug\\)\\)\\>"
    '(0 font-lock-preprocessor-face prepend))

   ;; DRACO DBC - One shot keywords that take no arguments but are
   ;; followed by parenthesis.
   (list
    ;; Match single keyword that is followed by 0 or more spaces,
    ;; followed by an opening paren.
    (concat "\\(Assert\\|Check\\|Ensure\\|Insist\\|Require"
            "\\|Remember\\|UT_EPILOG\\)\\>\\(\\s-*[(]\\)")
    '(1 font-lock-preprocessor-face prepend))

   ;; DRACO Unit testing - One shot keywords that take no arguments
   (list
    ;; Match single keyword
    "\\(FAILMSG\\|ITFAILS\\|PASSMSG\\)\\>"
    '(0 font-lock-preprocessor-face prepend))

   ;; Draco extra C++ typenames
   ;; (list
   ;;  ;; Match single keyword that is followed by at whitespace and
   ;;  ;; another word.
   ;;  "\\(size_t\\)\\s-+\\(\\sw+\\)"
   ;;  '(0 font-lock-type-face prepend)
   ;;  '(2 font-lock-variable-name-face prepend))

   ;; keywords that take a variable name as an argument
   (list
    (concat "\\([@\\\\]\\(param\\(?:\\s-*\\[\\(?:in\\|out\\|in,out\\)\\]\\)?"
	    "\\|a\\|namespace\\|relates\\(also\\)?"
	    "\\|var\\|def\\)\\)\\s-+\\(\\sw+\\)")
    '(1 font-lock-keyword-face prepend)
    '(4 font-lock-variable-name-face prepend))
   ;; keywords that take a type name as an argument
   (list
    (concat "\\([@\\\\]\\(class\\|struct\\|union\\|exception\\|enum"
	    "\\|throw\\|interface\\|protocol\\)\\)\\s-+\\(\\(\\sw\\|:\\)+\\)")
    '(1 font-lock-keyword-face prepend)
    '(3 font-lock-type-face prepend))
   ;; keywords that take a function name as an argument
   (list
    "\\([@\\\\]retval\\)\\s-+\\([^ \t\n]+\\)"
    '(1 font-lock-keyword-face prepend)
    '(2 font-lock-function-name-face prepend))
   ;; bold
   (list
    "\\([@\\\\]b\\)\\s-+\\([^ \t\n]+\\)"
    '(1 font-lock-keyword-face prepend)
    '(2 (quote bold) prepend))
   ;; code
   (list
    "\\([@\\\\][cp]\\)\\s-+\\([^ \t\n]+\\)"
    '(1 font-lock-keyword-face prepend)
    '(2 (quote underline) prepend))
   ;; italics/emphasised
   (list
    "\\([@\\\\]e\\(m\\)?\\)\\s-+\\([^ \t\n]+\\)"
    '(1 font-lock-keyword-face prepend)
    '(3 (quote italic) prepend))
   ;; keywords that take a list
   (list
    "\\([@\\\\]ingroup\\)\\s-+\\(\\(\\sw+\\s-*\\)+\\)\\s-*$"
    '(1 font-lock-keyword-face prepend)
    '(2 font-lock-string-face prepend))
   ;; one argument that can contain arbitrary non-whitespace stuff
   (list
    (concat "\\([@\\\\]\\(link\\|copydoc\\|xrefitem"
	    "\\|if\\(not\\)?\\|elseif\\)\\)"
	    "\\s-+\\([^ \t\n]+\\)")
    '(1 font-lock-keyword-face prepend)
    '(4 font-lock-string-face prepend))
   ;; one optional argument that can contain arbitrary non-whitespace stuff
   (list
    "\\([@\\\\]\\(cond\\|dir\\)\\(\\s-+[^ \t\n]+\\)?\\)"
    '(1 font-lock-keyword-face prepend)
    '(3 font-lock-string-face prepend t))
   ;; one optional argument with no space between
   (list
    "\\([@\\\\]\\(~\\)\\([^ \t\n]+\\)?\\)"
    '(1 font-lock-keyword-face prepend)
    '(3 font-lock-string-face prepend t))
   ;; one argument that has to be a filename
   (list
    (concat "\\([@\\\\]\\(example\\|\\(dont\\)?include\\|includelineno"
	    "\\|htmlinclude\\|verbinclude\\)\\)\\s-+"
	    "\\(\"?[~:\\/a-zA-Z0-9_. ]+\"?\\)")
    '(1 font-lock-keyword-face prepend)
    '(4 font-lock-string-face prepend))
   ;; dotfile <file> ["caption"]
   (list
    (concat "\\([@\\\\]dotfile\\)\\s-+"
	    "\\(\"?[~:\\/a-zA-Z0-9_. ]+\"?\\)\\(\\s-+\"[^\"]+\"\\)?")
    '(1 font-lock-keyword-face prepend)
    '(2 font-lock-string-face prepend)
    '(3 font-lock-string-face prepend t))
   ;; image <format> <file> ["caption"] [<sizeindication>=<size>]
   (list
    "\\([@\\\\]image\\)\\s-+\\(html\\|latex\\)\\s-+\\(\"?[~:\\/a-zA-Z0-9_. ]+\"?\\)\\(\\s-+\"[^\"]+\"\\)?\\(\\s-+\\sw+=[0-9]+\\sw+\\)?"
    '(1 font-lock-keyword-face prepend)
    '(2 font-lock-string-face prepend)
    '(3 font-lock-string-face prepend)
    '(4 font-lock-string-face prepend t)
    '(5 font-lock-string-face prepend t))
   ;; one argument that has to be a word
   (list
    (concat "\\([@\\\\]\\(addtogroup\\|defgroup\\|weakgroup"
	    "\\|page\\|anchor\\|ref\\|section\\|subsection"
	    "\\)\\)\\s-+\\(\\sw+\\)")
    '(1 font-lock-keyword-face prepend)
    '(3 font-lock-string-face prepend))
))

(defun draco-font-lock ()
  "Turn on font-lock for Draco keywords."
  ;; FIXME How do I turn *off* font-lock for Doxygen keywords?
  (interactive)
  (if (functionp 'font-lock-add-keywords)
      ;; Use new (proper?) font-lock-add-keywords function
      (font-lock-add-keywords nil draco-mode-keywords)
    ;; Use old-school way
    (let ((old (if (eq (car-safe font-lock-keywords) t)
		 (cdr font-lock-keywords)
	       font-lock-keywords)))
      (setq font-lock-keywords (append old draco-mode-keywords)))))

;; http://www.emacswiki.org/emacs/AddKeywords
;(unless (fboundp 'font-lock-add-keywords)
;  (defalias 'font-lock-add-keywords 'ignore))

;(font-lock-add-keywords 'c++-mode
;  '(("\\<\\(size_t\\):" 1 font-lock-keyword-face prepend)
;    ("\\<\\(and\\|or\\|not\\)\\>" . font-lock-keyword-face)))

;;============================================================
;; Default templates
;;============================================================

(defun draco-mode-invalid-style ()
  "Warn the user that he has set `draco-mode-style' to an invalid
style."
  (error (concat
	  "Invalid `draco-mode-style': "
	  draco-mode-style
	  ": must be one of \"Qt\" or \"C++\".")))

;; FIXME
;; The following was borrowed from "simple.el".
;; If anyone knows of a better/simpler way of doing this, please let me know.
(defconst draco-mode-comment-indent-function
  (lambda (skip)
    (save-excursion
      (beginning-of-line)
      (let ((eol (save-excursion (end-of-line) (point))))
	(and skip
	     (re-search-forward skip eol t)
	     (setq eol (match-beginning 0)))
	(goto-char eol)
	(skip-chars-backward " \t")
	(max comment-column (1+ (current-column))))))
  "Function to compute desired indentation for a comment.
This function is called with skip and with point at the beginning of
the comment's starting delimiter.")

;;======================================================================
;; Helper Functions
;;======================================================================

(defun if-exist-load-library (libname)
  "If libname exists in load-path then load it.
If libname does not exist then print a warning message
and continue."
  (interactive)
    (if (load libname)
	(message (concat "Done loading file " libname ".el(c)"))
      (message (concat "Warning: " libname ".el not found"))))

(defun draco-set-companion-file-alist ()
"Setup a default list of companion files.  The user may use
his/her own bindings by setting the variable
draco-companion-file-alist from the Options:Advanced:Group:Draco menu."
  (if (not draco-companion-file-alist)
      (setq draco-companion-file-alist
	    '(
	      ;; Pair headers with implementation files.
	      ("\\(.h\\)$" . ".c")
	      ("\\(.H\\)$" . ".C")
	      ("\\(.hh\\)$" . ".cc")
	      ("\\(.hh\\)$" . ".i.hh")
	      ("\\(.hh\\)$" . ".t.hh")
	      ("\\(.hxx\\)$" . ".cxx")
	      ("\\(.hpp\\)$" . ".cpp")
	      ("\\(.h\\)$" . ".cpp")
	      ;; Pair implementation files with headers.
	      ("\\(.c\\)$" . ".h")
	      ("\\(.C\\)$" . ".H")
	      ("\\(.cc\\)$" . ".hh")
	      ("\\(.i.hh\\)$" . ".hh")
	      ("\\(.t.hh\\)$" . ".hh")
	      ("\\(.cxx\\)$" . ".hxx")
	      ("\\(.cpp\\)$" . ".hpp")
	      ("\\(.cpp\\)$" . ".h")
	      ))))
(add-hook 'draco-mode-hook 'draco-set-companion-file-alist)

;; Remove whitespace found at the end of every line before saving.
(add-hook 'before-save-hook 'delete-trailing-whitespace)

;;======================================================================
;; Comment blocks
;;======================================================================

(defun draco-insert-function-doc ()
"Default function for inserting a comment block in front of a C++
function or method. Width is set by draco-code-comment-width.

//---------------------------------------------------------------------------//
/*!
 * \brief
 *
 * \param name description
 * \return description
 */
"
  (interactive)
  (beginning-of-line)
  (insert "//")
  (insert (make-string (- draco-code-comment-width 4) ?-))
  (insert "//\n")
  (insert "/*! \n")
  (insert " * \\brief \n")
  (insert " * \n")
  (insert " * \\param[in] name description\n")
  (insert " * \\param[inout] name description\n")
  (insert " * \\param[out] name description\n")
  (insert " * \\return description\n")
  (insert " */\n")
  (previous-line 3)
  (end-of-line)
)
(defun draco-insert-class-doc ()
"Function for inserting a class desicription boilerplate.  Width is set
by draco-code-comment-width.

//===========================================================================//
/*!
 * \class
 * \brief
//===========================================================================//
"
  (interactive)
  (insert "//")
  (insert (make-string (- draco-code-comment-width 4) ?=))
  (insert "//\n")
  (insert "/*!\n")
  (insert " * \\class \n")
  (insert " * \\brief \n")
  (insert "//")
  (insert (make-string (- draco-code-comment-width 4) ?=))
  (insert "//\n")
  (previous-line 2)
  (end-of-line)
)

(defun draco-insert-comment-divider ()
"Function for inserting a single line comment divider."
  (interactive)
  (beginning-of-line)
  (insert "//")
  (insert (make-string (- draco-code-comment-width 4) ?-))
  (insert "//\n")
)

(defun draco-save-and-kill-current-buffer ()
"Save the current buffer and then kill it."
  (interactive)
  (if (buffer-file-name (current-buffer))
      (save-buffer))
  (kill-buffer (buffer-name)))

(defun draco-toggle-previous-buffer ()
"Toggle to previous buffer."
  (interactive)
  (switch-to-buffer (other-buffer (buffer-name))))

(defun draco-find-companion-file ()
  "
Function to locate the corresponding .hh .i.hh or .cc file.
When a .hh file is in the current buffer and this function is run,
the corresponding .cc file will be loaded if it is available.
If it is not available, the script will look for a corresponding
.i.hh file to load.

The mapping between file types is stored in the emacs variable
draco-companion-file-alist."

  (interactive)
  (let ((companion draco-companion-file-alist)
	(pair-file ""))
    (catch 'found
      (while companion
	(if (string-match (car (car companion)) buffer-file-name)
	    (progn
	      ;; Found a match, now check to see if its the right one.
	      (setq pair-file (replace-match (cdr (car companion))
					     t t buffer-file-name))
	      (if (file-exists-p pair-file)
		  (progn
		    (message "found matching file, throwing 'found")
		    (throw 'found t))
		(setq companion (cdr companion))))
	  (message (concat "discarding car companion=" pair-file))
	  (setq companion (cdr companion)))))
    (if companion
	(find-file pair-file))))

;;----------------------------------------------------------------------
;; Comment Dividers
;;----------------------------------------------------------------------

(defun draco-f90-subroutine-divider ()
"Insert a comment block for f90

!-----------------------------------------------------------------------------!
!
!-----------------------------------------------------------------------------!
"
  (interactive)
  (beginning-of-line)
  (insert "!")
  (insert (make-string (- draco-code-comment-width 3) ?-))
  (insert "!\n! \n!")
  (insert (make-string (- draco-code-comment-width 3) ?-))
  (insert "!\n\n")
  (previous-line 3)
  (end-of-line)
)

(defun draco-f90-comment-divider ()
"Insert a single line prepended by an f90 comment mark.

!-----------------------------------------------------------------------------!
"
  (interactive)
  (beginning-of-line)
  (insert "!")
  (insert (make-string (- draco-code-comment-width 3) ?-))
  (insert "!\n")
  (end-of-line)
)

(defun draco-f77-subroutine-divider ()
"Insert a comment block.

c-----------------------------------------------------------------------------c
c
c-----------------------------------------------------------------------------c
"
  (interactive)
  (beginning-of-line)
  (insert "c")
  (insert (make-string (- draco-code-comment-width 3) ?-))
  (insert "c\nc \nc")
  (insert (make-string (- draco-code-comment-width 3) ?-))
  (insert "c\n\n")
  (previous-line 3)
  (end-of-line)
)

(defun draco-f77-comment-divider ()
"Insert a divider.

c-----------------------------------------------------------------------------c
"
  (interactive)
  (beginning-of-line)
  (insert "c")
  (insert (make-string (- draco-code-comment-width 3) ?-))
  (insert "c\n")
  (end-of-line)
)


(defun draco-latex-divider ()
"Insert a LaTeX comment block."
  (interactive)
  (beginning-of-line)
  (insert "%%---------------------------------------------------------------------------%%\n")
  (insert "%% \n")
  (insert "%%---------------------------------------------------------------------------%%\n\n")
  (previous-line 3)
  (end-of-line)
)

(defun draco-latex-comment-divider ()
"Insert a LaTeX divider."
  (interactive)
  (beginning-of-line)
  (insert "%%---------------------------------------------------------------------------%%\n")
  (end-of-line)
)

(defun draco-makefile-divider ()
"Insert a Makefile comment block."
  (interactive)
  (beginning-of-line)
  (insert "#")
  (insert (make-string (- draco-code-comment-width 2) ?-))
  (insert "#\n")
  (insert "# \n")
  (insert "#")
  (insert (make-string (- draco-code-comment-width 2) ?-))
  (insert "#\n")
;  (insert "#-------------------------------------------------------------------------------#\n")
;  (insert "# \n")
;  (insert "#-------------------------------------------------------------------------------#\n\n")
  (previous-line 2)
  (end-of-line)
)

(defun draco-makefile-comment-divider ()
"Insert a Makefile divider."
  (interactive)
  (beginning-of-line)

  (insert "#")
  (insert (make-string (- draco-code-comment-width 2) ?-))
  (insert "#\n")
;  (insert "#--------------------------------------------------------------------------------#\n")
  (end-of-line)
)

(defun draco-m4-divider ()
"Insert a M4 comment block."
  (interactive)
  (beginning-of-line)
  (insert "dnl ------------------------------------------------------------------------- dnl\n")
  (insert "dnl \n")
  (insert "dnl ------------------------------------------------------------------------- dnl\n\n")
  (previous-line 3)
  (end-of-line)
)

(defun draco-m4-comment-divider ()
"Insert a M4 divider."
  (interactive)
  (beginning-of-line)
  (insert
   "dnl ------------------------------------------------------------------------- dnl\n")
  (end-of-line)
)

(defun draco-c-comment-divider ()
"Insert a C style divider.

/*---------------------------------------------------------------------------*/
"
  (interactive)
  (beginning-of-line)
  (insert "/*")
  (insert (make-string (- draco-code-comment-width 4) ?-))
  (insert "*/\n")
  (end-of-line)
)

(defun draco-c-divider ()
"Insert a C style comment block.

/*---------------------------------------------------------------------------*/
/* text */
/*---------------------------------------------------------------------------*/
"
  (interactive)
  (beginning-of-line)
  (insert "/*")
  (insert (make-string (- draco-code-comment-width 4) ?-))
  (insert "*/\n/* \n/*")
  (insert (make-string (- draco-code-comment-width 4) ?-))
  (insert "*/\n\n")
  (previous-line 3)
  (end-of-line)
)

(defun draco-html-comment-divider ()
"Insert a HTML style divider."
  (interactive)
  (beginning-of-line)
  (insert "<!---------------------------------------------------------------------------->\n")
  (end-of-line)
)

(defun draco-cc-divider ()
"Insert a C style comment block:
//---------------------------------------------------------------------------//
//
//---------------------------------------------------------------------------//"
  (interactive)
  (beginning-of-line)
  (insert "//")
  (insert (make-string (- draco-code-comment-width 4) ?-))
  (insert "//\n// \n//")
  (insert (make-string (- draco-code-comment-width 4) ?-))
  (insert "//\n\n")
  (previous-line 3)
  (end-of-line)
)

(defun draco-elisp-comment-divider ()
"Insert an elisp divider."
  (interactive)
  (beginning-of-line)
  (insert ";;---------------------------------------------------------------------------;;\n")
  (end-of-line)
)

(defun draco-elisp-divider ()
"Insert an elisp comment block."
  (interactive)
  (beginning-of-line)
  (insert ";;---------------------------------------------------------------------------;;\n")
  (insert ";; \n")
  (insert ";;---------------------------------------------------------------------------;;\n\n")
  (previous-line 3)
  (end-of-line)
)


;;---------------------------------------------------------------------------;;
;; time stamp function
;;---------------------------------------------------------------------------;;
(defvar draco-user-name nil
  "Name of current user")
(defun draco-name-and-time ()
  "Insert a time stamp"
  (interactive)
  (require 'time-stamp)
  (time-stamp-toggle-active 1)

  (setq draco-user-name (concat fill-prefix
				"Author: "
				(user-full-name)
				", " (user-login-name)
				"@lanl.gov\n" ))
  (defvar time-stamp-format "%a, %y %b %d, %H:%M:%S %Z")
  (defvar draco-time-stamp (concat fill-prefix
				 "Date  : "
				 (time-stamp-string)))
  (insert draco-user-name)
  (insert draco-time-stamp))

;;============================================================
;; Minor mode implementation
;;============================================================

(defvar draco-mode nil
  "nil disables draco-mode, non-nil enables.")

(make-variable-buffer-local 'draco-mode)
(require 'draco-menu)

(defun draco-mode (&optional arg)
  ;; All of the following text shows up in the "mode help" (C-h m)
  "Minor mode for using/creating Draco C++ source code.

To see what version of doxymacs you are running, enter
`\\[draco-mode-version]'.

Key bindings:
\\{draco-mode-map}"
  (interactive "P")
  (setq draco-mode
        (if (null arg)
            ;; Toggle mode
            (not draco-mode)
          ;; Enable/Disable according to arg
          (> (prefix-numeric-value arg) 0)))
  (if draco-mode
      (progn
	(run-hooks 'draco-mode-hook)))) ;; gets the draco-menu

;;============================================================

;; Provide this symbol

(provide 'draco-mode)

;;; draco-mode.el ends here
