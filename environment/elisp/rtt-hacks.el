;;; rtt-hacks.el
;;; Geoffrey Furnish
;;; 28 June 1996
;;;
;;; Some lispery to help (X)Emacs users easily create and navigate
;;; files which conform to the coding practices of any of a variety of
;;; projects they may be working on.  Each project may support
;;; different coding styles and file nomenclatures.  Such differences
;;; will hopefully be autosensed by the provided functions.
;;;-------------------------------------------------------------------------;;;
;;; $Id$
;;;-------------------------------------------------------------------------;;;

;;; These alists provide a mapping between places you might be in the
;;; filesystem, and places you should look for template header and
;;; implementation files.  "hh" means "header" and "cc" means
;;; "implementation".  These terms are drawn from C++ parlance, but it
;;; is entirely conceivable that these abstractions will be applicable
;;; to other languages as well.
;;; 
;;; The car of each cons cell represents a regexp which will be
;;; applied to a directory name for detecting what project you are in,
;;; and the cdr of the cons cell represents the location of the
;;; corresponding template file, plugged into the same directory name,
;;; replacing the text matched by the car regexp.

(defvar rtt-hh-template-alist
  '(
    ;; POOMA support.
    ("/r1/src.*$" . "/r1/src/config/template.h")
    ("/r1/test.*$" . "/r1/src/config/template.h")
    ;; Generic Furnishian expectations.
    ("/src.*$" . "/src/config/template.hh")
    )
  "Places to look for the template header file.")

(defvar rtt-cc-template-alist
  '(
    ;; POOMA support.
    ("/r1/src.*$" . "/r1/src/config/template.cpp")
    ("/r1/test.*$" . "/r1/src/config/template.cpp")
    ;; Generic Furnishian expectations.
    ("/src.*$" . "/src/config/template.cc")
    )
  "Places to look for the template implimentation file.")

;;; These next few functions use the above alists to determine which
;;; template files should be used when creating new files.

(defun rtt-find-template-file (template-alist)
  "Find the template file in the specified alist, based on the current
   project." 

  (let ((dir (expand-file-name "."))
	(template template-alist)
	(pair-file "")
	(result ""))
    (catch 'found
      (while template
	(if (string-match (car (car template)) dir)
	    (progn
	      ;; Found a match, now check to see if its the right one.
	      (setq pair-file (replace-match (cdr (car template))
					     t t dir))
	      (if (file-exists-p pair-file)
		  (progn
		    (message "found matching file, throwing 'found")
		    (throw 'found t))
		(setq template (cdr template))))
	  (message "discarding car template")
	  (setq template (cdr template)))))
    (if template
	(setq result pair-file)
      (error "Template file not found."))))

(defun rtt-find-hh-template ()
  "Find the header file template for the current project."

  (rtt-find-template-file rtt-hh-template-alist))

(defun rtt-find-cc-template ()
  "Find the implementation file template for the current project."

  (rtt-find-template-file rtt-cc-template-alist))

;;; These next few functions are used for figuring out where you are in
;;; the source tree.  These may need additional generalization.

(defun rtt-get-pkgdir ()
  "Find which package directory we are in (subdir of src),
   prompting user to change if desired."

  (setq dir (expand-file-name "."))
  (setq parent (expand-file-name ".."))

  ;; Figure out the difference in length between dir and parent,
  ;; subtract one (for the "/"), and negate, in order to get that many
  ;; chars off the end of dir.

  (setq dlen (- (length dir) (length parent)))
  (setq xlen (* (- dlen 1) -1))

  (read-from-minibuffer "Package: " (substring dir xlen))
  )

(defun rtt-get-src-dir ()
  "Find path to src dir, if possible based on the tree we are in."

  (setq dir (expand-file-name "."))
  (setq x (string-match "/src" dir))

  (if (null x)
      (progn
	(setq dir (read-file-name "Project source directory: "))
	(setq x (string-match "/src" dir))
	(if (null x)
	      (error (concat "Can't locate .../src in path "
			     dir))
	  )
	)
    )

  (setq srcdir (substring dir 0 (+ x 4))))

(defun rtt-file-suffix (filename)
  "Figure out the suffix for the specified filename."

  (string-match "\\.[A-Za-z]*$" filename)
  (setq result (match-string 0 filename)))

;;; Now come functions to actually construct new files, using the
;;; appropriate template files.  Note that currently only
;;; rtt-new-class is sufficiently general.  The others should be
;;; upgraded. 

(defun rtt-new-cc-file (name)
  "Function to spontaneously setup a new C++ .cc file."

  (interactive "sBase Name: ")

  (setq prjsrc (rtt-get-src-dir))

  (setq nfile (concat name ".cc"))
  (setq tfile (concat prjsrc "/config/template.cc"))
  (if (file-exists-p nfile)
      (error "File %s already exists." nfile))

  (find-file nfile)
  (insert-file-contents tfile)

  ;; Now perform customizations

  (find-file nfile)
  (perform-replace "<date>" (current-time-string) nil nil nil )
  (goto-char (point-min))
  (perform-replace "<user>" (user-full-name) nil nil nil )
  (goto-char (point-min))
  (perform-replace "<basename>" name nil nil nil )
  (goto-char (point-min))

  ;; Now, hack out the <start> token, since we don't know what to do
  ;; for a genearal file.

  (perform-replace "<start>" "" nil nil nil ))

(defun rtt-new-hh-file (name)
  "Function to spontaneously setup a new C++ .hh file."

  (interactive "sClass name: ")

  ;; Figure out which subdir of src we are in.

  (setq pkg (rtt-get-pkgdir))

  ;; Make a "safe" version of the package name.  This is important
  ;; because we need to be able to stuff something into CPP macros for
  ;; use in include guards, so it better not have any funky characters
  ;; in it.  A particularly common problem is plus signs, which are
  ;; cute in package names, but bad news in identifiers.

  (setq spkg (mapconcat (function (lambda (x)
				    (cond
				     ((eq x ?\+) "")
				     ((eq x ?/) "_")
				     (t (format "%c" x)))))
			pkg "" ))

  (setq prjsrc (rtt-get-src-dir))
  (setq pkgpath (concat prjsrc "/" pkg "/"))
  (setq hhfile (expand-file-name (concat pkgpath name ".hh")))

  (setq tfile (concat prjsrc "/config/template.hh"))
  (if (file-exists-p hhfile)
      (error "File %s already exists." hhfile))

  (find-file hhfile)
  (insert-file-contents tfile)

  ;; Now perform customizations

  (perform-replace "<pkg>" pkg nil nil nil )
  (goto-char (point-min))
  (perform-replace "<spkg>" spkg nil nil nil )
  (goto-char (point-min))
  (perform-replace "<class>" name nil nil nil )
  (goto-char (point-min))
  (perform-replace "<date>" (current-time-string) nil nil nil )
  (goto-char (point-min))
  (perform-replace "<user>" (user-full-name) nil nil nil ))

(defun rtt-new-class (name)
  "A function to produce the boilerplate for a new C++ class.
  Template files are copied into the (user specified) directory, and a
  variety of boilerplate text tokens are filled in with appropriate
  values."

  (interactive "sClass name: ")

  ;; Figure out which subdir of src we are in.

  (setq pkg (rtt-get-pkgdir))

  ;; Make a "safe" version of the package name.  This is important
  ;; because we need to be able to stuff something into CPP macros for
  ;; use in include guards, so it better not have any funky characters
  ;; in it.  A particularly common problem is plus signs, which are
  ;; cute in package names, but bad news in identifiers.

  (setq spkg (mapconcat (function (lambda (x)
				    (cond
				     ((eq x ?\+) "")
				     ((eq x ?/) "_")
				     (t (format "%c" x)))))
			pkg "" ))

  ;; Now locate the template files.

  (setq thhfile (rtt-find-hh-template))
  (setq tccfile (rtt-find-cc-template))
  (if (not (and (file-exists-p thhfile)
		(file-exists-p tccfile)))
      (error "Cannot access class template files."))

  ;; Now figure out the suffixes of the template files.

  (setq hhsfx (rtt-file-suffix thhfile))
  (setq ccsfx (rtt-file-suffix tccfile))

  ;; Make sure we can create the new files.

  (setq prjsrc (rtt-get-src-dir))
  (setq pkgpath (concat prjsrc "/" pkg "/"))
  (setq hhfile (expand-file-name (concat pkgpath name hhsfx)))
  (setq ccfile (expand-file-name (concat pkgpath name ccsfx)))
  (if (or (file-exists-p hhfile)
	  (file-exists-p ccfile))
      (error "Cannot create class %s, files already exist." name))

  ;; Now try to copy the template files over.

  (copy-file thhfile hhfile)
  (copy-file tccfile ccfile)

  ;; set permissions to make the files write-able
  (shell-command (concat "/bin/chmod ug+w " hhfile))
  (shell-command (concat "/bin/chmod ug+w " ccfile))

  ;; Now load the header
  (find-file hhfile)
  
  ;; customize header
  (perform-replace "<pkg>" pkg nil nil nil )
  (goto-char (point-min))
  (perform-replace "<spkg>" spkg nil nil nil )
  (goto-char (point-min))
  (perform-replace "<class>" name nil nil nil )
  (goto-char (point-min))
  (perform-replace "<date>" (current-time-string) nil nil nil )
  (goto-char (point-min))
  (perform-replace "<user>" (user-full-name) nil nil nil )

  ;; Now load the cc file and customize.
  (find-file ccfile)

  ;; customize header
  (perform-replace "<date>" (current-time-string) nil nil nil )
  (goto-char (point-min))
  (perform-replace "<user>" (user-full-name) nil nil nil )
  (goto-char (point-min))
  (perform-replace "<basename>" name nil nil nil )
  (goto-char (point-min))

  ;; Now, change the <start> token to include the header.

  (perform-replace "<start>" "" nil nil nil )
  (insert "#include \"" pkg "/" name ".hh\"\n\n"))

;;; An alist for registering associations between file pairs.

(defvar rtt-companion-alist
  '(
    ;; Pair headers with implementation files.
    ("\\(.h\\)$" . ".c")
    ("\\(.H\\)$" . ".C")
    ("\\(.hh\\)$" . ".cc")
    ("\\(.hh\\)$" . ".i.hh")
    ("\\(.hxx\\)$" . ".cxx")
    ("\\(.hpp\\)$" . ".cpp")
    ("\\(.h\\)$" . ".cpp")
    ;; Pair implementation files with headers.
    ("\\(.c\\)$" . ".h")
    ("\\(.C\\)$" . ".H")
    ("\\(.cc\\)$" . ".hh")
    ("\\(.i.hh\\)$" . ".hh")
    ("\\(.cxx\\)$" . ".hxx")
    ("\\(.cpp\\)$" . ".hpp")
    ("\\(.cpp\\)$" . ".h")
    )
  "
Alist of possible file pairs.  The file type of the current buffer
is located within the list (first regexp).  The filename extension 
is loaded from the 2nd part of the pair and a filename with this 
new extension is loaded (if it exists).")

;; Locate companion files by searching through an alist, looking for a
;; match for the file name of the current buffer.  Once found,
;; construct the matching file name.  If that file exists, load it.
;; Otherwise, assume the match must be a "nonstandard" file name pair,
;; so continue looking through the alist for the next match on the
;; current buffer's filenmae.  Iterate until success or failure.

(defun rtt-find-companion-file ()
  "
Function to locate the corresponding .hh .i.hh or .cc file.
When a .hh file is in the current buffer and this function is run, 
the corresponding .cc file will be loaded if it is available.
If it is not available, the script will look for a corresponding 
.i.hh file to load. 

The mapping between file types is stored in the emacs variable
rtt-companion-alist."

  (interactive)

  (let ((companion rtt-companion-alist)
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

;;; Default function for inserting a comment block in front of a C++
;;; function or method.

(defun rtt-insert-function-doc ()
  (interactive)
  (beginning-of-line)
  (insert "//---------------------------------------------------------------------------//\n")
  (insert "/*! \n")
  (insert " * \\brief \n")
  (insert " * \n")
  (insert " * \\param name description\n")
  (insert " * \\return description\n")
  (insert " */\n")
  (previous-line 3)
  (end-of-line)
)

;;; Same as above, but for C functions.

(defun rtt-insert-C-function-doc ()
  (interactive)
  (insert "/*--------------------------------------------------------------------------*\\\n")
  (insert " * \n")
  (insert " *\n")
  (insert " * \n")
  (insert "\\*--------------------------------------------------------------------------*/\n")
  (previous-line 4)
  (end-of-line)
)

;;; Function for inserting a class desicription boilerplate.

(defun rtt-insert-class-doc ()
  (interactive)
  (insert "//===========================================================================//\n")
  (insert "/*!\n")
  (insert " * \\class \n")
  (insert " * \\brief \n")
  (insert "//===========================================================================//\n")
  (previous-line 2)
  (end-of-line)
)

;;; Function for inserting a single line comment divider.

(defun rtt-insert-comment-divider ()
  (interactive)
  (beginning-of-line)
  (insert "//---------------------------------------------------------------------------//\n")
)

;;; A standard style for Radtran.  Note, this should not sound
;;; overbearing.  This style is actually extremely lax, very much like
;;; the XEmacs universal style.  Only difference is we ask for 4
;;; character basic block indenting, and a couple of really minor
;;; formatting constraints.  This is not intended to be a super strict
;;; fascist regime.

(defconst rtt-c-style
  `(
    (c-basic-offset	. 4)
    (c-offsets-alist	. ((access-label      . -2 )
			   (comment-intro     . 0)
			   (inline-open       . 0)
			   (innamespace       . 0)
			   (substatement-open . 0)))))

;;; The following is taken from the cc-mode info document, and shows
;;; some things we might wish to investigate at some point.

;(defconst rtt-c-style
;       '((c-tab-always-indent        . t) 
;         (c-comment-only-line-offset . 4) 
;	 (c-basic-offset	. 4)
;         (c-hanging-braces-alist     . ((substatement-open after) 
;                                        (brace-list-open))) 
;         (c-hanging-colons-alist     . ((member-init-intro before) 
;                                        (inher-intro) 
;                                        (case-label after) 
;                                        (label after) 
;                                        (access-label after))) 
;         (c-cleanup-list             . (scope-operator 
;                                        empty-defun-braces 
;                                        defun-close-semi)) 
;         (c-offsets-alist            . ((arglist-close . c-lineup-arglist) 
;                                        (substatement-open . 0) 
;                                        (case-label        . 4) 
;                                        (block-open        . 0) 
;                                        (knr-argdecl-intro . -))) 
;         (c-echo-syntactic-information-p . t) 
;         ) 

;;
;; Make this style known to XEmacs.
;;

(require 'cc-mode)

(defun rtt-add-style ()
  (message "adding rtt style")
  (c-add-style "rtt" rtt-c-style t)
)

(require 'infer-cc-style)

;;; Install our style, and configure some convenience settings.

(defun rtt-c-mode-common-hook ()
  (rtt-add-style)
  (c-set-style "rtt")
  (infer-cc-style)
  (define-key c-mode-map "\C-m" 'newline-and-indent)
  (set-fill-column 77)
  )

;;; Make Python code have similar indenting as our C++ code.  Also,
;;; force Emacs to use spaces instead of tabs, since mixed mode causes
;;; problems with some editors.

(defun rtt-python-mode-hook ()
  (setq py-indent-offset 4)
  (setq indent-tabs-mode nil)
  )

;;; Configure the Radtran coding environment.  Mostly this means
;;; configuring editing modes for C++ and Python, but anything else we
;;; need will go here too.  Note in particular that we do not attempt
;;; to do any serious keystroke binding here, since people have so
;;; many different types of machines.

(defun rtt-setup ()	
  (interactive)
  (message "Configuring RadTran support.")

  ;; Set up C++ mode style variables for radtran team.
  (add-hook 'c-mode-common-hook 'rtt-c-mode-common-hook)
  (add-hook 'python-mode-hook   'rtt-python-mode-hook)

  ;; Enable mouse scroll wheel.
  ;;(if (require 'mwheel) (mwheel-install))  -- This fails on Theta so
  ;;use a work around:
  (if (load "mwheel" t t nil) (mwheel-install) (message "Cannot find mwheel.elc"))

  )

(provide 'rtt-hacks)
