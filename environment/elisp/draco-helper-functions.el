;;---------------------------------------------------------------------------;;
;; draco-helper-functions.el
;; Kelly Thompson
;; 20 January 2005
;;---------------------------------------------------------------------------;;

;;---------------------------------------------------------------------------;;
;; Draco mode helper functions:
;;
;; draco-get-local-dir-name (&optional ndirs)
;; draco-guess-package-name ()
;; draco-make-name-safe (name)
;; draco-query-namespace (guess)
;; draco-create-file-from-template( filename, template, &optional class-name )
;; draco-guess-names
;;---------------------------------------------------------------------------;;

;;---------------------------------------------------------------------------;;
;; Define some variables.

(defvar draco-package-name nil
  "The name of the current Draco package.  This variable is used when
generating new files (new class, new header, etc.) in a package
directory.")

(defvar draco-safe-package-name nil
  "This value is normally the same as draco-package-name.  However, if
draco-package-name contains invalid characters, then
draco-safe-package-name contains a modified package-name that has
appropriate substitutions for invalid characters.  For now the
following two transformations take place:

+ -> \"\"
/ -> \"_\"")

(defvar draco-namespace nil
  "This value contains the package local namespace name.  This
variable is used during the creating of new files within a package
directory. ")

(defvar draco-paper-name nil
  "This value contains the name of the LaTeX document that is being
created.")

;;---------------------------------------------------------------------------;;
;; Get name of current directory

(defun draco-get-local-dir-name (&optional ndirs)
  "Function to obtain the local directory name.

Argument \"ndirs\" specifies the number leafs to present.

The default value of ndirs == 1."

  (if (not ndirs)
      (setq ndirs 1))

  (setq dir (expand-file-name "."))
  (setq relative "..")
  (setq index 1)
  (while (< index ndirs)
    (setq relative (concat relative "/.."))
    (setq index (1+ index)))

  (setq parent (expand-file-name relative))

  ;; Figure out the difference in length between dir and parent,
  ;; subtract one (for the "/"), and negate, in order to get that many
  ;; chars off the end of dir.

  (setq dlen (- (length dir) (length parent)))
  (setq xlen (* (- dlen 1) -1))

 (substring dir xlen)
)

;;---------------------------------------------------------------------------;;
;; Query for the pkg name, with guessed default.

(defun draco-guess-package-name ()
  "Function to guess a package directory name."
;  (setq dir (expand-file-name "."))
;  (setq parent (expand-file-name ".."))
  (setq draco-guess-query-text "Package path: ")
  (setq draco-package-name (draco-get-local-dir-name 1))

  (if (string= draco-package-name "test")
      (progn
	(setq draco-package-name (substring (draco-get-local-dir-name
					     2) 0 -5))
	(setq draco-guess-query-text "Package name: ")))

  (if (string= draco-package-name "autodoc")
      (progn
	(setq draco-package-name (draco-get-local-dir-name 2))
	(setq draco-guess-query-text "Package doc path: ")))

  (read-from-minibuffer draco-guess-query-text draco-package-name)
)

;;---------------------------------------------------------------------------;;
;; Return a "safe" name from a provided name.

(defun draco-make-name-safe (name)
  "Function to get a safe name from a package name"
  (setq safe-name
	(mapconcat (function (lambda (x)
			       (cond
				((eq x ?\+) "x")
				((eq x ?/) "_")
				(t (format "%c" x)))))
		   name "" ))
)

;;---------------------------------------------------------------------------;;
;; Query for the desired namespace

(defun draco-query-namespace (guess)
  "Function to query the user for the desired namespace given a default guess"
  (read-from-minibuffer "Namespace for this package: " guess)
)

;;---------------------------------------------------------------------------;;
;; Create a new buffer, load the specified template, replace special
;; text with appropriate values.
(defun draco-create-buffer-from-template (filename
					  template-filename
					  &optional draco-alternate-class-name)
"Load the file specified by template-filename into a new buffer,
replace special text with approprate values.  The user must save the
buffer to actually create a new file.

This function requires that the following variables be
set and available:

draco-package-name
draco-safe-package-name
draco-namespace
draco-class-name
draco-paper-name (for LaTeX documents)

ARGUMENTS:

 - FILENAME is the name of the file on disk that is being created.

 - TEMPLATE-FILENAME is the name of the template file located in the
   directory specified by draco-templates-dir.

 - &optional DRACO-ALTERNATE-CLASS-NAME If this variable is provided,
   this routine will use its value for the class name instead of
   guessing the name.

 - &optional DRACO-DOC-TYPE defaults to nil (a document containing a
   program).  It can optionally be set to \"latex\" for the generation
   of new LaTeX documents from template files.

This function is typically called from draco-package,
draco-package-test, draco-package-doc, draco-class, draco-cc-head,
draco-cc-headin, draco-c-head, draco-c-headin, draco-cc-imp,
draco-cc-pt, and draco-cc-test.

Example usage is exhibited by functions defined in draco-new-files.el.
"
  ;; Require the file "filename" to be new (if it already exists, then
  ;; abort.
  (if (file-exists-p filename)
      (error "File %s already exists. Aborting new file command."
	     filename))

    ;; obtain draco-package-name, draco-safe-package-name, and draco-namespace
;    (draco-guess-names)

    ;; If the optional draco-alternate-class-name is provided use it.
    (if draco-alternate-class-name
	(setq draco-class-name draco-alternate-class-name))

    ;; If draco-class-name is not bound, then bind it to the filename
    ;; minus any extensions.
    (if (not (boundp 'draco-class-name))
	(setq draco-class-name filename))
;      (setq draco-class-name (substring filename 0 (string-match "\\." filename))))

    ;; Create a new buffer with the correct filename - user must save
    ;; the buffer for a new file to be created.
    (find-file filename)

    ;; Insert the text from the template file into the new buffer.
    (insert-file-contents template-filename)

    (perform-replace "<pkg>"       draco-package-name nil nil nil )
    (goto-char (point-min))
    (perform-replace "<spkg>"      draco-safe-package-name nil nil nil )
    (goto-char (point-min))
    (perform-replace "<tpkg>"      (concat draco-package-name "/test") nil nil nil )
    (goto-char (point-min))
    (perform-replace "<namespace>" draco-namespace nil nil nil )
    (goto-char (point-min))
    (perform-replace "<user>"      (user-full-name) nil nil nil )
    (goto-char (point-min))
    (perform-replace "<class>"     draco-class-name nil nil nil )
    (goto-char (point-min))
    (perform-replace "<basename>"  draco-class-name nil nil nil )
    (goto-char (point-min))
    (perform-replace "<date>"      (current-time-string) nil nil nil )
    (goto-char (point-min))
    (perform-replace "<papername>" draco-paper-name nil nil nil)
    (goto-char (point-min))
    (perform-replace "<start>" ""  nil nil nil )
    )

;;---------------------------------------------------------------------------;;
;; draco-guess-names

(defun draco-guess-names ( &optional draco-doc-type )
  "Function to guess package and namespace name.  Also, removes invalid
characters from these names. Defines the variables:

draco-package-name
draco-safe-package-name

draco-namespace  (if draco-doc-type is nil)
draco-paper-name (if draco-doc-type is non-nil)
                 must be either \"latex\" or \"bibtex\"

ARGUMENTS:

 - &optional DRACO-DOC-TYPE defaults to nil (a document containing a
   program).  It can optionally be set to \"latex\" for the generation
   of new LaTeX documents from template files.
"
  ;; set values.

  (setq draco-package-name
	(draco-guess-package-name))
  (setq draco-safe-package-name
	(draco-make-name-safe draco-package-name))

  (if draco-doc-type
      (if (string= draco-doc-type "latex")
	  ;; strip .tex from the provided memo name (if it exists)
	  (setq draco-paper-name (substring name 0 (string-match "\\.tex" name)))
	;; strip .tex from the provided memo name (if it exists)
	(setq draco-paper-name (substring name 0 (string-match "\\.bib" name))))

    ;; We get here if draco-doc-type == nil
    (setq draco-namespace
	  (draco-query-namespace (concat "rtt_" draco-safe-package-name))))

  ;; default draco-class-name to an empty string.

  )

;;---------------------------------------------------------------------------;;

(provide 'draco-helper-functions)
