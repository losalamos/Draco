;;; tme-db.el --- TME-db mode for GNU Emacs

;; Copyright (C) 1992, 1994, 1995, 1996, 1997 Free Software Foundation, Inc.

;; Author: Stefan Schoef <schoef@offis.uni-oldenburg.de>
;;	Bengt Martensson <bengt@mathematik.uni-Bremen.de>
;;	Mark Shapiro <shapiro@corto.inria.fr>
;;	Mike Newton <newton@gumby.cs.caltech.edu>
;;	Aaron Larson <alarson@src.honeywell.com>
;; Maintainer: Dirk Herrmann <D.Herrmann@tu-bs.de>
;; Keywords: TME-db, LaTeX, TeX

;; This file is part of GNU Emacs.

;; GNU Emacs is free software; you can redistribute it and/or modify
;; it under the terms of the GNU General Public License as published by
;; the Free Software Foundation; either version 2, or (at your option)
;; any later version.

;; GNU Emacs is distributed in the hope that it will be useful,
;; but WITHOUT ANY WARRANTY; without even the implied warranty of
;; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;; GNU General Public License for more details.

;; You should have received a copy of the GNU General Public License
;; along with GNU Emacs; see the file COPYING.  If not, write to the
;; Free Software Foundation, Inc., 59 Temple Place - Suite 330,
;; Boston, MA 02111-1307, USA.

;;; Commentary:

;;  Major mode for editing and validating TME-db files.

;;  Usage:
;;  See documentation for function tme-db-mode (or type "\M-x describe-mode"
;;  when you are in TME-db mode).

;;  Todo:
;;  Distribute texinfo file.

;;; Code:

(eval-when-compile
  (require 'compile))

;; User Options:

(defgroup tme-db nil
  "TME-db mode."
  :group 'tex
  :prefix "tme-db-")

(defgroup tme-db-autokey nil
  "Generates automatically a key from the author/editor and the title field"
  :group 'tme-db
  :prefix "tme-db-autokey-")

(defcustom tme-db-mode-hook nil
  "List of functions to call on entry to TME-db mode."
  :group 'tme-db
  :type 'hook)

(defcustom tme-db-field-delimiters 'braces
  "*Controls type of field delimiters used.
Set this to `braces' or `double-quotes' according to your personal
preferences.  This variable is buffer-local."
  :group 'tme-db
  :type '(choice (const braces)
		 (const double-quotes)))
(make-variable-buffer-local 'tme-db-field-delimiters)

(defcustom tme-db-entry-delimiters 'braces
  "*Controls type of entry delimiters used.
Set this to `braces' or `parentheses' according to your personal
preferences.  This variable is buffer-local."
  :group 'tme-db
  :type '(choice (const braces)
		 (const parentheses)))
(make-variable-buffer-local 'tme-db-entry-delimiters)

(defcustom tme-db-include-OPTcrossref '("InProceedings" "InCollection")
  "*All entries listed here will have an OPTcrossref field."
  :group 'tme-db
  :type '(repeat string))

(defcustom tme-db-include-OPTkey nil
  "*If non-nil, all entries will have an OPTkey field.
If this is a string, it will be used as the initial field text.
If this is a function, it will be called to generate the initial field text."
  :group 'tme-db
  :type '(choice (const :tag "None" nil)
		 (const :tag "Default" t)
		 (string :tag "Initial text")
		 (function :tag "Initialize Function" :value fun)))

(defcustom tme-db-user-optional-fields
  '(("annote" "Personal annotation (ignored)"))
  "*List of optional fields the user wants to have always present.
Entries should be of the same form as the OPTIONAL and
CROSSREF-OPTIONAL lists in `tme-db-entry-field-alist' (see documentation
of this variable for details)."
  :group 'tme-db
  :type '(repeat
	  (group (string :tag "Field")
		 (string :tag "Comment")
		 (option (group :inline t
				:extra-offset -4
				(choice :tag "Init" :value ""
					string
					function))))))

(defcustom tme-db-entry-format '(opts-or-alts numerical-fields)
  "*Controls type of formatting performed by `tme-db-clean-entry'.
It may be t, nil, or a list of symbols out of the following: 
opts-or-alts        Delete empty optional and alternative fields and
                      remove OPT and ALT prefixes from used fields.
numerical-fields    Delete delimiters around numeral fields.
page-dashes         Change double dashes in page field to single dash
                      (for scribe compatibility).
inherit-booktitle   If entry contains a crossref field and booktitle
                      field is empty, it is set to the contents of the
                      title field of the crossreferenced entry.
                      Caution: this will work only if buffer is
                       correctly sorted.
realign             Realign entries, so that field texts and perhaps equal
                      signs (depending on the value of
                      tme-db-align-at-equal-sign) begin in the same column.
last-comma          Add or delete comma on end of last field in entry,
                      according to value of `tme-db-comma-after-last-field'.
delimiters          Change delimiters according to variables
                      `tme-db-field-delimiters' and `tme-db-entry-delimiters'.
unify-case          Change case of entry and field names.

The value t means do all of the above formatting actions.
The value nil means do no formatting at all."
  :group 'tme-db
  :type '(choice (const :tag "None" nil)
		 (const :tag "All" t)
		 (set :menu-tag "Some"
		      (const opts-or-alts)
		      (const numerical-fields)
		      (const page-dashes)
		      (const inherit-booktitle)
		      (const realign)
		      (const last-comma)
		      (const delimiters)
		      (const unify-case))))

(defcustom tme-db-clean-entry-hook nil
  "*List of functions to call when entry has been cleaned.
Functions are called with point inside the cleaned entry, and the buffer
narrowed to just the entry."
  :group 'tme-db
  :type 'hook)

(defcustom tme-db-sort-ignore-string-entries t
  "*If non-nil, TME-db @String entries are not sort-significant.
That means they are ignored when determining ordering of the buffer
(e.g. sorting, locating alphabetical position for new entries, etc.).
This variable is buffer-local."
  :group 'tme-db
  :type 'boolean)
(make-variable-buffer-local 'tme-db-sort-ignore-string-entries)

(defcustom tme-db-maintain-sorted-entries nil
  "*If non-nil, TME-db mode maintains all TME-db entries in sorted order.
Setting this variable to nil will strip off some comfort (e.g. TAB
completion for reference keys in minibuffer, automatic detection of
duplicates) from TME-db mode.  See also `tme-db-sort-ignore-string-entries'.
This variable is buffer-local."
  :group 'tme-db
  :type 'boolean)
(make-variable-buffer-local 'tme-db-maintain-sorted-entries)

(defcustom tme-db-field-kill-ring-max 20
  "*Max length of `tme-db-field-kill-ring' before discarding oldest elements."
  :group 'tme-db
  :type 'integer)

(defcustom tme-db-entry-kill-ring-max 20
  "*Max length of `tme-db-entry-kill-ring' before discarding oldest elements."
  :group 'tme-db
  :type 'integer)

(defcustom tme-db-parse-keys-timeout 60
  "*Specifies interval for parsing buffers.
All TME-db buffers in emacs are parsed if emacs has been idle
`tme-db-parse-keys-timeout' seconds.  Only buffers which were modified
after last parsing and which are maintained in sorted order are parsed."
  :group 'tme-db
  :type 'integer)

(defvar tme-db-entry-field-alist
  '(
    ("Project" . (((("members" "Team members on the project")
		    ("name" "Project name")
		    ("objective" "Objective of project")
		    ("status" "Project open or closed"))
		   (("subname" "Subname of project")
		    ("begin" "Beginning date of project")
		    ("end" "Ending date of project")
		    ("subjects" "Subject category of project")
		    ("tasks" "Major tasks associated with project")
		    ("tasks" "Major tasks associated with project")
		    ("tasks" "Major tasks associated with project")
		    ("tasks" "Major tasks associated with project")
		    ("mileposts" "Mileposts associated with project")
		    ("codes" "Associated computer codes")
		    ("associations" "Associated projects and/or people")
		    ("logs" "Log notebooks with data about the project")
		    ("archives" "Archival storage of data associated with the project")
		    ("CVS" "CVS modules associated with the project")
		    ("note" "Any extra notes"))
		   (()()))))
    )
  "Defines reference types and their associated fields.
List of
(ENTRY-NAME (REQUIRED OPTIONAL) (CROSSREF-REQUIRED CROSSREF-OPTIONAL))
triples.
If the third element is nil, the first pair is always used.
If not, the second pair is used in the case of presence of a crossref
field and the third in the case of absence.
REQUIRED, OPTIONAL, CROSSREF-REQUIRED and CROSSREF-OPTIONAL are lists. 
Each element of these lists is a list of the form
(FIELD-NAME COMMENT-STRING INIT ALTERNATIVE-FLAG).
COMMENT-STRING, INIT, and ALTERNATIVE-FLAG are optional.
FIELD-NAME is the name of the field, COMMENT-STRING the comment to
appear in the echo area, INIT is either the initial content of the
field or a function, which is called to determine the initial content
of the field, and ALTERNATIVE-FLAG (either nil or t) marks if the
field is an alternative.  ALTERNATIVE-FLAG may be t only in the
REQUIRED or CROSSREF-REQUIRED lists.")

(defvar tme-db-comment-start "@Comment ")

(defcustom tme-db-add-entry-hook nil
  "List of functions to call when entry has been inserted."
  :group 'tme-db
  :type 'hook)

(defcustom tme-db-predefined-month-strings
  '(
    ("JAN") ("FEB") ("MAR") ("APR") ("MAY") ("JUN")
    ("JUL") ("AUG") ("SEP") ("OCT") ("NOV") ("DEC")
    )
  "Alist of month string definitions.
Should contain all strings used for months in the TME-db style files.
Each element is a list with just one element: the string."
  :group 'tme-db
  :type '(repeat
	  (list string)))

(defcustom tme-db-predefined-strings
  (append
   tme-db-predefined-month-strings
   '(
     ("imc") ("mc") ("fem") ("transport") ("radhydro")
     ("draco") ("milagro") ("milstone") ("ddimc") ("ddmc")
     ("dm") ("attila") ("dante") ("avatar") ("mcnp")
     ))
  "Alist of string definitions.
Should contain the strings defined in the TME-db style files.  Each
element is a list with just one element: the string."
  :group 'tme-db
  :type '(repeat
	  (list string)))

(defcustom tme-db-string-files nil
  "*List of TME-db files containing string definitions.
Those files must be specified using pathnames relative to the
directories specified in `tme-db-string-file-path'.  This variable is only
evaluated when TME-db mode is entered (i. e. when loading the TME-db
file)."
  :group 'tme-db
  :type '(repeat file))

(defvar tme-db-string-file-path (getenv "BIBINPUTS")
  "*Colon separated list of pathes to search for `tme-db-string-files'.")

(defcustom tme-db-help-message t
  "*If not nil print help messages in the echo area on entering a new field."
  :group 'tme-db
  :type 'boolean)

(defcustom tme-db-autokey-prefix-string ""
  "*String to use as a prefix for all generated keys.
See the documentation of function `tme-db-generate-autokey' for more detail."
  :group 'tme-db-autokey
  :type 'string)

(defcustom tme-db-autokey-names 1
  "*Number of names to use for the automatically generated reference key.
If this is variable is nil, all names are used.
Possibly more names are used according to `tme-db-autokey-names-stretch'.
See the documentation of function `tme-db-generate-autokey' for more detail."
  :group 'tme-db-autokey
  :type '(choice (const :tag "All" infty)
		 integer))

(defcustom tme-db-autokey-names-stretch 0
  "*Number of names that can additionally be used.
These names are used only, if all names are used then.
See the documentation of function `tme-db-generate-autokey' for details."
  :group 'tme-db-autokey
  :type 'integer)

(defcustom tme-db-autokey-additional-names ""
  "*String to prepend to the generated key if not all names could be used.
See the documentation of function `tme-db-generate-autokey' for details."
  :group 'tme-db-autokey
  :type 'string)

(defvar tme-db-autokey-transcriptions
  '(
    ;; language specific characters
    ("\\\\aa" "a")
    ("\\\\AA" "A")
    ("\\\"a\\|\\\\\\\"a\\|\\\\ae" "ae")
    ("\\\"A\\|\\\\\\\"A\\|\\\\AE" "Ae")
    ("\\\\i" "i")
    ("\\\\j" "j")
    ("\\\\l" "l")
    ("\\\\L" "L")
    ("\\\"o\\|\\\\\\\"o\\|\\\\o\\|\\\\oe" "oe")
    ("\\\"O\\|\\\\\\\"O\\|\\\\O\\|\\\\OE" "Oe")
    ("\\\"s\\|\\\\\\\"s" "ss")
    ("\\\"u\\|\\\\\\\"u" "ue")
    ("\\\"U\\|\\\\\\\"U" "Ue")
    ;; accents
    ("\\\\`\\|\\\\'\\|\\\\\\^\\|\\\\~\\|\\\\=\\|\\\\\\.\\|\\\\u\\|\\\\v\\|\\\\H\\|\\\\t\\|\\\\c\\|\\\\d\\|\\\\b" "")
    ;; braces
    ("{" "") ("}" ""))
  "Alist of (old-regexp new-string) pairs.
Used by the default values of `tme-db-autokey-name-change-strings' and
`tme-db-autokey-titleword-change-strings'.  Defaults to translating some
language specific characters to their ASCII transcriptions, and
removing any character accents.")

(defcustom tme-db-autokey-name-change-strings
  tme-db-autokey-transcriptions
  "Alist of (OLD-REGEXP NEW-STRING) pairs.
Any part of name matching a OLD-REGEXP is replaced by NEW-STRING.
Case is significant in OLD-REGEXP.  All regexps are tried in the
order in which they appear in the list, so be sure to avoid inifinite
loops here.
See the documentation of function `tme-db-generate-autokey' for details."
  :group 'tme-db-autokey
  :type '(repeat
	  (list (regexp :tag "Old")
		(string :tag "New"))))

(defcustom tme-db-autokey-name-case-convert 'downcase
  "*Function called for each name to perform case conversion.
See the documentation of function `tme-db-generate-autokey' for more detail."
  :group 'tme-db-autokey
  :type '(choice (const :tag "Preserve case" identity)
		 (const :tag "Downcase" downcase)
		 (const :tag "Capitalize" capitalize)
		 (const :tag "Upcase" upcase)
		 (function :tag "Conversion function")))

(defcustom tme-db-autokey-name-length 'infty
  "*Number of characters from name to incorporate into key.
If this is set to anything but a number, all characters are used.
See the documentation of function `tme-db-generate-autokey' for details."
  :group 'tme-db-autokey
  :type '(choice (const :tag "All" infty)
		 integer))

(defcustom tme-db-autokey-name-separator ""
  "*String that comes between any two names in the key.
See the documentation of function `tme-db-generate-autokey' for details."
  :group 'tme-db-autokey
  :type 'string)

(defcustom tme-db-autokey-year-length 2
  "*Number of rightmost digits from the year field yo incorporate into key.
See the documentation of function `tme-db-generate-autokey' for details."
  :group 'tme-db-autokey
  :type 'integer)

(defcustom tme-db-autokey-year-use-crossref-entry t
  "*If non-nil use year field from crossreferenced entry if necessary.
If this variable is non-nil and the current entry has no year, but a
valid crossref entry, the year field from the crossreferenced entry is
used.
See the documentation of function `tme-db-generate-autokey' for details."
  :group 'tme-db-autokey
  :type 'boolean)

(defcustom tme-db-autokey-titlewords 5
  "*Number of title words to use for the automatically generated reference key.
If this is set to anything but a number, all title words are used.
Possibly more words from the title are used according to
`tme-db-autokey-titlewords-stretch'.
See the documentation of function `tme-db-generate-autokey' for details."
  :group 'tme-db-autokey
  :type '(choice (const :tag "All" infty)
		 integer))

(defcustom tme-db-autokey-title-terminators
  '("\\." "!"  "\\?" ":" ";" "--")
  "*Regexp list defining the termination of the main part of the title.
Case of the regexps is ignored.
See the documentation of function `tme-db-generate-autokey' for details."
  :group 'tme-db-autokey
  :type '(repeat regexp))

(defcustom tme-db-autokey-titlewords-stretch 2
  "*Number of words that can additionally be used from the title.
These words are used only, if a sentence from the title can be ended then.
See the documentation of function `tme-db-generate-autokey' for details."
  :group 'tme-db-autokey
  :type 'integer)

(defcustom tme-db-autokey-titleword-ignore
  '("A" "An" "On" "The" "Eine?" "Der" "Die" "Das"
    "[^A-Z].*" ".*[^a-zA-Z0-9].*")
  "*Determines words from the title that are not to be used in the key.
Each item of the list is a regexp.  If a word of the title matchs a
regexp from that list, it is not included in the title part of the key.
See the documentation of function `tme-db-generate-autokey' for details."
  :group 'tme-db-autokey
  :type '(repeat regexp))

(defcustom tme-db-autokey-titleword-case-convert 'downcase
  "*Function called for each titleword to perform case conversion.
See the documentation of function `tme-db-generate-autokey' for more detail."
  :group 'tme-db-autokey
  :type '(choice (const :tag "Preserve case" identity)
		 (const	:tag "Downcase" downcase)
		 (const	:tag "Capitalize" capitalize)
		 (const	:tag "Upcase" upcase)
		 (function :tag "Conversion function")))

(defcustom tme-db-autokey-titleword-abbrevs nil
  "*Determines exceptions to the usual abbreviation mechanism.
An alist of (OLD-REGEXP NEW-STRING) pairs.  Case is ignored
in matching against OLD-REGEXP, and the first matching pair is used.
See the documentation of function `tme-db-generate-autokey' for details.")

(defcustom tme-db-autokey-titleword-change-strings
  tme-db-autokey-transcriptions
  "Alist of (OLD-REGEXP NEW-STRING) pairs.
Any part of title word matching a OLD-REGEXP is replaced by NEW-STRING.
Case is significant in OLD-REGEXP.  All regexps are tried in the
order in which they appear in the list, so be sure to avoid inifinite
loops here.
See the documentation of function `tme-db-generate-autokey' for details."
  :group 'tme-db-autokey
  :type '(repeat
	  (list (regexp :tag "Old")
		(string :tag "New"))))

(defcustom tme-db-autokey-titleword-length 5
  "*Number of characters from title words to incorporate into key.
If this is set to anything but a number, all characters are used.
See the documentation of function `tme-db-generate-autokey' for details."
  :group 'tme-db-autokey
  :type '(choice (const :tag "All" infty)
		 integer))

(defcustom tme-db-autokey-titleword-separator "_"
  "*String to be put between the title words.
See the documentation of function `tme-db-generate-autokey' for details."
  :group 'tme-db-autokey
  :type 'string)

(defcustom tme-db-autokey-name-year-separator ""
  "*String to be put between name part and year part of key.
See the documentation of function `tme-db-generate-autokey' for details."
  :group 'tme-db-autokey
  :type 'string)

(defcustom tme-db-autokey-year-title-separator ":_"
  "*String to be put between name part and year part of key.
See the documentation of function `tme-db-generate-autokey' for details."
  :group 'tme-db-autokey
  :type 'string)

(defcustom tme-db-autokey-edit-before-use t
  "*If non-nil, user is allowed to edit the generated key before it is used."
  :group 'tme-db-autokey
  :type 'boolean)

(defcustom tme-db-autokey-before-presentation-function nil
  "Function to call before the generated key is presented.
If non-nil this should be a single function, which is called before
the generated key is presented (in entry or, if
`tme-db-autokey-edit-before-use' is t, in minibuffer). This function
must take one argument (the automatically generated key), and must
return with a string (the key to use)."
  :group 'tme-db-autokey
  :type 'function)

(defcustom tme-db-entry-offset 0
  "*Offset for TME-db entries.
Added to the value of all other variables which determine colums."
  :group 'tme-db
  :type 'integer)

(defcustom tme-db-field-indentation 2
  "*Starting column for the name part in TME-db fields."
  :group 'tme-db
  :type 'integer)

(defcustom tme-db-text-indentation
  (+
   tme-db-field-indentation
   (length "organization = "))
  "*Starting column for the text part in TME-db fields.
Should be equal to the space needed for the longest name part."
  :group 'tme-db
  :type 'integer)

(defcustom tme-db-contline-indentation
  (+ tme-db-text-indentation 1)
  "*Starting column for continuation lines of TME-db fields."
  :group 'tme-db
  :type 'integer)

(defcustom tme-db-align-at-equal-sign nil
  "*If non-nil, align fields at equal sign instead of field text.
If non-nil, the column for the equal sign is
the value of `tme-db-text-indentation', minus 2."
  :group 'tme-db
  :type 'boolean)

(defcustom tme-db-comma-after-last-field nil
  "*If non-nil, a comma is put at end of last field in the entry template."
  :group 'tme-db
  :type 'boolean)

;; tme-db-font-lock-keywords is a user option as well, but since the
;; patterns used to define this variable are defined in a later
;; section of this file, it is defined later.

;; Special support taking care of variants
(if (boundp 'mark-active)
    (defun tme-db-mark-active ()
      ;; In Emacs mark-active indicates if mark is active.
      mark-active)
  (defun tme-db-mark-active ()
    ;; In XEmacs (mark) returns nil when not active.
    (if zmacs-regions (mark) (mark t))))

(if (fboundp 'run-with-idle-timer)
    ;; timer.el is distributed with Emacs
    (fset 'tme-db-run-with-idle-timer 'run-with-idle-timer)
  ;; timer.el is not distributed with XEmacs
  ;; Notice that this does not (yet) pass the arguments, but they
  ;; are not used (yet) in tme-db.el. Fix if needed.
  (defun tme-db-run-with-idle-timer (secs repeat function &rest args)
    (start-itimer "tme-db" function secs (if repeat secs nil) t)))


;; Syntax Table, Keybindings and TME-db Entry List
(defvar tme-db-mode-syntax-table
  (let ((st (make-syntax-table)))
    (modify-syntax-entry ?\" "\"" st)
    (modify-syntax-entry ?$ "$$  " st)
    (modify-syntax-entry ?% "<   " st)
    (modify-syntax-entry ?' "w   " st)
    (modify-syntax-entry ?@ "w   " st)
    (modify-syntax-entry ?\\ "\\" st)
    (modify-syntax-entry ?\f ">   " st)
    (modify-syntax-entry ?\n ">   " st)
    (modify-syntax-entry ?~ " " st)
    st))

(defvar tme-db-mode-map
  (let ((km (make-sparse-keymap)))
    (define-key km "\t" 'tme-db-find-text)
    (define-key km "\n" 'tme-db-next-field)
    (define-key km "\M-\t" 'tme-db-complete-string)
    (define-key km [(meta tab)] 'tme-db-complete-key)
    (define-key km "\C-c\"" 'tme-db-remove-delimiters)
    (define-key km "\C-c{" 'tme-db-remove-delimiters)
    (define-key km "\C-c}" 'tme-db-remove-delimiters)
    (define-key km "\C-c\C-c" 'tme-db-clean-entry)
    (define-key km "\C-c\C-q" 'tme-db-fill-entry)
    (define-key km "\C-c?" 'tme-db-print-help-message)
    (define-key km "\C-c\C-p" 'tme-db-pop-previous)
    (define-key km "\C-c\C-n" 'tme-db-pop-next)
    (define-key km "\C-c\C-k" 'tme-db-kill-field)
    (define-key km "\C-c\M-k" 'tme-db-copy-field-as-kill)
    (define-key km "\C-c\C-w" 'tme-db-kill-entry)
    (define-key km "\C-c\M-w" 'tme-db-copy-entry-as-kill)
    (define-key km "\C-c\C-y" 'tme-db-yank)
    (define-key km "\C-c\M-y" 'tme-db-yank-pop)
    (define-key km "\C-c\C-d" 'tme-db-empty-field)
    (define-key km "\C-c\C-f" 'tme-db-make-field)
    (define-key km "\C-c$" 'tme-db-ispell-abstract)
    (define-key km "\M-\C-a" 'tme-db-beginning-of-entry)
    (define-key km "\M-\C-e" 'tme-db-end-of-entry)
    (define-key km "\C-\M-l" 'tme-db-reposition-window)
    (define-key km "\C-\M-h" 'tme-db-mark-entry)
    (define-key km "\C-c\C-b" 'tme-db-entry)
    (define-key km "\C-c\C-t" 'tme-db-hide-entry-bodies)
    (define-key km "\C-c\C-rn" 'tme-db-narrow-to-entry)
    (define-key km "\C-c\C-rw" 'widen)
    (define-key km "\C-c\C-o" 'tme-db-remove-OPT-or-ALT)
    (define-key km "\C-c\C-e\C-p" 'tme-db-Project)
    (define-key km "\C-c\C-e\M-p" 'tme-db-Preamble)
    (define-key km "\C-c\C-e\C-s" 'tme-db-String)
    km))

(easy-menu-define
 tme-db-edit-menu tme-db-mode-map "TME-db-Edit Menu in TME-db mode"
 '("TDB-Edit"
   ("Moving inside an Entry"
    ["End of Field" tme-db-find-text t]
    ["Next Field" tme-db-next-field t]
    ["Beginning of Entry" tme-db-beginning-of-entry t]
    ["End of Entry" tme-db-end-of-entry t])
   ("Operating on Current Entry"
    ["Fill Entry" tme-db-fill-entry t]
    ["Clean Entry" tme-db-clean-entry t]
    "--"
    ["Kill Entry" tme-db-kill-entry t]
    ["Copy Entry to Kill Ring" tme-db-copy-entry-as-kill t]
    ["Paste Most Recently Killed Entry" tme-db-yank t]
    ["Paste Previously Killed Entry" tme-db-yank-pop t]
    "--"
    ["Ispell Entry" tme-db-ispell-entry t]
    ["Ispell Entry Abstract" tme-db-ispell-abstract t]
    ["Narrow to Entry" tme-db-narrow-to-entry t])
   ("Operating on Current Field"
    ["Remove Delimiters" tme-db-remove-delimiters t]
    ["Remove OPT or ALT Prefix" tme-db-remove-OPT-or-ALT t]
    ["Clear Field" tme-db-empty-field t]
    "--"
    ["Kill Field" tme-db-kill-field t]
    ["Copy Field to Kill Ring" tme-db-copy-field-as-kill t]
    ["Paste Most Recently Killed Field" tme-db-yank t]
    ["Paste Previously Killed Field" tme-db-yank-pop t]
    "--"
    ["Make New Field" tme-db-make-field t]
    "--"
    ["Snatch from Similar Following Field" tme-db-pop-next t]
    ["Snatch from Similar Preceding Field" tme-db-pop-previous t]
    "--"
    ["String Complete" tme-db-complete-string t]
    ["Key Complete" tme-db-complete-key t]
    "--"
    ["Help about Current Field" tme-db-print-help-message t])
   ("Operating on Buffer or Region"
    ;;["Validate Entries" tme-db-validate t]
    ["Sort Entries" tme-db-sort-buffer t]
    ["Reformat Entries" tme-db-reformat t]
    ["Hide Entry Bodies" tme-db-hide-entry-bodies t]
    ["Count Entries" tme-db-count-entries t])
   ("Miscellaneous"
    ["Convert Alien Buffer" tme-db-convert-alien t]
    ["Submit Bug Report" tme-db-submit-bug-report t])))

(easy-menu-define
 tme-db-entry-menu tme-db-mode-map "Entry-Types Menu in TME-db mode"
 (list "TDB-Entries"
       ["Project" tme-db-Project t]
       ["String" tme-db-String t]
       ["Preamble" tme-db-Preamble t]))


;; Bug Reporting

(defconst
  tme-db-maintainer-address "Dirk Herrmann <D.Herrmann@tu-bs.de>")
;; current maintainer


;; Internal Variables

(defvar tme-db-pop-previous-search-point nil)
;; Next point where tme-db-pop-previous starts looking for a similar
;; entry.

(defvar tme-db-pop-next-search-point nil)
;; Next point where tme-db-pop-next starts looking for a similar entry.

(defvar tme-db-field-kill-ring nil)
;; Ring of least recently killed fields. At most
;; tme-db-field-kill-ring-max items are kept here.

(defvar tme-db-field-kill-ring-yank-pointer nil)
;; The tail of tme-db-field-kill-ring whose car is the last item yanked.

(defvar tme-db-entry-kill-ring nil)
;; Ring of least recently killed entries. At most
;; tme-db-entry-kill-ring-max items are kept here.

(defvar tme-db-entry-kill-ring-yank-pointer nil)
;; The tail of tme-db-entry-kill-ring whose car is the last item yanked.

(defvar tme-db-last-kill-command nil)
;; Holds the type of the last kill command (either 'field or 'entry)

(defvar tme-db-strings nil)
;; Candidates for tme-db-complete-string. Initialized from
;; tme-db-predefined-strings and tme-db-string-files. This variable is
;; buffer-local.
(make-variable-buffer-local 'tme-db-strings)

(defvar tme-db-keys nil)
;; Candidates for TAB completion when entering a reference key using
;; the minibuffer. Also used for tme-db-complete-key. Initialized in
;; tme-db-mode and updated for each new entry. This variable is
;; buffer-local.
(make-variable-buffer-local 'tme-db-keys)

(defvar tme-db-buffer-last-parsed-tick nil)
;; Remembers the value returned by buffer-modified-tick when buffer
;; was parsed for keys the last time.
(make-variable-buffer-local 'tme-db-buffer-last-parsed-tick)

(defvar tme-db-parse-idle-timer nil)
;; Stores if timer is already installed

(defvar tme-db-progress-lastperc nil)
;; Holds the last reported percentage for the progress message

(defvar tme-db-progress-lastmes nil)
;; Holds the last reported progress message

(defvar tme-db-progress-interval nil)
;; Holds the chosen interval

(defvar tme-db-key-history nil)
;; Used to store the history list for reading keys

(defvar tme-db-entry-type-history nil)
;; Used to store the history list for reading entry types

(defvar tme-db-field-history nil)
;; Used to store the history list for reading field names

(defvar tme-db-reformat-previous-options nil)
;; Used to store the last reformat options given

(defvar tme-db-reformat-previous-labels nil)
;; Used to store the last reformat label option given


;; Functions to Parse the TME-db Entries

(defun tme-db-cfield (name text)
  ;; Create a regexp for a TME-db field of name NAME and text TEXT.
  (concat ",[ \t\n]*\\("
	  name
	  "\\)[ \t\n]*=[ \t\n]*\\("
	  text
	  "\\)"))
(defconst tme-db-name-in-cfield 1)
;; The regexp subexpression number of the name part in tme-db-cfield.

(defconst tme-db-text-in-cfield 2)
;; The regexp subexpression number of the text part in tme-db-cfield.

(defconst tme-db-reference-type "@[^\"#%'(),={} \t\n0-9][^\"#%'(),={} \t\n]*")
;; Regexp defining the type part of a TME-db reference entry (almost
;; the same as tme-db-field-name)

(defconst tme-db-reference-key "[][a-z0-9.:;?!`'/*@+=|()<>&_^$-]+")
;; Regexp defining the label part of a TME-db reference entry

(defconst tme-db-field-name "[^\"#%'(),={} \t\n0-9][^\"#%'(),={} \t\n]*")
;; Regexp defining the name part of a TME-db field (almost the same as
;; tme-db-reference-type)

(defconst tme-db-field-const "[][a-z0-9.:;?!`'/*@+=|<>&_^$-]+")
;; Regexp defining a tme-db field constant

(defconst tme-db-field-string-part-not-braced
  "[^{}]")
;; Match field string part without braces

(defconst tme-db-field-string-part-no-inner-braces
  (concat
   "{"
   tme-db-field-string-part-not-braced "*"
   "}"))
;; Match field string part with no inner braces

(defconst tme-db-field-string-part-1-inner-brace
  (concat
   "{"
   "\\("
     tme-db-field-string-part-not-braced 
     "\\|"
     "\\(" tme-db-field-string-part-no-inner-braces "\\)"
   "\\)*"
   "}"))
;; Match field string part with at most 1 inner brace

(defconst tme-db-field-string-part-2-inner-braces
  (concat
   "{"
   "\\("
     tme-db-field-string-part-not-braced
     "\\|"
     "\\(" tme-db-field-string-part-1-inner-brace "\\)"
   "\\)*"
   "}"))
;; Match field string part with at most 2 inner braces

(defconst tme-db-field-string-part-3-inner-braces
  (concat
   "{"
   "\\("
     tme-db-field-string-part-not-braced
     "\\|"
     "\\(" tme-db-field-string-part-2-inner-braces "\\)"
   "\\)*"
   "}"))
;; Match field string part with at most 3 inner braces

(defconst tme-db-field-string-braced
  tme-db-field-string-part-3-inner-braces)
;; Match braced field string with inner nesting level of braces at most 3

(defconst tme-db-field-string-quoted
  (concat
   "\""
   "\\("
     "[^\"\\]"      ;; every character except quote or backslash
     "\\|"
;;     "\\(" "\"[a-z-]" "\\)"  ;; a quote followed by a letter or dash 
;;     "\\|"
;; last two lines commented out until lines like
;;   author = "Stefan Sch"of"
;; are supported by TME-db
     "\\(" "\\\\\\(.\\|\n\\)"  "\\)" ;; a backslash followed by any character
   "\\)*"
   "\""))
;; Match quoted field string

(defconst tme-db-field-string
  (concat
   "\\(" tme-db-field-string-braced "\\)"
   "\\|"
   "\\(" tme-db-field-string-quoted "\\)"))
;; Match a braced or quoted string

(defconst tme-db-field-string-or-const
  (concat tme-db-field-const "\\|" tme-db-field-string))
;; Match either tme-db-field-string or tme-db-field-const.

(defconst tme-db-field-text
  (concat
    "\\(" tme-db-field-string-or-const "\\)"
    "\\([ \t\n]+#[ \t\n]+\\(" tme-db-field-string-or-const "\\)\\)*"))
;; Regexp defining the text part of a TME-db field: either a string,
;; or an empty string, or a constant followed by one or more # /
;; constant pairs.

(defconst tme-db-field
  (tme-db-cfield tme-db-field-name tme-db-field-text))
;; Regexp defining the format of a TME-db field.

(defconst tme-db-name-in-field tme-db-name-in-cfield)
;; The regexp subexpression number of the name part in TME-db-field.

(defconst tme-db-text-in-field tme-db-text-in-cfield)
;; The regexp subexpression number of the text part in TME-db-field.

(defconst tme-db-reference-head
  (concat "^[ \t]*\\("
	  tme-db-reference-type
	  "\\)[ \t]*[({][ \t]*\\("
	  tme-db-reference-key
	  "\\)"))
;; Regexp defining format of the header line of a TME-db reference
;; entry.

(defconst tme-db-reference-maybe-empty-head
  (concat tme-db-reference-head "?"))
;; Regexp defining format of the header line of a maybe empty
;; TME-db reference entry (without reference key).

(defconst tme-db-type-in-head 1)
;; The regexp subexpression number of the type part in
;; tme-db-reference-head.

(defconst tme-db-key-in-head 2)
;; The regexp subexpression number of the key part in
;; tme-db-reference-head.

(defconst tme-db-reference-infix (concat "[ \t\n]*" tme-db-field))
;; Regexp defining the (repeatable) infix of a tme-db reference

(defconst tme-db-reference-postfix "[ \t\n]*,?[ \t\n]*[})]")
;; Regexp defining the postfix of a tme-db reference

(defconst tme-db-key-in-reference tme-db-key-in-head)
;; The regexp subexpression number of the key part in
;; tme-db-reference.

(defconst tme-db-string
  (concat "^[ \t]*@string[ \t\n]*[({][ \t\n]*\\("
	  tme-db-reference-key
	  "\\)[ \t\n]*=[ \t\n]*\\("
	  tme-db-field-text
	  "\\)[ \t\n]*[})]"))
;; Regexp defining the format of a TME-db string entry.

(defconst tme-db-key-in-string 1)
;; The regexp subexpression of the name part in tme-db-string.

(defconst tme-db-text-in-string 2)
;; The regexp subexpression of the text part in tme-db-string.

(defvar tme-db-font-lock-keywords
  (list
   ;; reference type and reference label
   (list tme-db-reference-maybe-empty-head
         (list tme-db-type-in-head 'font-lock-function-name-face)
         (list tme-db-key-in-head 'font-lock-reference-face nil t))
   ;; comments
   (list 
    (concat "^\\([ \t]*" tme-db-comment-start ".*\\)$")
    1 'font-lock-comment-face)
   ;; optional field names (treated as comments)
   (list
    (concat "^[ \t]*\\(OPT" tme-db-field-name "\\)[ \t]*=")
    1 'font-lock-comment-face)
   ;; field names
   (list (concat "^[ \t]*\\(" tme-db-field-name "\\)[ \t]*=")
         1 'font-lock-variable-name-face)
   "*Default expressions to highlight in TME-db mode."))
;; now all needed patterns are defined


;; Helper Functions

(defun tme-db-delete-whitespace ()
  ;; Delete all whitespace starting at point
  (if (looking-at "[ \t\n]+")
      (delete-region (point) (match-end 0))))

(defun tme-db-current-line ()
  ;; this computes line number of point regardless whether the buffer
  ;; is narrowed
  (+ (count-lines 1 (point))
     (if (equal (current-column) 0) 1 0)))

(defun tme-db-member-of-regexp (string list)
  ;; Return non-nil if STRING is exactly matched by an element of
  ;; LIST. The value is actually the tail of LIST whose
  ;; car matches STRING.
  (let (case-fold-search)
    (while
        (and list (not (string-match (concat "^" (car list) "$") string)))
      (setq list (cdr list)))
    list))

(defun tme-db-assoc-of-regexp (string alist)
  ;; Return non-nil if STRING is exactly matched by the car of an
  ;; element of LIST (case ignored). The value is actually the element
  ;; of LIST whose car matches STRING.
  (let ((case-fold-search t))
    (while
        (and alist
             (not (string-match (concat "^" (car (car alist)) "$") string)))
      (setq alist (cdr alist)))
    (car alist)))

(defun tme-db-skip-to-valid-entry (&optional backward)
  ;; If not at beginning of valid TME-db entry, move to beginning of
  ;; the next valid one. With argument backward non-nil, move backward
  ;; to beginning of previous valid one. A valid entry is a
  ;; syntactical correct one with type contained in
  ;; tme-db-entry-field-alist or, if tme-db-sort-ignore-string-entries
  ;; is nil, a syntactical correct string entry.
  (let ((case-fold-search t)
        (valid-tme-db-entry
         (concat
          "@[ \t]*\\("
          (mapconcat
           (lambda (type)
             (concat "\\(" (car type) "\\)"))
           tme-db-entry-field-alist
           "\\|")
          "\\)"))
        found)
    (while (and (not found)
                (not (if backward
                         (bobp)
                       (eobp))))
      (let ((pnt (point)))
        (cond
         ((looking-at valid-tme-db-entry)
          (if (and
               (tme-db-search-reference nil nil t)
               (equal (match-beginning 0) pnt))
              (setq found t)))
         ((and (not tme-db-sort-ignore-string-entries)
               (looking-at tme-db-string))
          (setq found t)))
        (if found
            (goto-char pnt)
          (if backward
              (progn
                (goto-char (1- pnt))
                (if (re-search-backward "^[ \t]*\\(@\\)" nil 'move)
                    (goto-char (match-beginning 1))))
            (goto-char (1+ pnt))
            (if (re-search-forward "^[ \t]*@" nil 'move)
                (forward-char -1))))))))

(defun tme-db-map-entries (fun)
  ;; Call FUN for each TME-db entry starting with the current. Do this
  ;; to the end of the file. FUN is called with one argument, the key
  ;; of the entry, and with point inside the entry. If
  ;; tme-db-sort-ignore-string-entries is non-nil, FUN will not be
  ;; called for @String entries.
  (let ((case-fold-search t))
    (tme-db-beginning-of-entry)
    (while (re-search-forward tme-db-reference-maybe-empty-head nil t)
      (let ((pnt (point))
            (reference-type
             (downcase (buffer-substring-no-properties
                        (1+ (match-beginning tme-db-type-in-head))
                        (match-end tme-db-type-in-head))))
            (reference-key
             (if (match-beginning tme-db-key-in-head)
                 (buffer-substring-no-properties
                  (match-beginning tme-db-key-in-head)
                  (match-end tme-db-key-in-head))
               "")))
        (if (or
             (and
              (not tme-db-sort-ignore-string-entries)
              (string-equal "string" (downcase reference-type)))
             (assoc-ignore-case reference-type tme-db-entry-field-alist))
            (funcall fun reference-key))
        (goto-char pnt)
        (tme-db-end-of-entry)))))

(defun tme-db-progress-message (&optional flag interval)
  ;; echos a message about progress of current buffer
  ;; if flag is a string, the message is initialized (in this case a
  ;; value for INTERVAL may be given as well (if not this is set to 5))
  ;; if flag is done, the message is deinitialized
  ;; if flag is absent, a message is echoed if point was incremented
  ;; at least INTERVAL percent since last message was echoed
  (let* ((size (- (point-max) (point-min)))
         (perc (if (= size 0)
                   100
                 (/ (* 100 (- (point) (point-min))) size))))
    (if (or (and (not flag)
                 (>= perc
                     (+ tme-db-progress-interval tme-db-progress-lastperc)))
            (stringp flag))
        (progn
          (if (stringp flag)
              (progn
                (setq tme-db-progress-lastmes flag)
                (if interval
                    (setq tme-db-progress-interval interval)
                  (setq tme-db-progress-interval 5))))
          (setq tme-db-progress-lastperc perc)
          (message (concat tme-db-progress-lastmes " (%d%%)") perc))
      (if (equal flag 'done)
          (progn
            (message (concat tme-db-progress-lastmes " (done)"))
            (setq tme-db-progress-lastmes nil))))))


(defun tme-db-field-left-delimiter ()
  ;; returns a string dependent on tme-db-field-delimiters
  (if (equal tme-db-field-delimiters 'braces)
      "{"
    "\""))

(defun tme-db-field-right-delimiter ()
  ;; returns a string dependent on tme-db-field-delimiters
  (if (equal tme-db-field-delimiters 'braces)
      "}"
    "\""))

(defun tme-db-entry-left-delimiter ()
  ;; returns a string dependent on tme-db-field-delimiters
  (if (equal tme-db-entry-delimiters 'braces)
      "{"
    "("))

(defun tme-db-entry-right-delimiter ()
  ;; returns a string dependent on tme-db-field-delimiters
  (if (equal tme-db-entry-delimiters 'braces)
      "}"
    ")"))

(defun tme-db-search-reference
  (empty-head &optional bound noerror backward)
  ;; A helper function necessary since the failure stack size limit for
  ;; regexps was reduced in emacs 19.32.
  ;; It searches for a TME-db reference (maybe with an empty head if
  ;; EMPTY-HEAD is t).
  ;; BOUND and NOERROR are exactly as in re-search-forward. If
  ;; BACKWARD is non-nil, search is done in reverse direction. After
  ;; call to this function MATCH-BEGINNING and MATCH-END functions are
  ;; defined, but only for the head part of the reference (especially
  ;; (match-end 0) just gives the end of the head part).
  (let ((pnt (point))
        (prefix (if empty-head
                    tme-db-reference-maybe-empty-head
                  tme-db-reference-head))
        (infix tme-db-reference-infix)
        (postfix tme-db-reference-postfix))
    (if backward
        (let (found)
          (while (and (not found)
                      (re-search-backward prefix bound noerror))
            (setq found (tme-db-search-reference empty-head pnt t)))
          (if found
              (goto-char (match-beginning 0))
            (if (equal noerror nil)
                ;; yell
                (error "Search of TME-db reference failed."))
            (if (equal noerror t)
                ;; don't move
                (goto-char pnt))
            nil))
      (let ((limit (if bound bound (point-max)))
            md
            found)
        (while (and (not found)
                    (re-search-forward prefix bound noerror))
          (setq md (match-data))
          ;; save match-data of prefix regexp
          (let ((entry-closer
                 (if (save-excursion
                       (goto-char (match-end tme-db-type-in-head))
                       (looking-at "[ \t]*("))
                     ;; entry opened with parenthesis
                     ")"
                   "}")))
            (while (and
                    (looking-at infix)
                    (<= (match-end 0) limit))
              (goto-char (match-end 0)))
            ;; This matches the infix* part. The AND construction assures
            ;; that BOUND is respected.
            (if (and (looking-at postfix)
                     (string-equal
                      (buffer-substring-no-properties
                       (1- (match-end 0)) (match-end 0))
                      entry-closer)
                     (<= (match-end 0) limit))
                (progn
                  (re-search-forward postfix)
                  (setq found t)))))
        (if found
            (progn
              (store-match-data md)
              ;; to set match-beginning/end again
              (point))
          (if (equal noerror nil)
              ;; yell
              (error "Search of TME-db reference failed."))
          (if (equal noerror t)
              ;; don't move
              (goto-char pnt))
          nil)))))

(defun tme-db-flash-head ()
  ;; Flash at TME-db reference head before point, if exists.
  (let ((case-fold-search t)
        flash)
    (cond ((re-search-backward tme-db-reference-head nil t)
	   (goto-char (match-beginning tme-db-type-in-head))
	   (setq flash (match-end tme-db-key-in-reference)))
	  (t
	   (end-of-line)
	   (skip-chars-backward " \t")
	   (setq flash (point))
	   (beginning-of-line)
	   (skip-chars-forward " \t")))
    (if (pos-visible-in-window-p (point))
	(sit-for 1)
      (message "From: %s"
	       (buffer-substring (point) flash)))))

(defun tme-db-make-optional-field (e-t)
  "Makes an optional field named E-T in current TME-db entry."
  (if (consp e-t)
      (tme-db-make-field (cons (concat "OPT" (car e-t)) (cdr e-t)))
    (tme-db-make-field (concat "OPT" e-t))))

(defun tme-db-move-outside-of-entry ()
  ;; Make sure we are outside of a TME-db entry.
  (tme-db-end-of-entry)
  (skip-chars-forward " \t\n"))

(defun tme-db-beginning-of-first-entry ()
  ;; Go to the beginning of the first TME-db entry in buffer. Return
  ;; point.
  (goto-char (point-min))
  (if (re-search-forward "^[ \t]*@" nil 'move)
      (beginning-of-line))
  (point))

(defun tme-db-beginning-of-last-entry ()
  ;; Go to the beginning of the last TME-db entry in buffer.
  (goto-char (point-max))
  (if (re-search-backward "^[ \t]*@" nil 'move)
      (beginning-of-line))
  (point))

(defun tme-db-inside-field ()
  ;; Try to avoid point being at end of a TME-db field.
  (end-of-line)
  (skip-chars-backward " \t")
  (cond ((= (preceding-char) ?,)
	 (forward-char -2)))
  (cond ((or
          (= (preceding-char) ?})
          (= (preceding-char) ?\"))
         (forward-char -1))))

(defun tme-db-enclosing-field (&optional noerr)
  ;; Search for TME-db field enclosing point. Point moves to end of
  ;; field. Use match-beginning and match-end to parse the field. If
  ;; NOERR is non-nil, no error is signalled. In this case, t is
  ;; returned on success, nil otherwise.
  (let ((case-fold-search t)
        (old-point (point))
        (boe (tme-db-beginning-of-entry))
        (success t))
    (goto-char old-point)
    (if (not (re-search-backward
              (tme-db-cfield tme-db-field-name "")
              boe t))
        ;; Search possible beginning of field
        (progn
          (goto-char old-point)
          (if noerr
              (setq success nil)
            (error "Can't find enclosing TME-db field.")))
      (if (or (not (re-search-forward tme-db-field nil t))
              (< (match-end 0) old-point)
              (> (match-beginning 0) old-point))
          (progn
            (goto-char old-point)
            (if noerr
                (setq success nil)
              (error "Can't find enclosing TME-db field.")))))
    success))

(defun tme-db-enclosing-reference-maybe-empty-head ()
  ;; Search for TME-db reference enclosing point. Point moves to
  ;; end of reference. Beginning (but not end) of reference is given
  ;; by (match-beginning 0).
  (let ((case-fold-search t)
        (old-point (point)))
    (if (not
         (re-search-backward
          tme-db-reference-maybe-empty-head nil t))
        (progn
          (error "Can't find enclosing TME-db reference.")
          (goto-char old-point)))
    (goto-char (match-beginning tme-db-type-in-head))
    (if (not
         (tme-db-search-reference t nil t))
        (progn
          (error "Can't find enclosing TME-db reference.")
          (goto-char old-point)))))

(defun tme-db-insert-current-kill (n)
  (if (not tme-db-last-kill-command)
      (error "TME-db kill ring is empty.")
    (let* ((kr (if (equal tme-db-last-kill-command 'field)
                   'tme-db-field-kill-ring
                 'tme-db-entry-kill-ring))
           (kryp (if (equal tme-db-last-kill-command 'field)
                     'tme-db-field-kill-ring-yank-pointer
                   'tme-db-entry-kill-ring-yank-pointer))
           (ARGth-kill-element
            (nthcdr
             (mod (- n (length (eval kryp))) (length (eval kr)))
             (eval kr)))
           (current (car (set kryp ARGth-kill-element))))
      (cond
       ((equal tme-db-last-kill-command 'field)
        (let (tme-db-help-message)
          (tme-db-find-text nil t)
          (if (looking-at "[}\"]")
              (forward-char)))
        (set-mark (point))
        (message "Mark set")
        (tme-db-make-field (list (elt current 1) nil (elt current 2)) t))
       ((equal tme-db-last-kill-command 'entry)
        (if (not (eobp))
            (tme-db-beginning-of-entry))
        (set-mark (point))
        (message "Mark set")
        (insert (elt current 1)))
       (t
        (error
         "Unknown tag field: %s. Please submit a bug report."
         tme-db-last-kill-command))))))

(defun tme-db-format-entry ()
  ;; Helper function for tme-db-clean-entry. Formats current entry
  ;; according to variable tme-db-entry-format.
  (let ((case-fold-search t)
        (beg (point))
        (start (tme-db-beginning-of-entry))
        crossref-there
        alternatives-there
        non-empty-alternative)
    (let ((end (copy-marker (tme-db-end-of-entry))))
      (if (equal start (marker-position end))
          (error "Not on a known TME-db entry.")
        (goto-char start)
        (while (re-search-forward tme-db-field end t)
          ;; determine if reference has crossref entry and if at least
          ;; one alternative is non-empty
          (let ((begin-name (match-beginning tme-db-name-in-field))
                (end-name (match-end tme-db-name-in-field))
                (begin-text (match-beginning tme-db-text-in-field))
                (end-text (match-end tme-db-text-in-field)))
            (goto-char begin-name)
            (if (looking-at "ALT")
                (progn
                  (setq alternatives-there t)
                  (goto-char begin-text)
                  (if (not (looking-at "\\(\"\"\\)\\|\\({}\\)"))
                      (setq non-empty-alternative t))))
            (if (string-match
                 "\\(OPT\\)?crossref"
                 (buffer-substring-no-properties begin-name end-name))
                (progn
                  (setq
                   crossref-there
                   (buffer-substring-no-properties
                    (1+ begin-text) (1- end-text)))
                  (if (equal crossref-there "")
                      (setq crossref-there nil))))))
        (if (and alternatives-there
                 (not non-empty-alternative))
            (progn
              (goto-char beg)
              (error "All alternatives are empty.")))
        (goto-char start)
        (re-search-forward tme-db-reference-type end)
        (let* ((begin-type (1+ (match-beginning 0)))
               (end-type (match-end 0))
               (reference-type
                (downcase
                 (buffer-substring-no-properties begin-type end-type)))
               (entry-list
                (assoc-ignore-case reference-type tme-db-entry-field-alist))
               (req (elt (elt entry-list 1) 0))
               (creq (elt (elt entry-list 2) 0))
               (format (if (equal tme-db-entry-format t)
                           '(realign opts-or-alts numerical-fields
                                     last-comma page-dashes delimiters
                                     unify-case inherit-booktitle)
                         tme-db-entry-format))
               field-done)
          (if (memq 'unify-case format)
              (progn
                (delete-region begin-type end-type)
                (insert (car entry-list))))
          (if (memq 'delimiters format)
              (progn
                (goto-char end-type)
                (skip-chars-forward " \t\n")
                (delete-char 1)
                (insert (tme-db-entry-left-delimiter))))
          (goto-char start)
          (while (re-search-forward tme-db-field end t)
            (let* ((begin-field
                    (copy-marker (match-beginning 0)))
                   (end-field
                    (copy-marker (match-end 0)))
                   (begin-name
                    (copy-marker (match-beginning tme-db-name-in-field)))
                   (end-name
                    (copy-marker (match-end tme-db-name-in-field)))
                   (begin-text
                    (copy-marker (match-beginning tme-db-text-in-field)))
                   (end-text
                    (copy-marker (match-end tme-db-text-in-field)))
                   (field-name
                    (buffer-substring-no-properties
                     (if (string-match
                          "^OPT\\|ALT$"
                          (buffer-substring-no-properties
                           begin-name (+ begin-name 3)))
                         (+ begin-name 3)
                       begin-name)
                     end-name)))
              (cond
               ((and
                 (memq 'opts-or-alts format)
                 (progn (goto-char begin-name)
                        (looking-at "OPT\\|ALT")))
                (goto-char begin-text)
                (if (looking-at "\\(\"\"\\)\\|\\({}\\)")
                    ;; empty: delete whole field if really optional
                    ;; (missing crossref handled) or complain
                    (if (and
                         (progn
                           (goto-char begin-name)
                           (looking-at "OPT"))
                         (not crossref-there)
                         (assoc-ignore-case field-name req))
                        ;; field is not really optional
                        (progn
                          (goto-char begin-name)
                          (tme-db-remove-OPT-or-ALT)
                          (error
                           "Mandatory field ``%s'' is empty." field-name))
                      ;; field is optional
                      (delete-region begin-field end-field))
                  ;; otherwise: not empty, delete "OPT"
                  (goto-char begin-name)
                  (tme-db-remove-OPT-or-ALT)))
               ((and
                 (memq 'numerical-fields format)
                 (progn
                   (goto-char begin-text)
                   (looking-at "\\(\"[0-9]+\"\\)\\|\\({[0-9]+}\\)")))
                (goto-char end-text)
                (delete-char -1)
                (goto-char begin-text)
                (delete-char 1))
               (t
                (if (memq 'delimiters format)
                    (progn
                      (goto-char begin-text)
                      (while (and
                              (<= (point) end-text)
                              (re-search-forward
                               tme-db-field-string-or-const end-text t))
                        (let ((end (point)))
                          (goto-char (match-beginning 0))
                          (if (or
                               (and
                                (equal tme-db-field-delimiters 'double-quotes)
                                (looking-at tme-db-field-string-braced))
                               (and
                                (equal tme-db-field-delimiters 'braces)
                                (looking-at tme-db-field-string-quoted)))
                              (progn
                                (goto-char (match-end 0))
                                (delete-backward-char 1)
                                (insert (tme-db-field-right-delimiter))
                                (goto-char (match-beginning 0))
                                (delete-char 1)
                                (insert (tme-db-field-left-delimiter))))
                          (goto-char end)))))
                (if (and
                     (memq 'page-dashes format)
                     (string-match "^\\(OPT\\)?pages$" (downcase field-name))
                     (progn
                       (goto-char begin-text)
                       (looking-at
                        "\\([\"{][0-9]+\\)[ \t\n]*--?[ \t\n]*\\([0-9]+[\"}]\\)")))
                    (replace-match "\\1-\\2"))
                (if (and
                     (memq 'inherit-booktitle format)
                     (equal (downcase field-name) "booktitle")
                     (progn
                       (goto-char begin-text)
                       (looking-at "\\(\"\"\\)\\|\\({}\\)"))
                     crossref-there
                     (not (tme-db-find-entry-location crossref-there t)))
                    ;; booktitle field empty and crossref entry found
                    ;; --> insert title field of crossreferenced entry if there
                    (let ((end-of-crefd-entry (tme-db-end-of-entry)))
                      (tme-db-beginning-of-entry)
                      (if (re-search-forward
                           (tme-db-cfield "title" tme-db-field-text)
                           end-of-crefd-entry t)
                          (progn
                            (goto-char begin-text)
                            (forward-char)
                            (insert
                             (buffer-substring-no-properties 
                              (1+ (match-beginning tme-db-text-in-field))
                              (1- (match-end tme-db-text-in-field))))))))
                (if (progn
                      (goto-char begin-text)
                      (looking-at "\\(\"\"\\)\\|\\({}\\)"))
                    ;; if empty field, complain
                    (progn
                      (forward-char)
                      (if (or (and
                               crossref-there
                               (assoc-ignore-case
                                field-name creq))
                              (and
                               (not crossref-there)
                               (assoc-ignore-case
                                field-name req)))
                          (error
                           "Mandatory field ``%s'' is empty." field-name))))
                (if (memq 'unify-case format)
                    (let* ((fl
                            (car (cdr (assoc-ignore-case
                                       reference-type
                                       tme-db-entry-field-alist))))
                           (field-list
                            (append
                             (elt fl 0)
                             (elt fl 1)
                             tme-db-user-optional-fields))
                           (new-field-name
                            (car
                             (assoc-ignore-case field-name field-list))))
                      (goto-char begin-name)
                      (if new-field-name
                          (progn
                            (delete-region begin-name end-name)
                            (insert new-field-name))
                        (downcase-region begin-name end-name))))
                (setq field-done t)))
              (if (not field-done)
                  (goto-char begin-field)
                (setq field-done nil)
                (goto-char end-field))))
          (if (looking-at (tme-db-field-right-delimiter))
              (forward-char))
          (if (memq 'last-comma format)
              (cond ((and
                      tme-db-comma-after-last-field
                      (not (looking-at ",")))
                     (insert ","))
                    ((and
                      (not tme-db-comma-after-last-field)
                      (looking-at ","))
                     (delete-char 1))))
          (if (looking-at ",")
              (forward-char))
          (if (memq 'delimiters format)
              (progn
                (skip-chars-forward " \t\n")
                (delete-char 1)
                (insert (tme-db-entry-right-delimiter))))
          (if (memq 'realign format)
              (tme-db-fill-entry)))))))

(defun tme-db-autokey-change (string change-list)
  ;; Returns a string where some regexps are changed according to
  ;; change-list. Every item of change-list is an (old-regexp
  ;; new-string) pair.
  (let (case-fold-search
        (return-string string)
        (index 0)
        (len (length change-list))
        change-item)
    (while (< index len)
      (setq change-item (elt change-list index))
      (while (string-match (car change-item) return-string)
        (setq
         return-string
         (concat (substring return-string 0 (match-beginning 0))
                 (elt change-item 1)
                 (substring return-string (match-end 0)))))
      (setq index (1+ index)))
    return-string))

(defun tme-db-autokey-abbrev (string len)
  ;; Returns an abbreviation of string with at least len
  ;; characters. String is aborted only after a consonant or at the
  ;; word end. If len is not a number, string is returned unchanged.
  (cond ((or
          (not (numberp len))
          (<= (length string) len))
         string)
        ((equal len 0)
         "")
        (t
         (let* ((case-fold-search t)
                (abort-char
                 (string-match "[^aeiou]" string (1- len))))
           (if abort-char
               (substring string 0 (1+ abort-char))
             string)))))

(defun tme-db-autokey-get-namefield (min max)
  ;; returns the contents of the name field of the current entry
  ;; does some modifications based on `tme-db-autokey-name-change-strings'
  ;; and removes newlines unconditionally
  (goto-char min)
  (let ((case-fold-search t))
    (if (re-search-forward
         (tme-db-cfield "\\(author\\)\\|\\(editor\\)" tme-db-field-text)
         max t)
        (tme-db-autokey-change
         (buffer-substring-no-properties
          (1+ (match-beginning (+ tme-db-text-in-cfield 2)))
          (1- (match-end (+ tme-db-text-in-cfield 2))))
         (append tme-db-autokey-name-change-strings '(("\n" " "))))
      "")))

(defun tme-db-autokey-get-names (namefield)
  ;; gathers all names in namefield into a list
  (let ((case-fold-search t)
        names)
    (while (not (equal namefield ""))
      (let (name)
        (if (string-match "[ \t\n]and[ \t\n]" namefield)
            (setq name (substring namefield 0 (match-beginning 0))
                  namefield (substring namefield (match-end 0)))
          (setq name namefield
                namefield ""))
        (setq names (append names (list name)))))
    names))

(defun tme-db-autokey-demangle-name (fullname)
  ;; gets the `last' part from a well-formed name
  (let* (case-fold-search
         (lastname
          (if (string-match "," fullname)
              ;; name is of the form "von Last, First" or
              ;; "von Last, Jr, First"
              ;; --> take only the part before the comma
              (let ((von-last
                     (substring fullname 0 (match-beginning 0))))
                (if (string-match "^[a-z]" von-last)
                    ;; von-last has a "von" part --> take the "last" part
                    (if (string-match "[ \t][A-Z]" von-last)
                        (substring von-last (1+ (match-beginning 0)))
                      (error
                       "Name %s is incorrectly formed" fullname))
                  ;; von-last has no "von" part --> take all
                  von-last))
            ;; name is of the form "First von Last"
            (if (string-match "[ \t]" fullname)
                ;; more than one token
                (if (string-match "[ \t][a-z]" fullname)
                    ;; there is a "von" part
                    ;; --> take everything after that
                    (if (string-match
                         "[ \t][A-Z]" fullname (match-end 0))
                        (substring fullname (1+ (match-beginning 0)))
                      (error
                       "Name %s is incorrectly formed" fullname))
                  ;; there is no "von" part --> take only the last token
                  (if (string-match " [^ ]*$" fullname)
                      (substring fullname (1+ (match-beginning 0)))
                    (error "Name %s is incorrectly formed" fullname)))
              ;; only one token --> take it
              fullname)))
         (usename
          (if (string-match "[ \t]+" lastname)
              ;; lastname consists of two or more tokens
              ;; --> take only the first one
              (substring lastname 0 (match-beginning 0))
            lastname)))
    (funcall tme-db-autokey-name-case-convert usename)))

(defun tme-db-autokey-get-namelist (namefield)
  ;; gets namefield, performs abbreviations on the last parts, and
  ;; return results as a list
  (mapcar
   (lambda (fullname)
     (setq
      fullname (substring fullname (string-match "[^ \t]" fullname)))
     (tme-db-autokey-abbrev
      (tme-db-autokey-demangle-name fullname)
      tme-db-autokey-name-length))
   (tme-db-autokey-get-names namefield)))

(defun tme-db-autokey-get-yearfield (min max)
  ;; get year field from current or maybe crossreferenced entry
  (let ((case-fold-search t))
    (goto-char min)
    (if (re-search-forward
         (tme-db-cfield "year" tme-db-field-text) max t)
	(let ((year (buffer-substring-no-properties
		     (match-beginning tme-db-text-in-cfield)
		     (match-end tme-db-text-in-cfield))))
	  (string-match "[^0-9]*\\([0-9]+\\)" year)
	  (substring year (match-beginning 1) (match-end 1)))
      (if tme-db-autokey-year-use-crossref-entry
          (let ((crossref-field
                 (progn
                   (goto-char min)
                   (if (re-search-forward
                        (tme-db-cfield
                         "\\(OPT\\)?crossref" tme-db-field-text)
                        max t)
                       (buffer-substring-no-properties
                        (1+
                         (match-beginning (+ tme-db-text-in-cfield 1)))
                        (1-
                         (match-end (+ tme-db-text-in-cfield 1))))))))
            (if (not (tme-db-find-entry-location crossref-field t))
                (let ((end-of-crefd-entry (tme-db-end-of-entry)))
                  (tme-db-beginning-of-entry)
                  (if (re-search-forward
                       (tme-db-cfield "year" "[0-9]+")
                       end-of-crefd-entry t)
                      (buffer-substring-no-properties
                       (match-beginning tme-db-text-in-cfield)
                       (match-end tme-db-text-in-cfield))
                    ""))
              ""))
        ""))))

(defun tme-db-autokey-get-titlestring (min max)
  ;; get title field contents up to a terminator
  (let ((case-fold-search t))
    (let ((titlefield
           (progn
             (goto-char min)
             (if (re-search-forward
                  (tme-db-cfield "title" tme-db-field-text) max t)
                 (tme-db-autokey-change
                  (buffer-substring-no-properties
                   (match-beginning tme-db-text-in-cfield)
                   (match-end tme-db-text-in-cfield))
                  tme-db-autokey-titleword-change-strings)
               "")))
          (index 0)
          (numberofitems
           (length tme-db-autokey-title-terminators)))
      (while (< index numberofitems)
        (if (string-match
             (elt tme-db-autokey-title-terminators index) titlefield)
            (setq
             titlefield (substring titlefield 0 (match-beginning 0))))
        (setq index (1+ index)))
      titlefield)))

(defun tme-db-autokey-get-titles (titlestring)
  ;; gathers words from titlestring into a list. Ignores
  ;; specific words and uses only a specific amount of words.
  (let (case-fold-search
        titlewords
        titlewords-extra
        (counter 0))
    (while (and
            (not (equal titlestring ""))
            (or
             (not (numberp tme-db-autokey-titlewords))
             (< counter
                (+ tme-db-autokey-titlewords
                   tme-db-autokey-titlewords-stretch))))
      (if (string-match "\\b\\w+" titlestring)
          (let* ((end-match (match-end 0))
                 (titleword
		  (substring titlestring (match-beginning 0) end-match)))
	    (if (tme-db-member-of-regexp
		 titleword
		 tme-db-autokey-titleword-ignore)
		(setq counter (1- counter))
	      (setq 
	       titleword
	       (funcall tme-db-autokey-titleword-case-convert titleword))
	      (if (or (not (numberp tme-db-autokey-titlewords))
		      (< counter tme-db-autokey-titlewords))
		  (setq titlewords (append titlewords (list titleword)))
		(setq titlewords-extra
		      (append titlewords-extra (list titleword)))))
            (setq
             titlestring (substring titlestring end-match)))
        (setq titlestring ""))
      (setq counter (1+ counter)))
    (if (string-match "\\b\\w+" titlestring)
        titlewords
      (append titlewords titlewords-extra))))

(defun tme-db-autokey-get-titlelist (titlestring)
  ;; returns all words in titlestring as a list
  ;; does some abbreviation on the found words
  (mapcar
   (lambda (titleword)
     (let ((abbrev
            (tme-db-assoc-of-regexp
             titleword tme-db-autokey-titleword-abbrevs)))
       (if abbrev
           (elt abbrev 1)
         (tme-db-autokey-abbrev
          titleword
          tme-db-autokey-titleword-length))))
   (tme-db-autokey-get-titles titlestring)))

(defun tme-db-generate-autokey ()
  "Generates automatically a key from the author/editor and the title field.
This will only work for entries where each field begins on a separate line.
The generation algorithm works as follows:
 1. Use the value of `tme-db-autokey-prefix-string' as a prefix.
 2. If there is a non-empty author (preferred) or editor field,
    use it as the name part of the key.
 3. Change any substring found in
    `tme-db-autokey-name-change-strings' to the corresponding new
    one (see documentation of this variable for further detail).
 4. For every of at least first `tme-db-autokey-names' names in
    the name field, determine the last name. If there are maximal
    `tme-db-autokey-names' + `tme-db-autokey-names-stretch'
    names, all names are used.
 5. From every last name, take at least
    `tme-db-autokey-name-length' characters (abort only after a
    consonant or at a word end).
 6. Convert all last names according to the conversion function
    `tme-db-autokey-name-case-convert'.
 7. Build the name part of the key by concatenating all
    abbreviated last names with the string
    `tme-db-autokey-name-separator' between any two. If there are
    more names than are used in the name part, prepend the string
    contained in `tme-db-autokey-additional-names'.
 8. Build the year part of the key by truncating the contents of
    the year field to the rightmost `tme-db-autokey-year-length'
    digits (useful values are 2 and 4). If the year field is
    absent, but the entry has a valid crossref field and the
    variable `tme-db-autokey-year-use-crossref-entry' is non-nil,
    use the year field of the crossreferenced entry instead.
 9. For the title part of the key change the contents of the
    title field of the reference according to
    `tme-db-autokey-titleword-change-strings' to the
    corresponding new one (see documentation of this variable for
    further detail).
10. Abbreviate the result to the string up to (but not including)
    the first occurrence of a regexp matched by the items of
    `tme-db-autokey-title-terminators' and delete those words which
    appear in `tme-db-autokey-titleword-ignore'.
    Build the title part of the key by using at least the first
    `tme-db-autokey-titlewords' words from this
    abbreviated title. If the abbreviated title ends after
    maximal `tme-db-autokey-titlewords' +
    `tme-db-autokey-titlewords-stretch' words, all
    words from the abbreviated title are used.
11. Convert all used titlewords according to the conversion function
    `tme-db-autokey-titleword-case-convert'.
12. For every used title word that appears in
    `tme-db-autokey-titleword-abbrevs' use the corresponding
    abbreviation (see documentation of this variable for further
    detail).
13. From every title word not generated by an abbreviation, take
    at least `tme-db-autokey-titleword-length' characters (abort
    only after a consonant or at a word end).
14. Build the title part of the key by concatenating all
    abbreviated title words with the string
    `tme-db-autokey-titleword-separator' between any two.
15. At least, to get the key, concatenate
    `tme-db-autokey-prefix-string', the name part, the year part
    and the title part with `tme-db-autokey-name-year-separator'
    between the name part and the year part if both are non-empty
    and `tme-db-autokey-year-title-separator' between the year
    part and the title part if both are non-empty. If the year
    part is empty, but not the other two parts,
    `tme-db-autokey-year-title-separator' is used as well.
16. If the value of `tme-db-autokey-before-presentation-function'
    is non-nil, it must be a function taking one argument. This
    function is then called with the generated key as the
    argument. The return value of this function (a string) is
    used as the key.
17. If the value of `tme-db-autokey-edit-before-use' is non-nil,
    the key is then presented in the minibuffer to the user,
    where it can be edited. The key given by the user is then
    used.
"
  (let* ((pnt (point))
         (min (tme-db-beginning-of-entry))
         (max (tme-db-end-of-entry))
         (namefield (tme-db-autokey-get-namefield min max))
         (name-etal "")
         (namelist
          (let ((nl (tme-db-autokey-get-namelist namefield)))
            (if (or (not (numberp tme-db-autokey-names))
                    (<= (length nl)
                        (+ tme-db-autokey-names
                           tme-db-autokey-names-stretch)))
                nl
              (setq name-etal tme-db-autokey-additional-names)
              (let (nnl)
                (while (< (length nnl) tme-db-autokey-names)
                  (setq nnl (append nnl (list (car nl)))
                        nl (cdr nl)))
                nnl))))
         (namepart
          (concat
           (mapconcat (lambda (name) name)
                      namelist
                      tme-db-autokey-name-separator)
           name-etal))
         (yearfield (tme-db-autokey-get-yearfield min max))
         (yearpart
          (if (equal yearfield "")
              ""
            (substring
             yearfield
             (- (length yearfield) tme-db-autokey-year-length))))
         (titlestring (tme-db-autokey-get-titlestring min max))
         (titlelist (tme-db-autokey-get-titlelist titlestring))
         (titlepart
          (mapconcat
           (lambda (name) name)
           titlelist
           tme-db-autokey-titleword-separator))
         (autokey
          (concat
           tme-db-autokey-prefix-string
           namepart
           (if (not
                (or
                 (equal namepart "")
                 (equal yearpart "")))
               tme-db-autokey-name-year-separator)
           yearpart
           (if (not
                (or
                 (and
                  (equal namepart "")
                  (equal yearpart ""))
                 (equal titlepart "")))
               tme-db-autokey-year-title-separator)
           titlepart)))
    (if tme-db-autokey-before-presentation-function
        (setq
         autokey
         (funcall tme-db-autokey-before-presentation-function autokey)))
    (goto-char pnt)
    autokey))

(defun tme-db-parse-keys (add verbose &optional abortable)
  ;; Sets tme-db-keys to the keys used in the whole (possibly
  ;; restricted) buffer (either as entry keys or as crossref entries).
  ;; If ADD is non-nil adds the new keys to tme-db-keys instead of
  ;; simply resetting it. If VERBOSE is non-nil gives messages about
  ;; progress. If ABORTABLE is non-nil abort on user input.
  (if tme-db-maintain-sorted-entries
      (let ((case-fold-search t)
            (crossref-field
             (tme-db-cfield
              "crossref" (concat "[{\"]" tme-db-reference-key "[}\"]")))
            (labels (if add
                        tme-db-keys)))
        (save-excursion
          (goto-char (point-min))
          (if verbose
              (tme-db-progress-message
               (concat (buffer-name) ": parsing reference keys")))
          (if (catch 'userkey
                (tme-db-skip-to-valid-entry)
                (while (not (eobp))
                  (if (and
                       abortable
                       (input-pending-p))
                      (throw 'userkey t))
                  (if verbose
                      (tme-db-progress-message))
                  (let (label
                        label2)                     
                    (cond
                     ((looking-at tme-db-reference-head)
                      (setq
                       label
                       (buffer-substring-no-properties 
                        (match-beginning tme-db-key-in-head)
                        (match-end tme-db-key-in-head)))
                      (let ((p (point))
                            (m (tme-db-end-of-entry)))
                        (goto-char p)
                        (if (re-search-forward crossref-field m t)
                            (setq
                             label2
                             (buffer-substring-no-properties
                              (1+ (match-beginning tme-db-text-in-cfield))
                              (1- (match-end tme-db-text-in-cfield)))))
                        (goto-char p)))
                     ((looking-at tme-db-string)
                      (setq
                       label
                       (buffer-substring-no-properties
                        (match-beginning tme-db-key-in-string)
                        (match-end tme-db-key-in-string)))))
                    (forward-char)
                    (tme-db-skip-to-valid-entry)
                    (if (not (assoc label labels))
                        (setq labels
                              (cons (list label) labels)))
                    (if (and label2
                             (not (assoc label2 labels)))
                        (setq labels
                              (cons (list label2) labels))))))
              ;; user has aborted by typing a key --> return nil
              nil
            ;; successful operation --> return t
            (setq
             tme-db-buffer-last-parsed-tick (buffer-modified-tick)
             tme-db-keys labels)
            (if verbose
                (tme-db-progress-message 'done))
            t)))))

(defun tme-db-parse-buffers-stealthily ()
  ;; Called by tme-db-run-with-idle-timer. Whenever emacs has been idle
  ;; for tme-db-parse-keys-timeout seconds, all TME-db buffers (starting
  ;; with the current) are parsed.
  (let ((buffers (buffer-list)))
    (save-excursion
      (while (and buffers (not (input-pending-p)))
        (set-buffer (car buffers))
        (save-restriction
          (widen)
          (if (and
               (eq major-mode 'tme-db-mode)
               tme-db-maintain-sorted-entries
               (not
                (eq (buffer-modified-tick)
                    tme-db-buffer-last-parsed-tick)))
              (if (tme-db-parse-keys nil t t)
                  ;; successful operation --> remove buffer from list
                  (setq buffers (cdr buffers)))
            ;; buffer is no TME-db buffer or needs no parsing
            (setq buffers (cdr buffers))))))))

(defun tme-db-complete (string-list &optional complete-strings)
  ;; Complete word fragment before point to longest prefix of one
  ;; string defined in STRING-LIST. If point is not after the part of
  ;; a word, all strings are listed. If COMPLETE-STRINGS is non-nil,
  ;; add the strings defined in this buffer before cursor to
  ;; STRING-LIST and remove surrounding delimiters if complete string
  ;; could be expanded.
  (let* ((case-fold-search t)
         (end (point))
         (beg (save-excursion
                (re-search-backward "[ \t{\"]")
                (forward-char)
                (point)))
         (part-of-word (buffer-substring-no-properties beg end))
         (completions (copy-sequence string-list))
         (completion (save-excursion
                       (if complete-strings
                           (while (re-search-backward
                                   tme-db-string nil t)
                             (setq completions
                                   (cons
                                    (list
                                     (buffer-substring-no-properties
                                      (match-beginning tme-db-key-in-string)
                                      (match-end tme-db-key-in-string)))
                                    completions))))
                       (setq completions
                             (sort completions
                                   (lambda(x y)
                                     (string-lessp
                                      (car x)
                                      (car y)))))
                       (try-completion part-of-word completions))))
    (cond ((eq completion t)
           (if complete-strings
               ;; remove double-quotes or braces if field is no concatenation 
               (save-excursion
                 (tme-db-inside-field)
                 (tme-db-enclosing-field)
                 (let ((end (match-end tme-db-text-in-field)))
                   (goto-char (match-beginning tme-db-text-in-field))
                   (if (and
                        (looking-at tme-db-field-string)
                        (equal (match-end 0) end))
                       (tme-db-remove-delimiters))))))
          ((not completion)
           (error "Can't find completion for \"%s\"." part-of-word))
          ((not (string= part-of-word completion))
           (delete-region beg end)
           (insert completion)
           (if (and (assoc completion completions)
                    complete-strings)
               ;; remove double-quotes or braces if field is no concatenation
               (save-excursion
                 (tme-db-inside-field)
                 (tme-db-enclosing-field)
                 (let ((end (match-end tme-db-text-in-field)))
                   (goto-char (match-beginning tme-db-text-in-field))
                   (if (and
                        (looking-at tme-db-field-string)
                        (equal (match-end 0) end))
                       (tme-db-remove-delimiters))))))
          (t
           (message "Making completion list...")
           (let ((list (all-completions part-of-word completions)))
             (with-output-to-temp-buffer "*Completions*"
               (display-completion-list list)))
           (message "Making completion list...done")))))

(defun tme-db-do-auto-fill ()
  (let ((fill-prefix
         (make-string
          (+ tme-db-entry-offset tme-db-contline-indentation) ? )))
    (do-auto-fill)))

(defun tme-db-pop (arg direction)
  ;; generic function to be used by tme-db-pop-previous and tme-db-pop-next
  (let (tme-db-help-message)
    (tme-db-find-text nil))
  (save-excursion
    ;; parse current field
    (tme-db-inside-field)
    (tme-db-enclosing-field)
    (let ((case-fold-search t)
          (start-old-text (match-beginning tme-db-text-in-field))
	  (stop-old-text  (match-end tme-db-text-in-field))
	  (start-name (match-beginning tme-db-name-in-field))
	  (stop-name (match-end tme-db-name-in-field))
	  (new-text))
      (goto-char start-name)
      ;; construct regexp for field with same name as this one,
      ;; ignoring possible OPT's or ALT's
      (let ((matching-entry
	     (tme-db-cfield
	      (buffer-substring-no-properties
               (if (looking-at "OPT\\|ALT")
                   (+ (point) (length "OPT"))
                 (point))
               stop-name)
	      tme-db-field-text)))
	;; if executed several times in a row, start each search where
        ;; the last one was finished
	(cond ((eq last-command 'tme-db-pop)
               t
               )
	      (t
	       (tme-db-enclosing-reference-maybe-empty-head)
	       (setq
                tme-db-pop-previous-search-point (match-beginning 0)
                tme-db-pop-next-search-point (point))))
	(if (eq direction 'previous)
            (goto-char tme-db-pop-previous-search-point)
          (goto-char tme-db-pop-next-search-point))
        ;; Now search for arg'th previous/next similar field
	(cond
         ((if (eq direction 'previous)
              (re-search-backward matching-entry nil t arg)
            (re-search-forward matching-entry nil t arg))
          ;; Found a matching field. Remember boundaries.
	  (setq tme-db-pop-previous-search-point (match-beginning 0))
	  (setq tme-db-pop-next-search-point (match-end 0))
          (setq new-text
		(buffer-substring-no-properties
                 (match-beginning tme-db-text-in-field)
                 (match-end tme-db-text-in-field)))
          ;; change delimiters, if any changes needed
          (let ((start 0)
                old-open
                new-open
                old-close
                new-close)
            (if (equal tme-db-field-delimiters 'braces)
                (setq old-open ?\"
                      new-open ?\{
                      old-close ?\"
                      new-close ?\})
              (setq old-open ?\{
                    new-open ?\"
                    old-close ?\}
                    new-close ?\"))
            (while (string-match tme-db-field-string new-text start)
              (let ((beg (match-beginning 0))
                    (end (1- (match-end 0))))
                (if (and
                     (eq (aref new-text beg) old-open)
                     (eq (aref new-text end) old-close))
                    (progn
                      (aset new-text beg new-open)
                      (aset new-text end new-close))))
              (setq start (match-end 0))))
	  (tme-db-flash-head)
          ;; Go back to where we started, delete old text, and pop new.
	  (goto-char stop-old-text)
	  (delete-region start-old-text stop-old-text)
	  (insert new-text))
	 (t
          ;; search failed
	  (error (concat "No "
                         (if (eq direction 'previous)
                             "previous"
                           "next")
                         " matching TME-db field.")))))))
  (let (tme-db-help-message)
    (tme-db-find-text nil))
  (setq this-command 'tme-db-pop))


;; Interactive Functions:

;;;###autoload
(defun tme-db-mode ()
  "Major mode for editing TME-db files.

To submit a problem report, enter \\[tme-db-submit-bug-report] from a
TME-db mode buffer.  This automatically sets up a mail buffer with
version information already added.  You just need to add a description
of the problem, including a reproducable test case and send the
message.


General information on working with TME-db mode:

You should use commands as \\[tme-db-Book] to get a template for a
specific entry. You should then fill in all desired fields using
\\[tme-db-next-field] to jump from field to field. After having filled
in all desired fields in the entry, you should clean the new entry
with command \\[tme-db-clean-entry].

Some features of TME-db mode are available only by setting variable
tme-db-maintain-sorted-entries to t. However, then TME-db mode will
work with buffer containing only valid (syntactical correct) entries
and with entries being sorted. This is usually the case, if you have
created a buffer completely with TME-db mode and finished every new
entry with \\[tme-db-clean-entry].

For third party TME-db buffers, please call the function
`tme-db-convert-alien' to fully take advantage of all features of
TME-db mode.


Special information:

A command such as \\[tme-db-Book] will outline the fields for a TME-db book entry.

The optional fields start with the string OPT, and are thus ignored by TME-db.
Alternatives from which only one is required start with the string ALT.
The OPT or ALT string may be removed from a field with \\[tme-db-remove-OPT-or-ALT].
\\[tme-db-make-field] inserts a new field after the current one.
\\[tme-db-kill-field] kills the current field entirely.
\\[tme-db-yank] will yank the last recently killed field after the
current field.
\\[tme-db-remove-delimiters] removes the double-quotes or braces around the text of the current field.
 \\[tme-db-empty-field] replaces the text of the current field with the default \"\" or {}.

The command \\[tme-db-clean-entry] cleans the current entry, i.e. it removes OPT/ALT
from all non-empty optional or alternative fields, checks that no required
fields are empty, and does some formatting dependent on the value of
tme-db-entry-format.
Note: some functions in TME-db mode depend on entries being in a special 
format (all fields beginning on separate lines), so it is usually a bad 
idea to remove `realign' from tme-db-entry-format.

Use \\[tme-db-find-text] to position the cursor at the end of the current field.
Use \\[tme-db-next-field] to move to end of the next field.

The following may be of interest as well:

  Functions:
    tme-db-entry
    tme-db-kill-entry
    tme-db-yank-pop
    tme-db-pop-previous
    tme-db-pop-next
    tme-db-complete-string
    tme-db-complete-key
    tme-db-print-help-message
    tme-db-generate-autokey
    tme-db-beginning-of-entry
    tme-db-end-of-entry
    tme-db-reposition-window
    tme-db-mark-entry
    tme-db-ispell-abstract
    tme-db-ispell-entry
    tme-db-narrow-to-entry
    tme-db-hide-entry-bodies
    tme-db-sort-buffer
    tme-db-validate
    tme-db-count
    tme-db-fill-entry
    tme-db-reformat
    tme-db-convert-alien

  Variables:
    tme-db-field-delimiters
    tme-db-include-OPTcrossref
    tme-db-include-OPTkey
    tme-db-user-optional-fields
    tme-db-entry-format
    tme-db-sort-ignore-string-entries
    tme-db-maintain-sorted-entries
    tme-db-entry-field-alist
    tme-db-predefined-strings
    tme-db-string-files

---------------------------------------------------------
Entry to TME-db mode calls the value of `tme-db-mode-hook' if that value is
non-nil.

\\{tme-db-mode-map}
"
  (interactive)
  (kill-all-local-variables)
  (use-local-map tme-db-mode-map)
  (setq major-mode 'tme-db-mode)
  (setq mode-name "TME-db")
  (set-syntax-table tme-db-mode-syntax-table)
  (setq tme-db-strings tme-db-predefined-strings)
  (mapcar
   (lambda (filename)
     ;; collect pathnames
     (let* ((path (if tme-db-string-file-path
                      tme-db-string-file-path
                    "."))
            (dirs
             (mapcar
              (lambda (dirname)  ;; strips off trailing slashes
                (let ((len (length dirname)))
                  (if (equal (elt dirname (1- len)) "/")
                      (substring dirname 0 (1- (1- len)))
                    dirname)))
              (let (actdirs)
                (while (string-match ":" path)
                  (setq actdirs
                        (append actdirs
                                (list (substring path 0 (1- (match-end 0)))))
                        path (substring path (match-end 0))))
                (append actdirs (list path)))))
            (filename (if (string-match "\.tdb$" filename)
                          filename
                        (concat filename ".tdb")))
            fullfilename
            (item 0)
            (size (length dirs)))
       ;; test filenames
       (while (and
               (< item size)
               (not (file-readable-p
                     (setq fullfilename
                           (concat (elt dirs item) "/" filename)))))
         (setq item (1+ item)))
       (if (< item size)
           ;; file was found
           (let ((case-fold-search t)
                 (curbuf (current-buffer))
                 (bufname (make-temp-name ""))
                 (compl tme-db-strings))
             (create-file-buffer bufname)
             (set-buffer bufname)
             (insert-file-contents fullfilename)
             (goto-char (point-min))
             (while (re-search-forward tme-db-string nil t)
               (setq compl
                     (append compl
                             (list
                              (list (buffer-substring-no-properties
                                     (match-beginning tme-db-key-in-string)
                                     (match-end tme-db-key-in-string)))))))
             (kill-buffer bufname)
             (set-buffer curbuf)
             (setq tme-db-strings compl))
         (error
          "File %s not in paths defined by tme-db-string-file-path variable."
          filename))))
   tme-db-string-files)
  (if tme-db-maintain-sorted-entries
      (tme-db-run-with-idle-timer
       1 nil
       (lambda ()
         (tme-db-parse-keys nil t t))))
  ;; to get buffer parsed once if everything else (including things
  ;; installed in tme-db-mode-hook) has done its work
  (if (not tme-db-parse-idle-timer)
      (setq tme-db-parse-idle-timer
            (tme-db-run-with-idle-timer
             tme-db-parse-keys-timeout t
             'tme-db-parse-buffers-stealthily)))
  ;; Install stealthy parse function if not already installed
  (set (make-local-variable 'paragraph-start) "[ \f\n\t]*$")
  (set (make-local-variable 'comment-start) tme-db-comment-start)
  (set (make-local-variable 'comment-start-skip) tme-db-comment-start)
  (set (make-local-variable 'comment-column) 0)
  (set (make-local-variable 'normal-auto-fill-function)
       'tme-db-do-auto-fill)
  (set (make-local-variable 'font-lock-defaults)
       '(tme-db-font-lock-keywords
         nil t ((?$ . "\"")
                ;; Mathematical expressions should be fontified as strings
                (?\" . ".")
                ;; Quotes are field delimiters and quote-delimited
                ;; entries should be fontified in the same way as
                ;; brace-delimited ones
                )))
  (setq font-lock-mark-block-function
        (lambda ()
          (set-mark (tme-db-end-of-entry))
          (tme-db-beginning-of-entry)))
  (setq imenu-generic-expression
        (list (list nil tme-db-reference-head tme-db-key-in-head)))
  ;; XEmacs needs easy-menu-add, Emacs does not care
  (easy-menu-add tme-db-edit-menu)
  (easy-menu-add tme-db-entry-menu)
  (run-hooks 'tme-db-mode-hook))

(defun tme-db-submit-bug-report ()
  "Submit via mail a bug report on tme-db.el."
  (interactive)
  (if (y-or-n-p "Do you want to submit a bug report on TME-db mode? ")
      (progn
        (require 'reporter)
        (let ((reporter-prompt-for-summary-p t))
          (reporter-submit-bug-report
           tme-db-maintainer-address
           (concat "tme-db.el " "(emacs 19.35)")
           (list
            'system-configuration
            'system-configuration-options
            'tme-db-mode-hook
            'tme-db-parse-keys-timeout
            ;; possible general errors
            'tme-db-sort-ignore-string-entries
            'tme-db-maintain-sorted-entries
            'tme-db-entry-delimiters
            'tme-db-field-delimiters
            'tme-db-comma-after-last-field
            'tme-db-entry-offset
            'tme-db-field-indentation
            'tme-db-text-indentation
            'tme-db-contline-indentation
            'tme-db-align-at-equal-sign
            ;; possible sorting and parsing bugs
            'tme-db-entry-format
            'tme-db-add-entry-hook
            'tme-db-clean-entry-hook
            ;; possible cleaning error
            'tme-db-user-optional-fields
            ;; possible format error
            'tme-db-predefined-month-strings
            'tme-db-predefined-strings
            'tme-db-string-files
            'tme-db-string-file-path
            ;; possible format error
            'tme-db-font-lock-keywords
            ;; possible bugs regarding fontlocking
            'tme-db-autokey-prefix-string
            'tme-db-autokey-names
            'tme-db-autokey-names-stretch
            'tme-db-autokey-additional-names
            'tme-db-autokey-transcriptions
            'tme-db-autokey-name-change-strings
            'tme-db-autokey-name-case-convert
            'tme-db-autokey-name-length
            'tme-db-autokey-name-separator
            'tme-db-autokey-year-length
            'tme-db-autokey-year-use-crossref-entry
            'tme-db-autokey-titlewords
            'tme-db-autokey-title-terminators
            'tme-db-autokey-titlewords-stretch
            'tme-db-autokey-titleword-ignore
            'tme-db-autokey-titleword-case-convert
            'tme-db-autokey-titleword-abbrevs
            'tme-db-autokey-titleword-change-strings
            'tme-db-autokey-titleword-length
            'tme-db-autokey-titleword-separator
            'tme-db-autokey-name-year-separator
            'tme-db-autokey-year-title-separator
            'tme-db-autokey-edit-before-use
            'tme-db-autokey-before-presentation-function
            ;; possible bugs regarding automatic labels
            'tme-db-entry-field-alist
            ;; possible format error
            'tme-db-help-message
            'tme-db-include-OPTcrossref
            'tme-db-include-OPTkey
            'tme-db-field-kill-ring-max
            'tme-db-entry-kill-ring-max
            ;; user variables which shouldn't cause any errors
            )
           nil nil
           (concat "Hallo Dirk,
 
I want to report a bug on Emacs TME-db mode.
I've read the `Bugs' section in the `Emacs' info page, so I know how
to make a clear and unambiguous report. I have started a fresh Emacs
via `"invocation-name " --no-init-file --no-site-file', thereafter (in
case I'm reporting on a version of `tme-db.el' which is not part of
the standard emacs distribution) I loaded the questionable version
of `tme-db.el' with `M-x load-file', and then, to produce the buggy
behaviour, I did the following:")))
        (message nil))))

(defun tme-db-entry (entry-type)
  "Inserts a new TME-db entry.
After insertion it calls the functions in `tme-db-add-entry-hook'."
  (interactive (let* ((completion-ignore-case t)
		      (e-t (completing-read
                            "Entry Type: "
                            tme-db-entry-field-alist
                            nil t nil 'tme-db-entry-type-history)))
		 (list e-t)))
  (if (not tme-db-keys)
      (tme-db-parse-keys nil t))
  (let* (required
         optional
         (key
          (if tme-db-maintain-sorted-entries
              (completing-read
               (format "%s key: " entry-type)
               tme-db-keys nil nil nil 'tme-db-key-history)))
         (e (assoc-ignore-case entry-type tme-db-entry-field-alist))
         (r-n-o (elt e 1))
         (c-ref (elt e 2)))
    (if (not e)
        (error "TME-db entry type %s not defined." entry-type))
    (if (and
         (member entry-type tme-db-include-OPTcrossref)
         c-ref)
        (setq required (elt c-ref 0)
              optional (elt c-ref 1))
      (setq required (elt r-n-o 0)
            optional (elt r-n-o 1)))
    (if tme-db-maintain-sorted-entries
	(tme-db-find-entry-location key)
      (tme-db-move-outside-of-entry))
    (indent-to-column tme-db-entry-offset)
    (insert "@" entry-type (tme-db-entry-left-delimiter))
    (if key
	(insert key))
    (save-excursion
      (mapcar 'tme-db-make-field required)
      (if (member entry-type tme-db-include-OPTcrossref)
	  (tme-db-make-optional-field '("crossref")))
      (if tme-db-include-OPTkey
          (if (or
               (stringp tme-db-include-OPTkey)
               (fboundp tme-db-include-OPTkey))
              (tme-db-make-optional-field
               (list "key" nil tme-db-include-OPTkey))
            (tme-db-make-optional-field '("key"))))            
      (mapcar 'tme-db-make-optional-field optional)
      (mapcar 'tme-db-make-optional-field tme-db-user-optional-fields)
      (if tme-db-comma-after-last-field
          (insert ","))
      (insert "\n")
      (indent-to-column tme-db-entry-offset)
      (insert (tme-db-entry-right-delimiter) "\n\n"))
    (tme-db-next-field t)
    (run-hooks 'tme-db-add-entry-hook)))

(defun tme-db-print-help-message ()
  "Prints helpful information about current field in current TME-db entry."
  (interactive)
  (let* ((case-fold-search t)
         (pnt (point))
         (field-name
          (progn
            (condition-case errname
                (tme-db-enclosing-field)
              (search-failed
               (goto-char pnt)
               (error "Not on TME-db field.")))
            (let ((mb (match-beginning tme-db-name-in-field))
                  (me (match-end tme-db-name-in-field)))
              (goto-char mb)
              (buffer-substring-no-properties
               (if (looking-at "OPT\\|ALT")
                   (+ 3 mb)
                 mb)
               me))))
         (reference-type
          (progn
            (re-search-backward
             tme-db-reference-maybe-empty-head nil t)
            (buffer-substring-no-properties
             (1+ (match-beginning tme-db-type-in-head))
             (match-end tme-db-type-in-head))))
         (entry-list
          (assoc-ignore-case reference-type
                               tme-db-entry-field-alist))
         (c-r-list (elt entry-list 2))
         (req-opt-list
          (if (and
               (member reference-type tme-db-include-OPTcrossref)
               c-r-list)
              c-r-list
            (elt entry-list 1)))
         (list-of-entries (append
                           (elt req-opt-list 0)
                           (elt req-opt-list 1)
                           tme-db-user-optional-fields
                           (if (member
                                reference-type
                                tme-db-include-OPTcrossref)
                               '(("crossref"
                                  "Label of the crossreferenced entry")))
                           (if tme-db-include-OPTkey
                               '(("key"
                                  "Key used for label creation if author and editor fields are missing"))))))
    (goto-char pnt)
    (let ((comment (assoc-ignore-case field-name list-of-entries)))
      (if comment
          (message (elt comment 1))
        (message "NO COMMENT AVAILABLE")))))

(defun tme-db-make-field (e-t &optional called-by-yank)
  "Makes a field named E-T in current TME-db entry.
This function is for interactive and non-interactive purposes. To call
it interactively, just give it no arguments and enter the field name
using the minibuffer."
  (interactive "*P")
  (if (not e-t)
      (setq
       e-t
       (let* ((reference-type
               (save-excursion
                 (tme-db-enclosing-reference-maybe-empty-head)
                 (buffer-substring-no-properties
                  (1+ (match-beginning tme-db-type-in-head))
                  (match-end tme-db-type-in-head))))
              (fl
               (car (cdr (assoc-ignore-case
                          reference-type tme-db-entry-field-alist))))
              (field-list
               (append
                (elt fl 0) (elt fl 1) tme-db-user-optional-fields
                (if tme-db-include-OPTcrossref '(("crossref" nil)))
                (if tme-db-include-OPTkey '(("key" nil)))))
              (completion-ignore-case t))
         (completing-read
          "TME-db field name: " field-list
          nil nil nil tme-db-field-history))))
  (if (not (consp e-t))
      (setq e-t (list e-t)))
  (if (equal (length e-t) 1)
      (setq e-t (append e-t (list ""))))
  (if (equal (length e-t) 2)
      (setq e-t (append e-t (list ""))))
  (let ((name (if (elt e-t 3)
                  (concat "ALT" (car e-t))
                (car e-t))))
    (if (or (interactive-p) called-by-yank)
        (let (tme-db-help-message)
          (tme-db-find-text nil t)
          (if (looking-at "[}\"]")
              (forward-char))))
    (insert ",\n")
    (indent-to-column
     (+ tme-db-entry-offset tme-db-field-indentation))
    (insert name " ")
    (if tme-db-align-at-equal-sign
        (indent-to-column
         (+ tme-db-entry-offset (- tme-db-text-indentation 2))))
    (insert "= ")
    (if (not tme-db-align-at-equal-sign)
        (indent-to-column
         (+ tme-db-entry-offset tme-db-text-indentation)))
    (insert (if called-by-yank
                ""
              (tme-db-field-left-delimiter))
            (let ((init (elt e-t 2)))
              (cond
               ((stringp init)
                init)
               ((fboundp init)
                (funcall init))
               (t
                (error "%s is neither a string nor a function." init))))
            (if called-by-yank
                ""
              (tme-db-field-right-delimiter)))
    (if (interactive-p)
        (forward-char -1))))

(defun tme-db-beginning-of-entry ()
  "Move to beginning of TME-db entry.
If inside an entry, move to the beginning of it, otherwise move to the
beginning of the previous entry.
If called from a program, this function returns the new location of point."
  (interactive)
  (skip-chars-forward " \t")
  (if (looking-at "@")
      (forward-char))
  (re-search-backward "^[ \t]*@" nil 'move))

(defun tme-db-end-of-entry ()
  "Move to end of TME-db entry.
If inside an entry, move to the end of it, otherwise move to the end
of the previous entry.
If called from a program, this function returns the new location of point."
  (interactive)
  (let ((case-fold-search t)
        (valid-entry-head 
         (concat "[ \t]*@[ \t]*\\("
                 (mapconcat
                  (lambda (type)
                    (concat "\\(" (car type) "\\)"))
                  tme-db-entry-field-alist
                  "\\|")
                 "\\)"))
        (org (point))
        (pnt (tme-db-beginning-of-entry))
        err)
    (cond
     ((looking-at "[ \t]*@[ \t]*string[ \t\n]*[({]")
      (if (not (and
                (re-search-forward tme-db-string nil t)
                (equal (match-beginning 0) pnt)))
          (setq err t)))
     ((looking-at "[ \t]*@[ \t]*preamble[ \t\n]*")
      (goto-char (match-end 0))
      (cond
       ((looking-at "(")
        (if (not (re-search-forward ")[ \t]*\n\n" nil 'move))
            (setq err t)))
       ((looking-at "{")
        (if (not (re-search-forward "}[ \t]*\n\n" nil 'move))
            (setq err t)))
       (t
        (setq err t)))
      (if (not err)
          (progn
            (goto-char (match-beginning 0))
            (forward-char))))
     ((looking-at valid-entry-head)
      (tme-db-search-reference t nil t)
      (if (not (equal (match-beginning 0) pnt))
          (setq err t)))
     (t
      (if (interactive-p)
          (message "Not on a known TME-db entry."))
      (goto-char org)))
    (if err
        (progn
          (goto-char pnt)
          (error "Syntactical incorrect entry starts here."))))
  (point))
  
(defun tme-db-reposition-window (arg)
  "Make the current TME-db entry visible."
  (interactive "P")
  (save-excursion
    (goto-char
     (/ (+ (tme-db-beginning-of-entry) (tme-db-end-of-entry)) 2))
    (recenter arg)))

(defun tme-db-mark-entry ()
  "Put mark at beginning, point at end of current TME-db entry."
  (interactive)
  (set-mark (tme-db-beginning-of-entry))
  (tme-db-end-of-entry))

(defun tme-db-count-entries (&optional count-string-entries)
  "Count number of entries in current buffer or region.
With prefix argument it counts all entries, otherwise it counts all
except Strings.
If mark is active it counts entries in region, if not in whole buffer."
  (interactive "P")
  (let ((pnt (point))
        (start-point
         (if (tme-db-mark-active)
             (region-beginning)
           (tme-db-beginning-of-first-entry)))
        (end-point
         (if (tme-db-mark-active)
             (region-end)
           (point-max)))
        (number 0)
        (tme-db-sort-ignore-string-entries
         (not count-string-entries)))
    (save-restriction
      (narrow-to-region start-point end-point)
      (goto-char start-point)
      (tme-db-map-entries
       (lambda (current)
         (setq number (1+ number)))))
    (message (concat (if (tme-db-mark-active) "Region" "Buffer")
                     " contains %d entries.") number)
    (goto-char pnt)))

(defun tme-db-ispell-entry ()
  "Spell whole TME-db entry."
  (interactive)
  (ispell-region (tme-db-beginning-of-entry) (tme-db-end-of-entry)))

(defun tme-db-ispell-abstract ()
  "Spell abstract of TME-db entry."
  (interactive)
  (let ((case-fold-search t)
        (pnt (tme-db-end-of-entry)))
    (tme-db-beginning-of-entry)
    (if (not
         (re-search-forward
          (tme-db-cfield "abstract" tme-db-field-text) pnt t))
        (error "No abstract in entry.")))
  (ispell-region (match-beginning tme-db-text-in-cfield)
                 (match-end tme-db-text-in-cfield)))

(defun tme-db-narrow-to-entry ()
  "Narrow buffer to current TME-db entry."
  (interactive)
  (save-excursion
    (narrow-to-region
     (tme-db-beginning-of-entry) (tme-db-end-of-entry))))

(defun tme-db-hide-entry-bodies (&optional arg)
  "Hide all lines between first and last TME-db entries not beginning with @.
With argument, show all text."
  (interactive "P")
  (save-excursion
    (tme-db-beginning-of-first-entry)
    (let ((buffer-read-only nil))
      (if arg
	  (subst-char-in-region (point) (point-max) ?\r ?\n t)
        (while (not (eobp))
          (subst-char-in-region
           (point)
           (progn
             (re-search-forward "[\n\r]@" nil t)
             (forward-line -1)
             (point))
           ?\n ?\r t)
          (forward-line 1)))
      (setq selective-display (not arg)))))

(defun tme-db-sort-buffer ()
  "Sort TME-db buffer alphabetically by key.
Text outside of TME-db entries is not affected.  If
`tme-db-sort-ignore-string-entries' is non-nil, @String entries will be
ignored."
  (interactive)
  (save-restriction
    (narrow-to-region
     (tme-db-beginning-of-first-entry)
     (save-excursion
       (goto-char (point-max))
       (tme-db-end-of-entry)))
    (tme-db-skip-to-valid-entry)
    (sort-subr
     nil
     ;; NEXTREC function
     'tme-db-skip-to-valid-entry
     ;; ENDREC function
     'tme-db-end-of-entry
     ;; STARTKEY function
     (lambda ()
       (let ((case-fold-search t))
         (re-search-forward tme-db-reference-head)
         (buffer-substring-no-properties
          (match-beginning tme-db-key-in-head)
          (match-end tme-db-key-in-head)))))))
  
(defun tme-db-find-entry-location (entry-name &optional ignore-dups)
  "Looking for place to put the TME-db entry named ENTRY-NAME.
Performs a binary search (therefore, buffer is assumed to be in sorted
order, without duplicates (see \\[tme-db-validate]), if it is
not, tme-db-find-entry-location will fail). If entry-name is already
used as a reference key, an error is signaled. However, if optional
variable IGNORE-DUPS is non-nil, no error messages about duplicate
entries are signaled, but the error handling is assumed to be made in
the calling function.
The value is nil if an duplicate entry error occurred,
and t in all other cases."
  (let* ((case-fold-search t)
         (left
          (progn
            (tme-db-beginning-of-first-entry)
            (tme-db-skip-to-valid-entry)
            (tme-db-end-of-entry)))
         (right
          (progn
            (tme-db-beginning-of-last-entry)
            (tme-db-skip-to-valid-entry t)
            (point)))
         actual-point
         actual-key
         (done (>= left right))
         new
         dup)
    (while (not done)
      (setq actual-point (/ (+ left right) 2))
      (goto-char actual-point)
      (tme-db-skip-to-valid-entry t)
      (setq actual-key
            (progn
              (re-search-forward tme-db-reference-head)
              (buffer-substring-no-properties
               (match-beginning tme-db-key-in-head)
               (match-end tme-db-key-in-head))))
      (cond
       ((string-lessp entry-name actual-key)
        (setq new (tme-db-beginning-of-entry))
        (if (equal right new)
            (setq done t)
          (setq right new)))
       ((string-lessp actual-key entry-name)
        (setq new (tme-db-end-of-entry))
        (if (equal left new)
            (setq done t)
          (setq left new)))
       ((string-equal actual-key entry-name)
        (setq dup t
              done t)
        (if (not ignore-dups)
            (progn
              (tme-db-beginning-of-entry)
              (error "Entry with key `%s' already exists." entry-name))))))
    (if dup
        (progn
          (tme-db-beginning-of-entry)
          nil)
      (goto-char right)
      (setq actual-key
            (if (looking-at tme-db-reference-head)
                (buffer-substring-no-properties
                 (match-beginning tme-db-key-in-reference)
                 (match-end tme-db-key-in-reference))))
      (if (or
           (not actual-key) 
           (string-lessp actual-key entry-name)) 
          ;; buffer contains no valid entries or
          ;; greater than last entry --> append
          (progn
            (tme-db-end-of-entry)
            (if (not (bobp))
                (newline (forward-line 2)))
            (beginning-of-line))
        (goto-char right))
      t)))    

(defun tme-db-validate (&optional test-thoroughly)
  "Validate if buffer or region is syntactically correct.
Only known reference types are checked, so you can put comments
outside of entries.
With optional argument TEST-THOROUGHLY non-nil it checks for absence of
required fields and questionable month fields as well.
If mark is active, it validates current region, if not whole buffer.
Returns t if test was successful, nil otherwise."
  (interactive "P")
  (let (error-list
        syntax-error
        (case-fold-search t)
        (valid-tme-db-entry
         (concat
          "@[ \t]*\\(\\(string\\)\\|"
          (mapconcat
           (lambda (type)
             (concat "\\(" (car type) "\\)"))
           tme-db-entry-field-alist
           "\\|")
          "\\)"))
        (pnt (point))
        (start-point
         (if (tme-db-mark-active)
             (region-beginning)
           (tme-db-beginning-of-first-entry)))
        (end-point
         (if (tme-db-mark-active)
             (region-end)
           (point-max))))
    (save-restriction
      (narrow-to-region start-point end-point)
      ;; looking if entries fit syntactical structure
      (goto-char start-point)
      (tme-db-progress-message "Checking syntactical structure")
      (while (re-search-forward "^[ \t]*@" nil t)
        (tme-db-progress-message)
        (forward-char -1)
        (let ((p (point))
              (must-match
               (looking-at valid-tme-db-entry)))
          (if (not must-match)
              (forward-char)
            (let (tme-db-sort-ignore-string-entries)
              (tme-db-skip-to-valid-entry))
            (if (equal (point) p)
                (forward-char)
              (goto-char p)
              (setq
               error-list
               (cons (list
                      (tme-db-current-line)
                      "Syntax error (check esp. commas, braces, and quotes)") 
                     error-list))
              (forward-char)))))
      (tme-db-progress-message 'done)
      (if error-list
          (setq syntax-error t)
        ;; looking for correct sort order and duplicates (only if
        ;; there were no syntax errors)
        (if tme-db-maintain-sorted-entries
            (let (previous)
              (goto-char start-point)
              (tme-db-progress-message "Checking correct sort order")
              (tme-db-map-entries
               (lambda (current)
                 (tme-db-progress-message)
                 (cond ((or (not previous)
                            (string< previous current))
                        (setq previous current))
                       ((string-equal previous current)
                        (setq
                         error-list
                         (cons (list (tme-db-current-line)
                                     "Duplicate key with previous")
                               error-list)))
                       (t
                        (setq previous current
                              error-list
                              (cons (list (tme-db-current-line)
                                          "Entries out of order")
                                    error-list))))))
              (tme-db-progress-message 'done)))
        (if test-thoroughly
            (progn
              (goto-char start-point)
              (tme-db-progress-message
               "Checking required fields and month fields")
              (let ((tme-db-sort-ignore-string-entries t)
                    (questionable-month
                     (concat
                      "[{\"]\\("
                      (mapconcat
                       (lambda (mon)
                         (concat "\\(" (car mon) "\\)"))
                       tme-db-predefined-month-strings
                       "\\|")
                      "\\)[}\"]")))
                (tme-db-map-entries
                 (lambda (current)
                   (tme-db-progress-message)
                   (let* ((beg (tme-db-beginning-of-entry))
                          (end (tme-db-end-of-entry))
                          (entry-list
                           (progn
                             (goto-char beg)
                             (tme-db-search-reference nil end)
                             (assoc-ignore-case
                              (buffer-substring-no-properties
                               (1+ (match-beginning tme-db-type-in-head))
                               (match-end tme-db-type-in-head))
                              tme-db-entry-field-alist)))
                          (req (copy-sequence (elt (elt entry-list 1) 0)))
                          (creq (copy-sequence (elt (elt entry-list 2) 0)))
                          crossref-there)
                     (goto-char beg)
                     (while (re-search-forward tme-db-field end t)
                       (let ((field-name
                              (buffer-substring-no-properties
                               (match-beginning tme-db-name-in-field)
                               (match-end tme-db-name-in-field))))
                         (if (and (equal (downcase field-name) "month")
                                  (string-match
                                   questionable-month
                                   (buffer-substring-no-properties
                                    (match-beginning tme-db-text-in-field)
                                    (match-end tme-db-text-in-field))))
                             (setq
                              error-list
                              (cons
                               (list
                                (tme-db-current-line)
                                "Questionable month field (delimited string)")
                               error-list)))
                         (setq
                          req
                          (delete (assoc-ignore-case field-name req) req)
                          creq
                          (delete (assoc-ignore-case field-name creq) creq))
                         (if (equal (downcase field-name) "crossref")
                             (setq crossref-there t))))
                     (if crossref-there
                         (setq req creq))
                     (if (or (> (length req) 1)
                             (and (= (length req) 1)
                                  (not (elt (car req) 3))))
                         ;; two (or more) fields missed or one field
                         ;; missed and this isn't flagged alternative
                         ;; (notice that this fails if there are more
                         ;; than two alternatives in a TME-db entry,
                         ;; which isn't the case momentarily)
                         (setq
                          error-list
                          (cons
                           (list (save-excursion
                                   (tme-db-beginning-of-entry)
                                   (tme-db-current-line))
                                 (concat
                                  "Required field \""
                                  (car (car req))
                                  "\" missing"))
                           error-list)))))))
              (tme-db-progress-message 'done)))))
    (goto-char pnt)
    (if error-list
        (let ((bufnam (buffer-name))
              (dir default-directory))
          (setq error-list
                (sort error-list
                      (lambda (a b)
                        (< (car a) (car b)))))
          (let ((pop-up-windows t))
            (pop-to-buffer nil t))
          (switch-to-buffer
           (get-buffer-create "*TME-db validation errors*") t)
          ;; don't use switch-to-buffer-other-window, since this
          ;; doesn't allow the second parameter NORECORD
          (setq default-directory dir)
          (toggle-read-only -1)
          (compilation-mode)
          (delete-region (point-min) (point-max))
          (goto-char (point-min))
          (insert
           "TME-db mode command `tme-db-validate'\n"
           (if syntax-error
               "Maybe undetected errors due to syntax errors. Correct and validate again."
             "")
           "\n")
          (while error-list
            (insert
             bufnam ":" (number-to-string (elt (car error-list) 0))
             ": " (elt (car error-list) 1) "\n")
            (setq error-list (cdr error-list)))
          (compilation-parse-errors nil nil)
          (setq compilation-old-error-list compilation-error-list)
          ;; this is necessary to avoid reparsing of buffer if you
          ;; switch to compilation buffer and enter
          ;; `compile-goto-error'
          (set-buffer-modified-p nil)
          (toggle-read-only 1)
          (goto-char (point-min))
          (other-window -1)
          ;; return nil
          nil)
      (if (tme-db-mark-active)
          (message "Region is syntactically correct")
        (message "Buffer is syntactically correct"))
      t)))

(defun tme-db-next-field (arg)
  "Finds end of text of next TME-db field; with arg, to its beginning."
  (interactive "P")
  (tme-db-inside-field)
  (let ((start (point)))
    (condition-case ()
	(progn
	  (tme-db-enclosing-field)
	  (goto-char (match-end 0))
	  (forward-char 2))
      (error
       (goto-char start)
       (end-of-line)
       (forward-char))))
  (tme-db-find-text arg t))

(defun tme-db-find-text (arg &optional as-if-interactive silent)
  "Go to end of text of current field; with ARG, go to beginning."
  (interactive "P")
  (tme-db-inside-field)
  (if (tme-db-enclosing-field (or (interactive-p) as-if-interactive))
      (progn
        (if arg
            (progn
              (goto-char (match-beginning tme-db-text-in-field))
              (if (looking-at "[{\"]")
                  (forward-char)))
          (goto-char (match-end tme-db-text-in-field))
          (if (or
               (= (preceding-char) ?})
               (= (preceding-char) ?\"))
              (forward-char -1)))
        (if tme-db-help-message
            (tme-db-print-help-message)))
    (beginning-of-line)
    (cond
     ((looking-at tme-db-string)
      (goto-char
       (if arg
           (match-beginning tme-db-text-in-string)
         (match-end tme-db-text-in-string))))
     ((looking-at tme-db-reference-maybe-empty-head)
      (goto-char
       (if arg
           (match-beginning tme-db-key-in-head)
         (match-end 0))))
     (t
      (if (not silent)
          (error "Not on TME-db field."))))))

(defun tme-db-remove-OPT-or-ALT ()
  "Removes the string starting optional/alternative fields.
Aligns text and goes thereafter to end of text."
  (interactive)
  (let ((case-fold-search t))
    (tme-db-inside-field)
    (tme-db-enclosing-field)
    (save-excursion
      (goto-char (match-beginning tme-db-name-in-field))
      (if (looking-at "OPT\\|ALT")
          (progn
            (delete-char (length "OPT"))
            ;; make field non-OPT
            (search-forward "=")
            (forward-char -1)
            (delete-horizontal-space)
            (if tme-db-align-at-equal-sign
                (indent-to-column (- tme-db-text-indentation 2))
              (insert " "))
            (search-forward "=")
            (delete-horizontal-space)
            (if tme-db-align-at-equal-sign
                (insert " ")
              (indent-to-column tme-db-text-indentation)))))
    (tme-db-inside-field)))

(defun tme-db-remove-delimiters ()
  "Removes \"\" or {} around string."
  (interactive)
  (let ((case-fold-search t))
    (save-excursion
      (tme-db-inside-field)
      (tme-db-enclosing-field)
      (let ((start (match-beginning tme-db-text-in-field))
            (stop (match-end tme-db-text-in-field)))
        (goto-char start)
        (while (re-search-forward tme-db-field-string stop t)
          (let ((beg (copy-marker (match-beginning 0)))
                (end (copy-marker (match-end 0))))
            (goto-char beg)
            (if (looking-at "[{\"]")
                (delete-char 1))
            (goto-char end)
            (forward-char -1)
            (if (looking-at "[}\"]")
                (delete-char 1))))))))

(defun tme-db-kill-field (&optional copy-only)
  "Kills the entire enclosing TME-db field.
With prefix arg, copy the current field to `tme-db-field-kill-ring,'
but do not actually kill it."
  (interactive "P")
  (let ((pnt (point))
        (case-fold-search t))
    (tme-db-inside-field)
    (tme-db-enclosing-field)
    (let ((the-end (match-end 0))
          (the-beginning (match-beginning 0)))
      (goto-char the-end)
      (skip-chars-forward " \t\n,")
      (setq
       tme-db-field-kill-ring
       (cons
        (list
         'field
         (buffer-substring-no-properties
          (match-beginning tme-db-name-in-field)
          (match-end tme-db-name-in-field))
         (buffer-substring-no-properties
          (match-beginning tme-db-text-in-field)
          (match-end tme-db-text-in-field)))
        tme-db-field-kill-ring))
      (if (> (length tme-db-field-kill-ring) tme-db-field-kill-ring-max)
          (setcdr
           (nthcdr (1- tme-db-field-kill-ring-max) tme-db-field-kill-ring)
           nil))
      (setq tme-db-field-kill-ring-yank-pointer tme-db-field-kill-ring)
      (if copy-only
          (goto-char pnt)
        (delete-region the-beginning the-end)
        (let (tme-db-help-message)
          (tme-db-find-text nil t t)))))
  (setq tme-db-last-kill-command 'field))

(defun tme-db-copy-field-as-kill ()
  (interactive)
  (tme-db-kill-field t))

(defun tme-db-kill-entry (&optional copy-only)
  "Kill the entire enclosing TME-db reference entry.
With prefix arg copy the current reference entry to
`tme-db-entry-kill-ring', but do not actually kill it."
  (interactive "P")
  (let ((pnt (point))
        (case-fold-search t)
        (beg (tme-db-beginning-of-entry))
        (end
         (progn
           (tme-db-end-of-entry)
           (if (re-search-forward
                tme-db-reference-maybe-empty-head nil 'move)
               (goto-char (match-beginning 0)))
           (point))))
    (setq
     tme-db-entry-kill-ring
     (cons
      (list 'entry (buffer-substring-no-properties beg end))
      tme-db-entry-kill-ring))
    (if (> (length tme-db-entry-kill-ring) tme-db-entry-kill-ring-max)
        (setcdr
         (nthcdr (1- tme-db-entry-kill-ring-max) tme-db-entry-kill-ring)
         nil))
    (setq tme-db-entry-kill-ring-yank-pointer tme-db-entry-kill-ring)
    (if copy-only
        (goto-char pnt)
      (delete-region beg end)))
  (setq tme-db-last-kill-command 'entry))

(defun tme-db-copy-entry-as-kill ()
  (interactive)
  (tme-db-kill-entry t))

(defun tme-db-yank (&optional n)
  "Reinsert the last TME-db item.
More precisely, reinsert the field or entry killed or yanked most recently.
With argument N, reinsert the Nth most recently killed TME-db item.
See also the command \\[tme-db-yank-pop]]."
  (interactive "*p")
  (tme-db-insert-current-kill (1- n))
  (setq this-command 'tme-db-yank))    

(defun tme-db-yank-pop (n)
  "Replace just-yanked killed TME-db item with a different.
This command is allowed only immediately after a `tme-db-yank' or a
`tme-db-yank-pop'.
At such a time, the region contains a reinserted previously killed
TME-db item. `tme-db-yank-pop' deletes that item and inserts in its
place a different killed TME-db item.

With no argument, the previous kill is inserted.
With argument N, insert the Nth previous kill.
If N is negative, this is a more recent kill.

The sequence of kills wraps around, so that after the oldest one
comes the newest one."
  (interactive "*p")
  (if (not (eq last-command 'tme-db-yank))
      (error "Previous command was not a TME-db yank"))
  (setq this-command 'tme-db-yank)
  (let ((inhibit-read-only t))
    (delete-region (point) (mark t))
    (tme-db-insert-current-kill n)))

(defun tme-db-empty-field ()
  "Delete the text part of the current field, replace with empty text."
  (interactive)
  (tme-db-inside-field)
  (tme-db-enclosing-field)
  (goto-char (match-beginning tme-db-text-in-field))
  (delete-region (point) (match-end tme-db-text-in-field))
  (insert (concat (tme-db-field-left-delimiter)
                  (tme-db-field-right-delimiter)) )
  (tme-db-find-text t))

(defun tme-db-pop-previous (arg)
  "Replace text of current field with the similar field in previous entry.
With arg, goes up ARG entries. Repeated, goes up so many times. May be
intermixed with \\[tme-db-pop-next] (tme-db-pop-next)."
  (interactive "p")
  (tme-db-pop arg 'previous))

(defun tme-db-pop-next (arg)
  "Replace text of current field with the text of similar field in next entry.
With arg, goes down ARG entries. Repeated, goes down so many times. May be
intermixed with \\[tme-db-pop-previous] (tme-db-pop-previous)."
  (interactive "p")
  (tme-db-pop arg 'next))

(defun tme-db-clean-entry (&optional new-label called-by-reformat)
  "Finish editing the current TME-db entry and clean it up.
Checks that no required fields are empty and formats entry dependent
on the value of tme-db-entry-format.
If label of entry is empty or a prefix argument is given, calculate a
new entry label (note: this only will work if fields in entry begin on
separate lines prior to calling tme-db-clean-entry or if 'realign is
contained in tme-db-entry-format).
Don't call this on `string' or `preamble' entries.
At end of the cleaning process, the functions in
tme-db-clean-entry-hook are called with region narrowed to entry."
  (interactive "P")
  (tme-db-format-entry)
  (let* ((case-fold-search t)
         (eob (tme-db-end-of-entry))
         (key (progn
                (tme-db-beginning-of-entry)
                (if (re-search-forward
                     tme-db-reference-head eob t)
                    (buffer-substring-no-properties
                     (match-beginning tme-db-key-in-head)
                     (match-end tme-db-key-in-head))))))
    (if (or
         new-label
         (not key))
        (progn
          (let ((autokey
                 (if tme-db-autokey-edit-before-use
                     (read-from-minibuffer
                      "Key to use: " (tme-db-generate-autokey) nil nil
                      'tme-db-key-history)
                   (tme-db-generate-autokey))))
            (tme-db-beginning-of-entry)
            (re-search-forward tme-db-reference-maybe-empty-head)
            (if (match-beginning tme-db-key-in-head)
                (delete-region (match-beginning tme-db-key-in-head)
                               (match-end tme-db-key-in-head)))
            (insert autokey)
            (let* ((start (tme-db-beginning-of-entry))
                   (end (progn
                          (tme-db-end-of-entry)
                          (if (re-search-forward
                               tme-db-reference-maybe-empty-head nil 'move)
                              (goto-char (match-beginning 0)))
                          (point)))
                   (entry (buffer-substring start end)))
              (delete-region start end)
              (let ((success
                     (or
                      called-by-reformat
                      (not tme-db-maintain-sorted-entries)
                      (tme-db-find-entry-location autokey t))))
                (insert entry)
                (forward-char -1)
                (tme-db-beginning-of-entry)
                (re-search-forward tme-db-reference-head)
                (if (not success)
                    (error
                     "New inserted reference yields duplicate key."))))))))
  (if (not called-by-reformat)
      (save-excursion
        (save-restriction
          (narrow-to-region
           (tme-db-beginning-of-entry) (tme-db-end-of-entry))
          (tme-db-parse-keys t nil)
          (run-hooks 'tme-db-clean-entry-hook)))))

(defun tme-db-fill-entry ()
  "Fill current entry.
Realigns entry, so that every field starts on a separate line. Field
names appear in column `tme-db-field-indentation', field text starts in
column tme-db-text-indentation and continuation lines start here, too.
If `tme-db-align-at-equal-sign' is non-nil, align equal signs also."
  (interactive "*")
  (let ((pnt (copy-marker (point)))
        (end (copy-marker (tme-db-end-of-entry))))
    (tme-db-beginning-of-entry)
    (tme-db-delete-whitespace)
    (indent-to-column tme-db-entry-offset)
    (while (re-search-forward tme-db-field end t)
      (let* ((begin-field
              (copy-marker (match-beginning 0)))
             (end-field
              (copy-marker (match-end 0)))
             (begin-name
              (copy-marker (match-beginning tme-db-name-in-field)))
             (end-name
              (copy-marker (match-end tme-db-name-in-field))))
        (goto-char begin-field)
        (forward-char)
        (tme-db-delete-whitespace)
        (open-line 1)
        (forward-char)
        (indent-to-column
         (+ tme-db-entry-offset tme-db-field-indentation))
        (re-search-forward "[ \t\n]*=" end)
        (replace-match "=")
        (forward-char -1)
        (if tme-db-align-at-equal-sign
            (indent-to-column
             (+ tme-db-entry-offset (- tme-db-text-indentation 2)))
          (insert " "))
        (forward-char)
        (tme-db-delete-whitespace)
        (if tme-db-align-at-equal-sign
            (insert " ")
          (indent-to-column tme-db-text-indentation))
        (while (re-search-forward "[ \t\n]+" end-field 'move)
          (replace-match " "))
        (tme-db-do-auto-fill)))
    (if (looking-at ",")
        (forward-char))
    (tme-db-delete-whitespace)
    (open-line 1)
    (forward-char)
    (indent-to-column tme-db-entry-offset)
    (goto-char pnt)))

(defun tme-db-reformat (&optional additional-options called-by-convert-alien)
  "Reformat all TME-db entries in buffer or region.
With prefix argument, read options for reformatting from minibuffer.
With C-u C-u prefix argument, reuse previous answers (if any) again.
If mark is active it reformats entries in region, if not in whole buffer."
  (interactive "*P")
  (let* ((pnt (point))
         (use-previous-options
          (and (equal (prefix-numeric-value additional-options) 16)
               (or tme-db-reformat-previous-options
                   tme-db-reformat-previous-labels)))
         (tme-db-entry-format
          (if additional-options
              (if use-previous-options
                  tme-db-reformat-previous-options
                (setq
                 tme-db-reformat-previous-options
                 (delq
                  nil
                  (list
                   (if (or
                        called-by-convert-alien
                        (y-or-n-p
                         "Realign entries (recommended for files not created by TME-db mode)? "))
                       'realign)
                   (if (y-or-n-p
                        "Remove empty optional and alternative fields? ")
                       'opts-or-alts)
                   (if (y-or-n-p
                        "Remove delimiters around pure numerical fields? ")
                       'numerical-fields)
                   (if (y-or-n-p (concat
                                  (if tme-db-comma-after-last-field
                                      "Insert"
                                    "Remove")
                                  " comma at end of entry? "))
                       'last-comma)
                   (if (y-or-n-p
                        "Replace double page dashes by single ones? ")
                       'page-dashes)
                   (if (y-or-n-p
                        "Force delimiters? ")
                       'delimiters)
                   (if (y-or-n-p
                        "Unify case of entry types and field names? ")
                       'unify-case)))))
            '(realign)))
         (labels
             (if additional-options
                 (if use-previous-options
                     tme-db-reformat-previous-labels
                   (setq
                    tme-db-reformat-previous-labels
                    (y-or-n-p "Generate automatically new reference labels? ")))))
         tme-db-autokey-edit-before-use
         (tme-db-sort-ignore-string-entries t)
         (start-point
          (if (tme-db-mark-active)
              (region-beginning)
            (progn
              (tme-db-beginning-of-first-entry)
              (tme-db-skip-to-valid-entry)
              (point))))
         (end-point
          (if (tme-db-mark-active)
              (region-end)
            (point-max)))
         (valid-tme-db-entry
          (concat
           "[ \t\n]+\\(@[ \t]*\\("
           (mapconcat
            (lambda (type)
              (concat "\\(" (car type) "\\)"))
            tme-db-entry-field-alist
            "\\|")
           "\\)\\)")))
    (save-restriction
      (narrow-to-region start-point end-point)
      (if (memq 'realign tme-db-entry-format)
          (progn
            (goto-char (point-min))
            (while (re-search-forward valid-tme-db-entry nil t)
              (replace-match "\n\\1"))))
      (goto-char start-point)
      (tme-db-progress-message "Formatting" 1)
      (tme-db-map-entries
       (lambda (current)
         (tme-db-progress-message)
         (tme-db-clean-entry labels labels)
         (if (memq 'realign tme-db-entry-format)
             (progn
               (tme-db-end-of-entry)
               (tme-db-delete-whitespace)
               (open-line 2)))))
      (tme-db-progress-message 'done))
    (if (and
         labels
         tme-db-maintain-sorted-entries
         (not called-by-convert-alien))
        (progn
          (tme-db-sort-buffer)
          (setq tme-db-keys nil)
          (tme-db-parse-keys nil t t)))
    (goto-char pnt)))

(defun tme-db-convert-alien (&optional do-additional-reformatting)
  "Converts an alien TME-db buffer to be fully usable by TME-db mode.
If a file doesn't confirm with some standards used by TME-db mode,
some of the high-level features of TME-db mode won't be available.
This function tries to convert current buffer to confirm with these standards.
With prefix argument DO-ADDITIONAL-REFORMATTING
non-nil, read options for reformatting entries from minibuffer."
  (interactive "*P")
  (message "Starting to validate buffer...")
  (sit-for 1 nil t)
  (goto-char (point-min))
  (while (re-search-forward "[ \t\n]+@" nil t)
    (replace-match "\n@"))
  (message
   "If errors occur, correct them and call `tme-db-convert-alien' again")
  (sit-for 5 nil t)
  (if (let ((tme-db-mark-active)
            tme-db-maintain-sorted-entries)
        (tme-db-validate))
      (progn
        (message "Starting to reformat entries...")
        (sit-for 2 nil t)
        (tme-db-reformat do-additional-reformatting t)
        (if tme-db-maintain-sorted-entries
            (progn
              (message "Starting to sort buffer...")
              (tme-db-sort-buffer)))
        (goto-char (point-max))
        (message "Buffer is now parsable. Please save it."))))

(defun tme-db-complete-string ()
  "Complete word fragment before point to longest prefix of a defined string.
If point is not after the part of a word, all strings are listed.
Remove surrounding delimiters if complete string could be expanded."
  (interactive "*")
  (tme-db-complete tme-db-strings t))

(defun tme-db-complete-key ()
  "Complete word fragment before point to longest prefix of a defined key.
If point is not after the part of a word, all keys are listed. This
function is most useful in completing crossref entries."
  (interactive "*")
  (if (not tme-db-keys)
      (tme-db-parse-keys nil t))
  (tme-db-complete tme-db-keys))

(defun tme-db-Article ()
  (interactive)
  (tme-db-entry "Article"))

(defun tme-db-Project ()
  (interactive)
  (tme-db-entry "Project"))

(defun tme-db-String ()
  (interactive)
  (if (not tme-db-keys)
      (tme-db-parse-keys nil t))
  (let ((key
         (if (and
              tme-db-maintain-sorted-entries
              (not tme-db-sort-ignore-string-entries))
             (completing-read
              "String key: " tme-db-keys nil nil nil 'tme-db-key-history))))
    (if (and
         tme-db-maintain-sorted-entries
         (not tme-db-sort-ignore-string-entries))
	(tme-db-find-entry-location key)
      (tme-db-move-outside-of-entry))
    (indent-to-column tme-db-entry-offset)
    (insert
     (concat
      "@String"
      (tme-db-entry-left-delimiter)
      (if (and
           tme-db-maintain-sorted-entries
           (not tme-db-sort-ignore-string-entries))
          key)
      " = "
      (tme-db-field-left-delimiter)
      (tme-db-field-right-delimiter)
      (tme-db-entry-right-delimiter)
      "\n"))
  (forward-line -1)
  (forward-char
   (if (and
        tme-db-maintain-sorted-entries
        (not tme-db-sort-ignore-string-entries))
       (+ (length "@String{") (length key) (length " = {"))
     (length "@String{")))))

(defun tme-db-Preamble ()
  (interactive)
  (tme-db-move-outside-of-entry)
  (indent-to-column tme-db-entry-offset)
  (insert
   "@Preamble"
   (tme-db-entry-left-delimiter)
   (tme-db-entry-right-delimiter)
   "\n")
  (forward-line -1)
  (forward-char 10))


;; Make TME-db a Feature

(provide 'tme-db)

;;; tme-db.el ends here
