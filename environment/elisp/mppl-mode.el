;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; MPPL mode for GNU Emacs, version 18.57.2
;;; written by Maurice J. LeBrun (mjl@dino.ph.utexas.edu)
;;; Version:         $Revision$
;;; Last Modified:   $Date$
;;;
;;; Based on Fortran mode for GNU Emacs by Michael D. Prange.
;;; Also code adapted from c-mode.el, c++-mode.el.
;;;
;;; MPPL is a Fortran preprocessor with the power of Ratfor and m4 rolled
;;; together into one consistent language.  The loop constructs are keyword
;;; delimited, similar in spirit to Fortran 90 (the syntax, however, is
;;; different).  
;;;
;;; MPPL is capable of processing basically vanilla Fortran code, but the
;;; mode defined here assumes that the code is in free-format, with no
;;; reliance on Fortran fixed format except that line numbers may start in
;;; the first column without affecting indentation.  If you want fixed
;;; format, use Fortran mode.  Similarly, unlike Fortran mode, "continue"
;;; is not understood as a do-loop delimiter; use "enddo" instead.
;;;
;;; To get started place the following command in your .emacs:
;;;
;;; (autoload 'mppl-mode   "mppl-mode" nil t)
;;;
;;; If you use "#" for mppl comments rather than "!", also put the following:
;;;
;;; (setq mppl-comment-char "#")
;;;
;;; Files with extension ".p", ".m", and ".i" are automatically recognized as
;;; MPPL files.
;;;
;;; Note this is just an initial effort, and there are a lot of improvements
;;; that could be made.  Feel free to send suggestions to me at
;;; mjl@dino.ph.utexas.edu.
;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;; Base settings.

(defvar mppl-comment-char "!"
  "*Character which is used for all mppl comments.")

(defvar mppl-do-indent 4
  "*Extra indentation applied to `do' blocks.")

(defvar mppl-if-indent 4
  "*Extra indentation applied to `if' blocks.")

(defvar mppl-bracket-indent 4
  "*Extra indentation applied to macro definition blocks.")

(defvar mppl-continuation-indent 4
  "*Extra indentation applied to `continuation' lines.")

(defvar mppl-minimum-statement-indent 0
  "*Minimum indentation for mppl statements.")

(defvar mppl-comment-indent-style 'fixed
  "*nil forces comment lines not to be touched,
'fixed produces fixed comment indentation to comment-column,
and 'relative indents to current mppl indentation plus comment-column.")

(defvar mppl-comment-line-column 2
  "*Indentation for text in comment lines.")

(defvar mppl-comment-indent-char " "
  "*Character to be inserted for MPPL comment indentation.
Normally a space.")

;; Utility functions for finding comments and continuations.

(setq comment-introducer-regexp "[!#]")

(setq comment-line-regexp
      (concat "[ \t]*" comment-introducer-regexp))

(setq trailing-comment-regexp
      (concat "[ \t]*\\(" comment-introducer-regexp ".*\\)?$"))

(setq continued-line-regexp
      (concat ".*[-+*(,&|~><=\]" trailing-comment-regexp))

(setq indented-line-regexp
      (concat ".*[-+*(,&|~><=\[]" trailing-comment-regexp))

;; Rest of mode setup info.

(defvar mppl-startup-message t
  "*Non-nil displays a startup message when mppl-mode is first called.")

(defconst mppl-mode-version "$Revision$")

(defvar mppl-mode-syntax-table nil
  "Syntax table in use in mppl-mode buffers.")

;; Set up syntax entries.
;; Probably this still isn't quite right.

(if mppl-mode-syntax-table
    ()
  (setq mppl-mode-syntax-table (make-syntax-table))
  (modify-syntax-entry ?\; "w"	mppl-mode-syntax-table)
  (modify-syntax-entry ?+  "."	mppl-mode-syntax-table)
  (modify-syntax-entry ?-  "."	mppl-mode-syntax-table)
  (modify-syntax-entry ?*  "."	mppl-mode-syntax-table)
  (modify-syntax-entry ?/  "."	mppl-mode-syntax-table)
  (modify-syntax-entry ?\' "\"" mppl-mode-syntax-table)
  (modify-syntax-entry ?\" "\"" mppl-mode-syntax-table)
  (modify-syntax-entry ?\\ "/"	mppl-mode-syntax-table)
  (modify-syntax-entry ?.  "w"	mppl-mode-syntax-table)
  (modify-syntax-entry ?!  "<"	mppl-mode-syntax-table)
  (modify-syntax-entry ?#  "<"	mppl-mode-syntax-table)
  (modify-syntax-entry ?\n ">"	mppl-mode-syntax-table)
  )

(defvar mppl-mode-map () 
  "Keymap used in mppl mode.")

(if mppl-mode-map
    ()
  (setq mppl-mode-map (make-sparse-keymap))
  (define-key mppl-mode-map "[" 	'electric-mppl-bracket)
  (define-key mppl-mode-map "]" 	'electric-mppl-bracket)
  (define-key mppl-mode-map "\C-c\C-c"  'mppl-comment-region)
  (define-key mppl-mode-map "\C-c\C-u"  'mppl-uncomment-region)
  (define-key mppl-mode-map "\e\C-a"	'beginning-of-mppl-subprogram)
  (define-key mppl-mode-map "\e\C-e"	'end-of-mppl-subprogram)
  (define-key mppl-mode-map "\e;"	'mppl-indent-comment)
  (define-key mppl-mode-map "\e\C-q"	'mppl-indent-subprogram)
  (define-key mppl-mode-map "\t"	'mppl-indent-line))

  (define-key mppl-mode-map "\e\C-h"	'mark-mppl-subprogram)
  (define-key mppl-mode-map "\C-c\C-p"	'mppl-previous-statement)
  (define-key mppl-mode-map "\C-c\C-n"	'mppl-next-statement)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; mppl-mode 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun mppl-mode ()
  "Major mode for editing mppl code.
Tab indents the current mppl line correctly. 

Variables controlling indentation style and extra features:

 mppl-do-indent
    Extra indentation within do blocks.  (default 4)
 mppl-if-indent
    Extra indentation within if blocks.  (default 4)
 mppl-bracket-indent
    Extra indentation within macro definition blocks.  (default 4)
 mppl-continuation-indent
    Extra indentation appled to continuation statements.  (default 4)
 mppl-comment-line-column
    Amount of indentation for text within full-line comments. (default 6)
 mppl-comment-indent-style
    nil    means don't change indentation of text in full-line comments,
    fixed  means indent that text at column mppl-comment-line-column
    relative  means indent at mppl-comment-line-column beyond the
 	      indentation for a line of code.
    Default value is fixed.
 mppl-comment-char
    Character (\"!\" or \"#\") to be inserted for a full-line or in-line
    comment (default \"!\").  Both characters are always checked for when
    determining if an existing line is a comment.
 mppl-comment-indent-char
    Character to be inserted instead of space for full-line comment
    indentation.  (default is a space)
 mppl-minimum-statement-indent
    Minimum indentation for mppl statements. (default 0)
 mppl-startup-message
    Set to nil to inhibit message first time mppl-mode is used.

Turning on MPPL mode calls the value of the variable mppl-mode-hook 
with no args, if that value is non-nil.
\\{mppl-mode-map}"

;; MPPL mode startup code

  (interactive)
  (kill-all-local-variables)
  (if mppl-startup-message
      (message "Emacs MPPL mode %s." mppl-mode-version))
  (setq mppl-startup-message nil)
  (set-syntax-table mppl-mode-syntax-table)
  (set (make-local-variable 'indent-line-function) 'mppl-indent-line)
  (set (make-local-variable 'comment-indent-hook) 'mppl-comment-hook)
  (require 'filladapt)

;; Some of this needs to be changed

  (set (make-local-variable 'comment-line-start)
       (concat mppl-comment-char " "))

  (set (make-local-variable 'comment-start)
       (concat mppl-comment-char " "))

  (set (make-local-variable 'mppl-comment-region)
       (concat mppl-comment-char " "))

  (set (make-local-variable 'require-final-newline) t)
  (set (make-local-variable 'indent-tabs-mode) nil)

  (use-local-map mppl-mode-map)
  (setq mode-name "MPPL")
  (setq major-mode 'mppl-mode)
  (run-hooks 'mppl-mode-hook))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; calculate-mppl-indent
;;
;; This function is the real meat of the indentation code, and is necessarily
;; fairly complex.  It calculates the indentation column for the current
;; line, taking into account indentation for if/endif, do/enddo, [], and line
;; continuation.
;;
;; Line numbers MUST start at column 1 if you want them to not affect the
;; indentation level.
;;
;; Continuation markers at column 5 are not recognized, period!  Use MPPL
;; end-of-line continuation markers instead.  Except for "/" (which
;; unfortunately can be the tail end of a DATA statement), trailing operators
;; signify an implicit continuation (a trailing "\" is explicit
;; continuation).
;;
;; Continuation is a bit tricky, because if the current line is a
;; continuation there may be multiple continuations, and we don't want
;; multiple indents.  Also, we need to "undo" the continued line indent after
;; the last continued line.
;;
;; Also need to allow for trailing comments.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun calculate-mppl-indent ()
  "Calculates the mppl indent column based on previous lines."
  (let (icol first-statement (case-fold-search t))

;; Excursion: check previous line for: if-then, do, "[", and continuation.

    (save-excursion
      (setq first-statement (mppl-previous-statement))
      (if first-statement
	  (setq icol mppl-minimum-statement-indent)
	(progn
	  (if (= (point) (point-min))
	      (setq icol mppl-minimum-statement-indent)
	    (setq icol (mppl-current-line-indentation)))
	  (skip-chars-forward "0-9")
	  (skip-chars-forward " \t")
	  (cond

;; if

	   ((looking-at "if[ \t]*(")
	    (if (or
		 (looking-at ".*)[ \t]*then\\b[ \t]*[^ \t(=a-z0-9]")

;; Multi-line if-then -- search forward for "then".

		 (let (then-test)
		   (while
		       (and
			(= (forward-line 1) 0)
			(looking-at "     [^ 0]")
			(not
			 (setq then-test
			       (looking-at
				".*then\\b[ \t]*[^ \t(=a-z0-9]")
			       )
			 )))
		   then-test))
		(setq icol (+ icol mppl-if-indent))))

;; Open bracket, not followed by a close bracket.

	   ((looking-at ".*\\[[^]]*$")
	    (setq icol (+ icol mppl-bracket-indent)))

;; else or else if

	   ((looking-at "\\(else\\|elseif\\)\\b")
	    (setq icol (+ icol mppl-if-indent)))

;; do

	   ((looking-at "do\\b")
	    (setq icol (+ icol mppl-do-indent)))

;; Continuation marker found at end of line.  Check if the previous line is a
;; continuation, and if NOT only then indent.

	   ((looking-at continued-line-regexp)
	    (progn
	      (setq first-statement (mppl-previous-statement))
	      (if first-statement
		  (setq icol (+ icol mppl-continuation-indent))
		(progn
		  (if (= (point) (point-min))
		      (setq icol (+ icol mppl-continuation-indent))
		    (skip-chars-forward "0-9")
		    (skip-chars-forward " \t")
		    (cond
		     ((not (looking-at continued-line-regexp))
		      (setq icol (+ icol mppl-continuation-indent)))))))
	      ))
	   ))))

;; Excursion: check current line for: enddo, endif, "[", "]".

    (save-excursion
      (beginning-of-line)
      (cond
       ((looking-at "[ \t]*$"))
       ((looking-at comment-line-regexp)
	(cond
	 ((eq mppl-comment-indent-style 'relative)
	  (setq icol (+ icol mppl-comment-line-column)))
	 ((eq mppl-comment-indent-style 'fixed)
	  (setq icol mppl-comment-line-column))))
       (first-statement)
       (t
	(skip-chars-forward "0-9")
	(skip-chars-forward " \t")
	(cond

;; end if

	 ((looking-at "end[ \t]*if\\b")
	  (setq icol (- icol mppl-if-indent)))

;; Open brackets by themselves on a line are probably part of a continuation.
;; It looks better in this case to defer indentation until the following
;; line.

	 ((looking-at
	   (concat "[ \t]*\\[" trailing-comment-regexp))
	  (setq icol (- icol mppl-bracket-indent)))

;; Close bracket.

	 ((looking-at "].*")
	  (setq icol (- icol mppl-bracket-indent)))

;; else if

	 ((looking-at "\\(else\\|elseif\\)\\b")
	  (setq icol (- icol mppl-if-indent)))

;; end do

	 ((looking-at "end[ \t]*do\\b")
	  (setq icol (- icol mppl-do-indent)))
	 ))))

;; Excursion: check current line for non-continuation.

    (save-excursion
      (beginning-of-line)
      (cond
       ((looking-at "[ \t]*$"))
       ((looking-at comment-line-regexp)
	(cond
	 ((eq mppl-comment-indent-style 'relative)
	  (setq icol (+ icol mppl-comment-line-column)))
	 ((eq mppl-comment-indent-style 'fixed)
	  (setq icol mppl-comment-line-column))))
       (first-statement)
       (t
	(skip-chars-forward "0-9")
	(skip-chars-forward " \t")
	(cond

;; If we are on a new statement after a continued line, need to back up one.
;; I compare to indented-line-regexp instead of continued-line-regexp.  The
;; former adds an open bracket at the end of the line.  This way, when the
;; open bracket appears by itself on a line, indentation gets deferred until
;; after the open bracket, resulting in macro enclosures of the form:
;;
;;	ifelse(A, B,
;;	[
;;	    george
;;	])
;;
;; If the open bracket is instead put at the end of the line containing the
;; ifelse, the indentation will be unchanged.

	 (t
	  (progn
	    (setq first-statement (mppl-previous-statement))
	    (if (not first-statement)
		(progn
		  (skip-chars-forward "0-9")
		  (skip-chars-forward " \t")
		  (if (not (looking-at indented-line-regexp))

;; Current line is not a continuation.  Back up to the last non-continued
;; line, and set the indentation accordingly.
;;
;; Right now I take the easy way out and back up one indentation level
;; instead, if the previous line was a continued statement.  This will be
;; wrong if for some reason the preceding line was indented by other than the
;; usual amount.  

		      (progn
			(setq first-statement (mppl-previous-statement))
			(if (not first-statement)
			    (progn
			      (skip-chars-forward "0-9")
			      (skip-chars-forward " \t")
			      (cond
			       ((looking-at continued-line-regexp)
				(setq icol
				      (- icol mppl-continuation-indent))))
			      ))
			))
		  ))
	    ))))))

    (max mppl-minimum-statement-indent icol)))


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; electric-mppl-bracket
;;
;; automatically goes back one bracket indentation level when a close
;; bracket "]" is entered, and momentarily jumps back to the matching
;; open bracket.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun electric-mppl-bracket (arg)
  "Insert character and correct line's indentation."
  (interactive "P")
  (let (insertpos)
    (if (and (not arg)
	     (eolp)
	     (or (save-excursion
		   (skip-chars-backward " \t")
		   (bolp)) ))
	(progn
	  (insert last-command-char)
	  (mppl-indent-line)
	  (save-excursion
	    (if insertpos (goto-char (1+ insertpos)))
	    (delete-char -1))))
    (if insertpos
	(save-excursion
	  (goto-char insertpos)
	  (self-insert-command (prefix-numeric-value arg)))
      (self-insert-command (prefix-numeric-value arg)))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; current-line-indentation
;;
;; Calculates the indentation of the current line.  Usually used on
;; excursions to past lines in determining the indentation to use next.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun mppl-current-line-indentation () 
  "Indentation of current line.
This is the column position of the first non-whitespace character aside from
the line number.  For comment lines, returns indentation of the first
non-indentation text within the comment."

  (save-excursion
    (beginning-of-line)
    (cond ((looking-at comment-line-regexp)
	   (goto-char (match-end 0))
	   (skip-chars-forward
	     (if (stringp mppl-comment-indent-char)
		 mppl-comment-indent-char
	         (char-to-string mppl-comment-indent-char)))))
    ;; Move past line number or whitespace.
    (skip-chars-forward "0-9")
    (skip-chars-forward " \t")
    (current-column)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; mppl-comment-region, mppl-uncomment-region
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun mppl-comment-region (beg end)
  "Comment out all lines in a region between mark and current point by
inserting comment-start in front of each line."
  (interactive "*r")
  (save-excursion
    (save-restriction
      (narrow-to-region
       (progn (goto-char beg) (beginning-of-line) (point))
       (progn (goto-char end) (or (bolp) (forward-line 1)) (point)))
      (goto-char (point-min))
      (while (not (eobp))
	(insert comment-start)
	(forward-line 1)))))

(defun mppl-uncomment-region (beg end)
  "Uncomment all lines in region between mark and current point by deleting
the leading comment from each line, if any."
  (interactive "*r")
  (save-excursion
    (save-restriction
      (narrow-to-region
       (progn (goto-char beg) (beginning-of-line) (point))
       (progn (goto-char end) (forward-line 1) (point)))
      (goto-char (point-min))
      (while (not (eobp))
	(if (looking-at comment-start)
	    (delete-region (match-beginning 0) (match-end 0)))
	(forward-line 1)))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; mppl-indent-to-column
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun mppl-indent-to-column (col)
  "Indents current line with spaces to column COL.
Note that a line which has a number as the first non-whitespace character is
a numbered line."
  (save-excursion
    (beginning-of-line)
    (if (looking-at comment-line-regexp)
	(if mppl-comment-indent-style
	    (let ((char (if (stringp mppl-comment-indent-char)
			    (aref mppl-comment-indent-char 0)
			    mppl-comment-indent-char)))
	      (goto-char (match-end 0))
	      (delete-horizontal-regexp (concat " \t" (char-to-string char)))
	      (insert-char char (- col (current-column)))))
      (cond ((eobp))
	    ((looking-at "[0-9]+")
	       (let ((extra-space (- 5 (- (match-end 0) (point)))))
		 (if (< extra-space 0)
		     (message "Warning: line number exceeds 5-digit limit."))
		 (skip-chars-forward "0-9"))))

;; Point is now after any line number.
;; Put body of statement where specified.

      (delete-horizontal-space)
      (indent-to col))))

(defun delete-horizontal-regexp (chars)
  "Delete all characters in CHARS around point.
CHARS is like the inside of a [...] in a regular expression
except that ] is never special and \ quotes ^, - or \."
  (interactive "*s")
  (skip-chars-backward chars)
  (delete-region (point) (progn (skip-chars-forward chars) (point))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; mppl-indent-line
;;
;; Note: the fortran.el code has a recursive path between fortran-indent-line
;; and fortran-indent-comment that I couldn't figure out.  It seems though
;; that this routine didn't need some of the comment handling code so I left
;; it out and it seems to work fine now.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun mppl-indent-line ()
  "Indents current mppl line based on its contents and on previous lines."
  (interactive)
  (let ((cfi (calculate-mppl-indent)))
    (save-excursion
      (beginning-of-line)
      (if (not (= cfi (mppl-current-line-indentation)))
	  (mppl-indent-to-column cfi)))

;; Never leave point in left margin.

    (if (< (current-column) cfi)
	(move-to-column cfi))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; mppl-indent-subprogram
;;
;; This really does work!  Leaves trailing blanks on empty lines, though.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun mppl-indent-subprogram ()
  "Properly indents the MPPL subprogram which contains point."
  (interactive)
  (save-excursion
    (mark-mppl-subprogram)
    (message "Indenting subprogram...")
    (indent-region (point) (mark) nil))
  (message "Indenting subprogram...done."))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Comment-line indentation handling code.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun mppl-comment-hook ()
  (save-excursion
    (skip-chars-backward " \t")
    (max (+ 1 (current-column))
	 comment-column)))

(defun mppl-indent-comment ()
  "Create comment on current line.
Existing comments are not affected.  If the line has no comment, an inline
comment is inserted at the column given by \\[comment-start].  Otherwise, a
separate-line comment is inserted, on this line or on a new line inserted
before this line if this line is not blank."
  (interactive)
  (beginning-of-line)

  (cond ((looking-at comment-line-regexp)
;	 (mppl-indent-line)
	 )

;; No existing comment so insert inline comment, unless line is now blank.

	((and comment-start (not (looking-at "^[ \t]*$")))
	 (end-of-line)
	 (delete-horizontal-space)
	 (indent-to (mppl-comment-hook))
	 (insert comment-start))

;; Else insert separate-line comment, making a new line if nec.

	(t
	 (if (looking-at "^[ \t]*$")
	     (delete-horizontal-space)
	   (beginning-of-line)
	   (insert "\n")
	   (forward-char -1))
	 (insert comment-line-start)
	 (insert-char (if (stringp mppl-comment-indent-char)
			  (aref mppl-comment-indent-char 0)
			  mppl-comment-indent-char)
		      (- (calculate-mppl-indent) (current-column))))))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Routines to deal with the beginning/end of a subroutine, as marked by the
;; END statement.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun beginning-of-mppl-subprogram ()
  "Moves point to the beginning of the current mppl subprogram."
  (interactive)
  (let ((case-fold-search t))
    (beginning-of-line -1)
    (re-search-backward "^[ \t0-9]*end\\b[ \t]*[^ \t=(a-z]" nil 'move)
    (if (looking-at "^[ \t0-9]*end\\b[ \t]*[^ \t=(a-z]")
	(forward-line 1))))

(defun end-of-mppl-subprogram ()
  "Moves point to the end of the current mppl subprogram."
  (interactive)
  (let ((case-fold-search t))
    (beginning-of-line 2)
    (re-search-forward "^[ \t0-9]*end\\b[ \t]*[^ \t=(a-z]" nil 'move)
    (goto-char (match-beginning 0))
    (forward-line 1)))

(defun mark-mppl-subprogram ()
  "Put mark at end of mppl subprogram, point at beginning. 
The marks are pushed."
  (interactive)
  (end-of-mppl-subprogram)
  (push-mark (point))
  (beginning-of-mppl-subprogram))
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Go to previous/next mppl statement.  These don't handle continued lines at
;; present -- I just hacked out the old fortran continuation line handling.
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(defun mppl-previous-statement ()
  "Moves point to beginning of the previous mppl statement.
Returns 'first-statement if that statement is the first
non-comment MPPL statement in the file, and nil otherwise.
At present does not handle continued lines."
  (interactive)
  (beginning-of-line)
  (while (and (setq not-first-statement (= (forward-line -1) 0))
	      (or (looking-at comment-line-regexp)
		  (looking-at "[ \t]*$"))))
  (cond 
   ((not not-first-statement)
    'first-statement)))

(defun mppl-next-statement ()
  "Moves point to beginning of the next mppl statement.
Returns 'last-statement if that statement is the last
non-comment MPPL statement in the file, and nil otherwise.
At present does not handle continued lines."
  (interactive)
  (let (not-last-statement)
    (beginning-of-line)
    (while (and (setq not-last-statement (= (forward-line 1) 0))
 		(or (looking-at comment-line-regexp)
 		    (looking-at "[ \t]*$"))))
    (if (not not-last-statement)
 	'last-statement)))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(autoload 'mppl-mode "mppl" "Edit file FILENAME in mppl-mode." t)
(setq auto-mode-alist (cons '("\\.i$" . mppl-mode) auto-mode-alist))
(setq auto-mode-alist (cons '("\\.p$" . mppl-mode) auto-mode-alist))
(setq auto-mode-alist (cons '("\\.m$" . mppl-mode) auto-mode-alist))
