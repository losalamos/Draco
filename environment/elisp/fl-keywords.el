;; $Id$
;; 
;; Add additional keywords for font-lock modes.
;;
;;======================================================================

(require 'font-lock)

(make-face 'cvs-unknown-face)
(make-face 'cvs-modified-face)

(defconst cvs-font-lock-keywords
  '(("^In directory \\(.+\\)$" 1 cvs-header-face)
    ("^[* ] \\w+ +\\(ci\\)" 1 cvs-status-face)
;    ("^[* ] \\(Conflict\\|Merged\\)" 1 cvs-status-face)
    ("^[* ] \\w+ +\\(ci +\\)?\\(.+\\)$" 2 cvs-filename-face)
    ("^[* ] \\(Unknown\\)" 1 cvs-unknown-face)
    ("^[* ] \\(Modified\\)" 1 cvs-modified-face)
    ))

; Should probably make some faces for working with MPPL.  Maybe
; different faces for MPPL vs. native Fortran features, etc.

(make-face 'tpc-stuff-face)
(make-face 'mppl-stuff-face)

(defconst mppl-font-lock-keywords
  '(
; Fortran keywords
    ("[ 	]+\\(if\\|then\\|else\\|endif\\|goto\\|do\\|enddo\\|continue\\|subroutine\\|entry\\|return\\|end\\|program\\|call\\|external\\|read\\|write\\)" . font-lock-keyword-face)
; Fortran types
    ("[ 	]+\\(integer\\|real\\|logical\\)\\(*[0-9]+\\)?" . font-lock-type-face)
; MPPL keywords
    ("[ 	]+\\(define\\|include\\|eval\\)[ 	]+" . font-lock-keyword-face)
; Function names
    ("[ 	]+external[ 	]+\\(\\([A-Za-z_0-9]+\\)\\(,[ 	]*\\)?\\)+" 1 font-lock-function-name-face)
    ("[ 	]+subroutine[ 	]+\\([A-Za-z_0-9]+\\)" 1 font-lock-function-name-face)
    ("[ 	]+entry[ 	]+\\([A-Za-z_0-9]+\\)" 1 font-lock-function-name-face)
    ("[ 	]+call[ 	]+\\([A-Za-z_0-9]+\\)" 1 font-lock-function-name-face)
; TPC things
    ("[ 	]+\\(declare_[A-Za-z_]+\\)" 1 tpc-stuff-face)
    ("[ 	]+\\(include_[A-Za-z_]+\\)" 1 tpc-stuff-face)
    ("[ 	]+\\(shortloop\\|quit\\|transmit\\|REAL\\|dbug_enter\\|dbug_return\\)" . tpc-stuff-face)
    ))

(defconst  tex-font-lock-keywords-2 (list 
  '("\\(\\\\[a-zA-Z@]+\\)" 1 font-lock-keyword-face t) 
  '("\\\\\\(begin\\|end\\){\\([^}]+\\)}" 2 font-lock-function-name-face t)
;  '("\\(\\\\\\w+\\){\\(\\w+\\)}" 2 font-lock-function-name-face t) 
;  '("\\(\\\\\\w+\\){\\(\\w+\\)}{\\(\\w+\\)}" 3 font-lock-function-name-face t)
;  '("\\(\\\\\\w+\\){\\(\\w+\\)}{\\(\\w+\\)}{\\(\\w+\\)}" 4 
;    font-lock-function-name-face t) 
  '("\\([^\\\\]\\|^\\)\\$+\\([^$]*\\)\\$+" 2 font-lock-doc-string-face t) 
  '("{\\\\\\(em\\|it\\|sl\\|sf\\|tt\\)\\([^}]+\\)}" 2 italic t)
;  '("{\\\\it\\([^}]+\\)}" 1 italic t)
  '("{\\\\bf\\([^}]+\\)}" 1 bold t) 
  '("^[ 	]*\\\\def[\\\\@]\\(\\w+\\)\\W" 1 font-lock-function-name-face
t) 
  ;; make argument of \bm (my boldmath macro) bold
  '("\\\\bms*\\({\\([^}]*\\)}\\| \\(.\\)\\|\\(\\\\[a-zA-Z]+\\)\\)" 1 
	bold t) 
  '("\\\\section{\\([^}]+\\)}" 1 font-lock-type-face t)
  '("\\\\subsection{\\([^}]+\\)}" 1 bold-italic t)
  '("\\\\ref{\\([^}]+\\)}" 1 font-lock-string-face t)
  '("\\\\cite{\\([^}]+\\)}" 1 font-lock-string-face t)
  '("\\\\label{\\([^}]+\\)}" 1 font-lock-string-face t)
  '("\\\\bibitem{\\([^}]+\\)}" 1 font-lock-string-face t)
  )
  "More fancy version of tex-font-lock-keywords for highlighting in LaTeX mode"
)

(make-face 'gts-macros-face)
(make-face 'font-lock-kull-macros-face)

(let ((preprocessor '(
	 ("^#[ 	]*\\(undef\\|define\\|ifdef\\|ifndef\\)\\s +\\(\\(\\sw\\|\\s_\\)+\\)" 2 
	  font-lock-function-name-face t)
     ("^\\(#\\s *\\w+\\)\\>" 1 font-lock-doc-string-face t)
	 ("^#[ 	]*\\(undef\\|define\\)\\s +\\(\\(\\sw\\|\\s_\\)+\\)" 2 
	  font-lock-function-name-face t)	 
	 ("^#[ 	]*include[ 	]+\\(<[^>\"
]+>\\)" 1 font-lock-string-face t))))

  (defconst c-font-lock-keywords-3
   (append
    preprocessor
    '(
      ;; key words
      ("\\(^\\|[^_]\\)\\<\\(typedef\\|struct\\|union\\|enum\\|\\)\\>\\([^_]\\|$\\)" 2 bold-italic)
	  ("\\(^\\|[^_]\\)\\<\\(register\\|unsigned\\|extern\\|static\\|auto\\)\\>\\([^_]\\|$\\)" 2 bold)
      ("\\(^\\|[^_]\\)\\<\\(return\\|goto\\|if\\|else if\\|else\\|case\\|default\\|switch\\|break\\|continue\\|while\\|do\\|for\\)\\>\\([^_]\\|$\\)" 2 font-lock-keyword-face)
      ;; function decls 
	  ("^\\(\\(\\sw\\|\\s_\\|[:~*]\\)+[ 	]+\\)?\\(\\(\\sw\\|\\s_\\|[:~*]\\)+[ 	]+\\)?\\(\\(\\sw\\|\\s_\\|[:~*]\\)+[ 	]+\\)?\\([*&]+[ 	]*\\)?\\(\\(\\sw\\|\\s_\\|[:~*]\\)+\\)[ 	]*(" 8 font-lock-function-name-face t)
      ("^\\(\\(\\w\\|[$_]\\)+[ \t]*::[ \t]*\\)?\\(\\(\\w\\|[$_]\\)+\\|operator.*\\)\\s *\\(\\(\\w\\|[$_]\\)+\\s *((\\|(\\)[^)]*)+" . font-lockfunction-name-face)
      ;; datatype -- black magic regular expression
	  ("\\(^\\|[({,]\\)\\s *\\(\\(typedef\\|const\\|register\\|volatile\\|unsigned\\|extern\\|static\\|auto\\)\\s +\\)*\\(\\(\\w\\|\\s_\\)+\\)\\>\\([ 	*]*&\\s *\\|\\s +\\**\\|\\**\\s +\\)\\<[a-zA-Z_]\\(\\w\\|\\s_\\)*\\>" 4 font-lock-type-face)
   ; strings should be done automatically but ocassionally 
   ; font-lock uses the doc-string-face rather than plain string-face
   ("\"\\(\\(\\\\\"\\|[^\"]\\)+\\)\"" 1 font-lock-string-face t)
;   ("'\\(\\?.\\)'" 1 font-lock-string-face t)
	 "Yet another form of highlighting for c-mode")
	))
 
  (defconst c++-font-lock-keywords-gmf
    (append
     preprocessor
     '(
;;; key words
       ("\\(^\\|[^_]\\)\\<\\(template\\|typedef\\|struct\\|union\\|class\\|enum\\|public\\|private\\|protected\\)\\>\\([^_]\\|$\\)" 2 font-lock-keyword-face)
       ("\\(^\\|[^_]\\)\\<\\(inline\\|virtual\\|friend\\|const\\|register\\|volatile\\|unsigned\\|extern\\|static\\|auto\\)\\>\\([^_]\\|$\\)" 2 font-lock-keyword-face)
       ("\\(^\\|[^_]\\)\\<\\(return\\|goto\\|if\\|else if\\|else\\|case\\|default\\|switch\\|break\\|continue\\|while\\|do\\|for\\|throw\\|try\\|catch\\|delete\\|new\\)\\>\\([^_]\\|$\\)" 2 font-lock-keyword-face)
       ("\\(^\\|[^_]\\)\\<\\(explicit\\|mutable\\|typename\\|using\\|namespace\\)\\>\\([^_]\\|$\\)" 2 font-lock-keyword-face)
;;; casts.
       ("\\(^\\|[^_]\\)\\<\\(static_cast\\|dynamic_cast\\|reinterpret_cast\\|const_cast\\)\\>\\([^_]\\|$\\)" 2 font-lock-keyword-face)
       ("\\(^\\|[^_]\\)\\<\\(static_cast\\|dynamic_cast\\|reinterpret_cast\\|const_cast\\)\\><\\([A-Za-z0-9:_* 	]+\\)>" 3 font-lock-type-face)

       ;; function decls 
;;;	  ("^\\(inline\\s +\\)?\\(\\(\\sw\\|\\s_\\|::\\|[~*]\\)+[ 	]+\\)?\\(\\(\\sw\\|\\s_\\|::\\|[~*]\\)+[ 	]+\\)?\\(\\(\\sw\\|\\s_\\|::\\|[~*]\\)+[ 	]+\\)?\\([*&]+[ 	]*\\)?\\(\\(\\sw\\|\\s_\\|::\\|[~*]\\)+\\|operator[><+-*/=~^!]+\\)[ 	]*(" 9 font-lock-function-name-face t)
;;;	  ("^\\(\\s *\\(virtual\\|friend\\)\\s +\\)?[ 	]*\\([A-Za-z0-9_<>, *&]+\\)\\**[ 	]+\\**\\([A-Za-z0-9_]+\\|operator[][<>\\+\\-\\*\\=/~!^]+\\)[ 	]*(" 3 font-lock-type-face)
;;;	  ("^\\(\\s *\\(virtual\\|friend\\)\\s +\\)?[ 	]*\\([A-Za-z0-9_<>, *&]+\\)\\**[ 	]+\\**\\([A-Za-z0-9_]+\\|operator[][<>\\+\\-\\*\\=/~!^]+\\)[ 	]*(" 4 font-lock-function-name-face)

       ;; These are supposed to do function prototypes.  The first is
       ;; supposed to do the return type, and the second the function
       ;; name.

       ;;           type qualifiers                                               return type                                 function name or operator name

      ("^\\(\\s *\\(virtual\\|inline\\|friend\\|const\\|volatile\\)\\s +\\)*[ 	]*\\([A-Za-z_][A-Za-z0-9_:]*\\(\\s *<\\s *[A-Za-z_0-9][A-Za-z0-9_<>]*\\(,\\s *[A-Za-z_][A-Za-z0-9_<>]*\\)*\\s *>\\)?\\)&?\\**[ 	]+\\**\\([A-Za-z_][A-Za-z0-9_]+\\|operator[^A-Za-z0-9_]+\\)[ 	]*(" 3 font-lock-type-face)
      ("^\\(\\s *\\(virtual\\|inline\\|friend\\|const\\|volatile\\)\\s +\\)*[ 	]*\\([A-Za-z_][A-Za-z0-9_:]*\\(\\s *<\\s *[A-Za-z_0-9][A-Za-z0-9_<>]*\\(,\\s *[A-Za-z_][A-Za-z0-9_<>]*\\)*\\s *>\\)?\\)&?\\**[ 	]+\\**\\([A-Za-z_][A-Za-z0-9_]+\\|operator[^A-Za-z0-9_]+\\)[ 	]*(" 6 font-lock-function-name-face)
      ;; return type was "[A-Za-z0-9_<>,*&]+"
      ;; change to be same as function argument types..
      ;; "[A-Za-z_][A-Za-z0-9_:]*\\(\\s *<\\s *[A-Za-z_][A-Za-z0-9_<>]*\\(,\\s *[A-Za-z_][A-Za-z0-9_<>]*\\)*\\s *>\\)?"

      ;; Do the "=0" on pure virtual function declarations.
      (")\\s *\\(\\s +const\\s +\\)?\\(=\\s *0\\)\\s *;" 2 font-lock-doc-string-face)

; This next one was provided, I'm generalizing to what is show in next line.
;     ("^class\\s +\\(\\(\\sw\\|\\s_\\)+\\)\\>" 1 font-lock-function-name-face)
      ("\\(class\\|new\\|public\\|protected\\|private\\)\\s +\\(\\(\\sw\\|\\s_\\)+\\(<[A-Za-z0-9_, *&]+>\\)?\\)" 2 font-lock-type-face)

;      ("^\\(\\(\\w\\|[$_]\\)+[ \t]*::[ \t]*\\)?\\(\\(\\w\\|[$_]\\)+\\|operator.*\\)\\s *\\(\\(\\w\\|[$_]\\)+\\s *((\\|(\\)[^)]*)+" . font-lock-function-name-face)
      ;; datatype -- black magic regular expression
; This next one is the original, the following are attempts to improve it.
;	  ("\\(^\\|[({,]\\)\\s *\\(\\(virtual\\|friend\\|typedef\\|inline\\|const\\|register\\|volatile\\|unsigned\\|extern\\|static\\|auto\\)\\s +\\)*\\(\\(\\w\\|\\s_\\)+\\)\\>\\([ 	*]*&\\s *\\|\\s +\\**\\|\\**\\s +\\)\\<[a-zA-Z_]\\(\\w\\|\\s_\\)*\\>" 4 font-lock-type-face)

; next two lines were enabled.
      ;; This next one is supposed to do function argument types and variable declarations.
      ;;   \n | ( | , | const                 typename         
      ("\\(^\\|[(,]\\|const\\s \\)\\s *\\([A-Za-z_][A-Za-z0-9_:]*\\(\\s *<\\s *\\([A-Za-z_0-9]+\\|[A-Za-z_][A-Za-z0-9_<>]*\\)\\(,\\s *[A-Za-z_][A-Za-z0-9_<>]*\\)*\\s *>\\)?\\)[ 	*&]+[A-Za-z_][A-Za-z0-9_]*\\s *\\([=),;[]\\|//\\|/\\*\\)" 2 font-lock-type-face)

      ;; This next one is supposed to do function return types.
;	  ("^\\s *\\(const\\s +\\)?\\([A-Za-z_][A-Za-z0-9_<>]+\\)\\s +[A-Za-z0-9_<>:]+\\s *(" 2 font-lock-type-face)

      ;; This next one is supposed to be a template specifier
      ;; catch-all. This would be useful for doing template base
      ;; class names in member initialization lists, for example.

;;;	  ("[A-Za-z0-9_]+<[A-Za-z0-9_ 	,*&]+>" . font-lock-type-face)
;;;	  ("[A-Za-z0-9_]+<[A-Za-z0-9_ 	,*&<>]+>" . font-lock-type-face)
      ("[A-Za-z_][A-Za-z0-9_]*<\\s *[A-Za-z_0-9][A-Za-z0-9_ 	,*<>]*>" . font-lock-type-face)
      ;; The above should be generalized to be like the other
      ;; template type recognition strings, as above in function
      ;; argument type and return type detection.

      ("\\([A-Za-z0-9_<>,]+\\)::\\(~?[A-Za-z0-9_]+\\)" 1 font-lock-type-face)
;	  ("\\([A-Za-z0-9_<>,]+\\)::\\(~?[A-Za-z0-9_]+\\|operator[][<>\\+\\-\\*\\=/~!^]+\\)[ 	]*(" 2 font-lock-function-name-face)
      ("\\([A-Za-z0-9_<>,]+\\)::\\(~?[A-Za-z0-9_]+\\|operator[^A-Za-z0-9_]+\\)[ 	]*(" 2 font-lock-function-name-face)


; this next one didn't work
;	  ("\\(^\\|[({,]\\)\\s *\\(\\(virtual\\|friend\\|typedef\\|inline\\|const\\|register\\|volatile\\|unsigned\\|extern\\|static\\|auto\\)\\s +\\)*\\([A-Z]+\\)\\>\\([ 	*]*&\\s *\\|\\s +\\**\\|\\**\\s +\\)\\<[a-zA-Z_]\\(\\w\\|\\s_\\)*\\>" 4 font-lock-type-face)
   ; strings should be done automatically but ocassionally 
   ; font-lock uses the doc-string-face rather than plain string-face
;   ("\"\\(\\(\\\\\"\\|[^\"]\\)+\\)\"" 1 font-lock-string-face t)

; This one is supposed to get return types of member functions
;     ("\\([A-Za-z0-9_]+\\(<[A-Za-z_, <>*&]+>\\)?\\)\\**&*\\s +\\**&*\\([A-Za-z0-9_]+\\)\\(<[A-Za-z_, <>*&]+>\\)?::\\([A-Za-z0-9_]+\\|operator[][<>\\+\\-\\*\\=/~!^]+\\)\\s *(" 1 font-lock-type-face)
      ("\\([A-Za-z0-9_]+\\(<[A-Za-z_, <>*&]+>\\)?\\)\\**&*\\s +\\**&*\\([A-Za-z0-9_]+\\)\\(<[A-Za-z_, <>*&]+>\\)?::\\([A-Za-z0-9_]+\\|operator[^\\sw]+\\)\\s *(" 1 font-lock-type-face)

      ;; GTS macros

      ("\\<\\(Assert\\|Catch\\|Throw\\|Tcl_Assert\\)\\>" 1 gts-macros-face)

      ;; Kull macros (contract programming model).

      ("\\<\\(Check\\|Require\\|Ensure\\|Remember\\|Assert\\|Insist\\)\\>" 1 font-lock-kull-macros-face)
      ))
    "Geoff's C++ font lock keywords--guaranteed to make you live longer!"
    )
  )


(defconst compilation-font-lock-keywords '(
   ("^[-_.\"A-Za-z0-9/+]+\\(:\\|, line \\)[0-9]+: \\([wW]arning:\\).*$" .
    font-lock-keyword-face)
   ("^[-_.\"A-Za-z0-9/+]+\\(: *\\|, line \\)[0-9]+:.*$" . font-lock-function-name-face)
   ("^[^:\n]+-[a-zA-Z][^:\n]+$" . font-lock-doc-string-face)
   ("\\(^[-_.\"A-Za-z0-9/+]+\\)\\(: *\\|, line \\)[0-9]+" 1 font-lock-string-face t)
   ("^[-_.\"A-Za-z0-9/+]+\\(: *[0-9]+\\|, line [0-9]+\\)" 1 bold t)
   )
  "Keywords for fontifying compilation buffers"
)

;(defconst makefile-font-lock-keywords '(
;   ("^#.*$" . font-lock-comment-face)
;   ("[^$]#.*$" . font-lock-comment-face)
;   ;; rules
;   ("^\\([^ \t\n]*%?[^ \t\n]*[ \t]*::?\\)[ \t]" 1 font-lock-type-face t)
;   ("^\\(\\.[A-Za-z][A-Za-z]?\\..[ \t]*::?\\)" 1 font-lock-type-face t)
;   ("^[^ \t\n]+[ \t]*:;?\\(.*\\)$" 1 font-lock-doc-string-face t)
;   ;; variable definition
;   ("^[_A-Za-z0-9]+[ \t]*\+?=" . font-lock-function-name-face)
;   ("\\( \\|:=\\)[_A-Za-z0-9]+[ \t]*\\+=" . font-lock-function-name-face)
;   ;; variable references
;   ("\\(\\$\\$?\\([^ \t\n{(]\\|[{(][^ \t\n)}]+[)}]\\)\\)" 
;	1 font-lock-keyword-face t)
;   ("^include " . font-lock-string-face))
;  "Keywords for fontifying makefiles"
;)

;; GMF: Took the above an modified to the following.  Seems to be more
;; regular for the kinds of makefiles I write.

(defconst makefile-font-lock-keywords '(
   ("^#.*$" . font-lock-comment-face)
   ("[^$]#.*$" . font-lock-comment-face)
   ;; rules
   ("^\\([^ \t\n]*%?[^ \t\n]*[ \t]*::?\\)[ \t\n]" 1 font-lock-type-face t)
   ("^\\(\\.[A-Za-z][A-Za-z]?\\..[ \t]*::?\\)" 1 font-lock-type-face t)
   ("^[^ \t\n]+[ \t]*:;?\\(.*\\)$" 1 font-lock-doc-string-face t)
   ;; variable definition
;   ("^[_A-Za-z0-9+]+[ \t]*\\+?=" . font-lock-function-name-face)
;   ("\\( \\|:=\\)[_A-Za-z0-9]+[ \t]*\\+=" . font-lock-function-name-face)
   ("^[_A-Za-z0-9+]+[ \t]*\\+?=" . font-lock-doc-string-face)
   ("\\( \\|:=\\)[_A-Za-z0-9]+[ \t]*\\+=" . font-lock-doc-string-face)
   ;; variable references
;   ("\\(\\$\\$?\\([^ \t\n{(]\\|[{(][^ \t\n)}]+[)}]\\)\\)" 
;	1 font-lock-keyword-face t)
   ("\\(\\$\\$?\\([^ \t\n{(]\\|[{(][^ \t\n)}]+[)}]\\)\\)" 
	1 font-lock-type-face t)
;;   ("^include " . font-lock-string-face)

;; GNU make special stuff.

;    ("^\\(include\\|ifdef\\|ifndef\\|ifeq\\|ifneq\\|else\\|endif\\)[ 	\n]+" 1 font-lock-string-face)
    ("^\\(include\\|ifdef\\|ifndef\\|ifeq\\|ifneq\\|else\\|endif\\)[ 	\n]+" 1 font-lock-keyword-face)
;    ("^ifdef[ 	]+\\([A-Za-z0-9_+]+\\)" 1 font-lock-function-name-face)
    ("^ifdef[ 	]+\\([A-Za-z0-9_+]+\\)" 1 font-lock-doc-string-face)
;    ("\\$(\\(wildcard\\|patsubst\\|shell\\)" 1 font-lock-doc-string-face)
    ("\\$(\\(wildcard\\|patsubst\\|shell\\)[ 	]" 1 font-lock-function-name-face)
   )
  "Keywords for fontifying makefiles"
)

(defconst pascal-font-lock-keywords '(
  ("(\\*\\(.*\\)\\*)" 1 font-lock-comment-face t)
   ("{\\([^}]*\\)}" 1 font-lock-comment-face t)
   ;; Doesn't work when there are strings in comments....
   ;; ("'[^']*'" nil string)
   ("^#.*$" font-lock-doc-string-face)
   ("^[ \t]*\\(procedure\\|function\\)[ \t]+\\w+[^ \t(;]*" font-lock-function-name-face)
   ("\\<\\(program\\|begin\\|end\\)\\>" bold-italic)
   ("\\<\\(external\\|forward\\)\\>" bold)
   ("\\<\\(label\\|const\\|type\\|var\\)\\>" font-lock-keyword-face)
   ("\\<\\(record\\|array\\|file\\)\\>"  font-lock-type-face)
   ("\\<\\(of\\|to\\|for\\|if\\|then\\|else\\|case\\|while\\|do\\|until\\|and\\|or\\|not\\|with\\|repeat\\)\\>" font-lock-keyword-face)
   )
  "Keywords for fontifying pascal programs"
)

(make-face 'dired-permissions-face)
(make-face 'dired-links-face)
(make-face 'dired-owner-name-face)
(make-face 'dired-group-name-face)
(make-face 'dired-size-face)
(make-face 'dired-date-time-face)

(make-face 'dired-default-filename-face)
(make-face 'dired-directory-name-face)

(defconst dired-font-lock-keywords-2
  '(
    ("^  [dl]?\\([-rwx]+\\)" 1 dired-permissions-face)
    ("^\\s +[^ \\n]+\\s +\\([0-9]+\\)" 1 dired-links-face)
    ("^  [-ldrwx]+\\s +[0-9]+\\s +\\([A-Za-z_]+\\)" 1 dired-owner-name-face)
    ("^  [-ldrwx]+\\s +[0-9]+\\s +[A-Za-z_]+\\s +\\([A-Za-z_]+\\)" 1 dired-group-name-face)
    ("^\\s +[^ \\n]+\\s +[^ ]+\\s +[^ ]+\\s +[^ ]+\\s +\\([0-9]+\\)" 1 dired-size-face)
    ("^\\s +[^ \\n]+\\s +[^ ]+\\s +[^ ]+\\s +[^ ]+\\s +[0-9]+\\s +\\([A-Za-z]+\\s +[0-9]+\\s +[0-9:]+\\)" 1 dired-date-time-face)
    ("^\\s +[^ \\n]+\\s +[^ ]+\\s +[^ ]+\\s +[^ ]+\\s +[0-9]+\\s +\\([A-Za-z]+\\s +[0-9]+\\s +[0-9:]+\\)\\s +\\(.*\\)" 2 dired-default-filename-face)
    )
  "Keywords for fontifying dired buffers"
  )

(defconst tcl-font-lock-keywords
  '(
; Plain Tcl keywords.
    ("[ 	:[{]\\(tell\\|open\\|eof\\|pwd\\|glob\\|list\\|pid\\|exec\\|time\\|unknown\\|eval\\|lrange\\|lsearch\\|gets\\|case\\|lappend\\|proc\\|break\\|llength\\|auto_execok\\|return\\|linsert\\|error\\|catch\\|info\\|split\\|array\\|if\\|auto_mkindex\\|concat\\|join\\|lreplace\\|source\\|global\\|switch\\|auto_reset\\|close\\|for\\|cd\\|auto_load\\|file\\|append\\|format\\|read\\|set\\|scan\\|trace\\|seek\\|while\\|flush\\|continue\\|uplevel\\|foreach\\|rename\\|regexp\\|upvar\\|expr\\|unset\\|regsub\\|history\\|exit\\|puts\\|incr\\|lindex\\|lsort\\|string\\)[ 	:]" 1 font-lock-keyword-face)
; Plain Tk widget names.
    ("[ 	:]\\(frame\\|button\\|label\\|menu\\|menubutton\\|checkbutton\\|radiobutton\\|scrollbar\\|scale\\|toplevel\\)[ 	]" 
     1 font-lock-type-face)
    ("\\<pack\\|blt_table\\>" . font-lock-keyword-face)
    ("-\\(textvariable\\|command\\|value\\|variable\\|text\\)" 1 font-lock-comment-face)
    ("-\\(showvalue\\|from\\|to\\|orient\\|label\\)" 1 font-lock-comment-face)
    ("-\\(side\\|expand\\|fill\\|after\\|before\\|anchor\\)" 1 font-lock-comment-face)
; [incr Tcl] keywords.
    ("\\<inherit\\|itcl_class\\|constructor\\|destructor\\|method\\|public\\|protected\\>" 
     . font-lock-keyword-face)
; [Incr Tcl] misc.
    ("inherit \\([A-Za-z_]+\\)" 1 font-lock-type-face)
    ("itcl_class \\([A-Za-z_]+\\)" 1 font-lock-type-face)
    ("method \\([A-Za-z_]+\\)" 1 font-lock-function-name-face)
; Variable substitutions
    ("$\\([A-Za-z0-9_\\.]+\\)" 1 font-lock-doc-string-face)
; Add new regexps here.
    )
  "Keywords for fontifying Tcl files"
  )


(setq c-font-lock-keywords c-font-lock-keywords-2)
(setq c++-font-lock-keywords c++-font-lock-keywords-gmf)
(setq lisp-font-lock-keywords lisp-font-lock-keywords-2)
(setq dired-font-lock-keywords dired-font-lock-keywords-2)
(setq java-font-lock-keywords java-font-lock-keywords-3)

(provide 'fl-keywords)
