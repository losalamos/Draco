;; infer-cc-style.el
;; Geoffrey Furnish
;; 9 July 1997
;;
;; Facilitate work on multiple C++ projects with differing editing
;; styles, by providing a way to set the cc-mode style based on
;; inferences drawn from studying the contents of the buffer.  This is
;; based directly on the infer-file-mode code.
;;
;; Probably we need to generalize this to support heuristics based
;; both on file contents, and on file name.  But that will have to
;; wait until later.

(require 'pooma-hacks)

(defvar infer-cc-style-alist
  '(
    ("The POOMA Framework" . "POOMA")
    )
  "Alist of file content patterns and corresponding cc-mode style names.")

(defvar infer-cc-style-limit 1000
  "Number of characters to search when inferring cc-mode style based
on file contents.  Nil means search entire buffer.  See
infer-cc-style-alist.")

(defun infer-cc-style ()
  "Infer cc-mode style based on contents of buffer.  To use, add this
function to your c-mode-common-hook, or such like.  To customize, edit
infer-cc-style-alist."
  (let ((style-alist infer-cc-style-alist))
    (catch 'found
      (while style-alist
	(if (save-excursion
	      (goto-char (point-min))
	      (re-search-forward (car (car style-alist)) 
				 infer-cc-style-limit t))
	    (throw 'found t)
	  (setq style-alist (cdr style-alist)))))
    (if style-alist
	(c-set-style (cdr (car style-alist))))))

(provide 'infer-cc-style)
