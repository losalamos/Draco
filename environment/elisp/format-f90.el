;;; File  : format-f90.el
;;; Author: Kelly Thompson, kgt@lanl.gov
;;; Date  : Wednesday, Sep 05, 2018, 16:01 pm

;;; Description:

;;; With format-f90.sh, load a file in Emacs, apply standard Draco
;;; formatting rules, save and exit Emacs via a non-interactive batch
;;; process.

;;; Use with format-f90.sh.
;;;----------------------------------------------------------------------------

(defun collapse-blank-lines
  (start end)
  (interactive "r")
  (replace-regexp "^\n\\{2,\\}" "\n" nil start end)
)

(defun emacs-format-f90-sources ()
   "Format the whole buffer accourding to Draco F90 rules."
   ;; For more information about f90-mode settins in Emacs:
   ;; 1. M-x f90-mode
   ;; 2. C-h m
   (set-fill-column 80)
   (turn-on-auto-fill)
   (setq f90-beginning-ampersand nil)
   (setq f90-associate-indent 0)
   ;; (setq f90-do-indent 3)
   ;; (setq f90-if-indent 3)
   ;; (setq f90-type-indent 3)
   ;; (setq f90-program-indent 2)
   ;; (setq f90-critical-indent 2)
   ;; (setq f90-continuation-indent 5)
   ;; (setq f90-indented-comment-re "!")
   ;; (setq f90-break-delimiters "[-+*/><=,% \t]")
   ;; (setq f90-break-before-delimiters t)
   ;; (setq f90-auto-keyword-case 'downcase-word)
   (indent-region (point-min) (point-max) nil)
   (untabify (point-min) (point-max))
   (collapse-blank-lines (point-min) (point-max))
   (save-buffer)
)
