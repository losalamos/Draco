; -*- Mode: emacs -*-
; ------------------------------------------------------------------------------
; File:   f90-format.el
; Date:   Wednesday, Nov 07, 2018, 11:11 am
; Author: Kelly Thompson
; Note:   Copyright (C) 2016-2018, Los Alamos National Security, LLC.
;         All rights are reserved.
; ------------------------------------------------------------------------------
; If emacs is available, enforce f90 formatting (indentation, etc) rules.
; ------------------------------------------------------------------------------

(defun collapse-blank-lines
  (start end)
  (interactive "r")
  (replace-regexp "^\n\\{2,\\}" "\n" nil start end)
)

; C-u 81 M-x goto-long-line
; C-e
; <space>
; [repeat]
(defun goto-long-line (len)
  "Go to the first line that is at least LEN characters long.
Use a prefix arg to provide LEN.
Plain `C-u' (no number) uses `fill-column' as LEN."
  (interactive "P")
  (setq len  (if (consp len) fill-column (prefix-numeric-value len)))
  (let ((start-line                 (line-number-at-pos))
        (len-found                  0)
        (found                      nil)
        (inhibit-field-text-motion  t))
    (while (and (not found) (not (eobp)))
      (forward-line 1)
      (setq found  (< len (setq len-found  (- (line-end-position) (point))))))
    (if found
        (when (interactive-p)
          (message "Line %d: %d chars" (line-number-at-pos) len-found))
      (goto-line start-line)
      (message "Not found"))))

;; (defun wrap-long-lines ()
;;   "Find all lines > 80 columns in the current buffer and wrap them."
;;   (interactive "p")
;;   (set-fill-column 80)
;;   (goto-long-line)
;; foob
;;   (end-of-line)
;;   (f90-do-auto-fill)
;; )

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
   (untabify (point-min) (point-max))
   (delete-trailing-whitespace)
   (indent-region (point-min) (point-max) nil)
   (collapse-blank-lines (point-min) (point-max))
   (save-buffer)
)

; ------------------------------------------------------------------------------
; end f90-format.el
; ------------------------------------------------------------------------------
