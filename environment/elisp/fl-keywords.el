;; $Id$
;; 
;; Add additional keywords for font-lock modes.
;;
;;======================================================================

(require 'font-lock)

(make-face 'font-lock-kull-macros-face)

(defun add-draco-dbc-font-lock-keywords ()
  (interactive)
  (setq font-lock-keywords
	(append
	 '((
            "\\<\\(Check\\|Ensure\\|Insist\\|Assert\\|Remember\\|Require\\)\\>\\([ ]*\\s(\\)"
	    1 font-lock-kull-macros-face t))
	 c++-font-lock-keywords-3)))

(provide 'fl-keywords)

