;; $Id$
;; 
;; Add additional keywords for font-lock modes.
;;
;;======================================================================

(require 'font-lock)

(defun add-draco-dbc-font-lock-keywords ()
  (interactive)
  (setq font-lock-keywords
	(append
	 '(("\\<\\(Check\\|Ensure\\|Insist\\|Require\\)\\>"
	    1 font-lock-warning-face t))
	 c++-font-lock-keywords-3)))

(provide 'fl-keywords)

