;; $Id$
;; 
;; Add additional keywords for font-lock modes.
;;
;;======================================================================

(require 'font-lock)

(defun add-draco-dbc-font-lock-keywords ()
  (setq font-lock-keywords
	(append
	 '(("\\<\\(Check\\|Ensure\\|Insist\\|Require\\)\\>"
	    0 font-lock-warning-face t))
	 c++-font-lock-keywords-2)))

(provide 'fl-keywords)

