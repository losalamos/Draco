;; $Id$
;; 
;; Add additional keywords for font-lock modes.
;;
;;======================================================================

(require 'font-lock)

(make-face 'font-lock-kull-macros-face)

(defconst c++-font-lock-keywords-gmf
  '(
    ("\\<\\(Check\\|Require\\|Ensure\\|Remember\\|Assert\\|Insist\\)\\>" 1 font-lock-kull-macros-face)
    ))

(setq c++-font-lock-keywords (cons c++-font-lock-keywords c++-font-lock-keywords-gmf))

(provide 'fl-keywords)
