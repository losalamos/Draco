; nml-mode.el

(defconst nml-font-lock-keywords
  '(
    ("\\(block\\):" 1 font-lock-keyword-face)
    ("\\(boxstart\\|boxend\\|page\\)" 1 font-lock-keyword-face)
    ("block:[ 	]*\\([A-Za-z_]+\\)" 1 font-lock-function-name-face)
    ("\\(#.*\\)" 1 font-lock-comment-face)
    ("^[ 	]*\\([A-Za-z_]+\\)[ 	]+\\([A-Za-z0-9_]+\\)[ 	]+\\([A-Za-z0-9_.-]+\\)" 1 font-lock-type-face)
    ("\"\\(.*\\)\"" 1 font-lock-string-face)
    )
  "Regexps for showing off namelist description files."
  )

(defun nml-mode ()
  "A major mode for editing Namelist Description files."
  (interactive)
  (setq major-mode 'nml-mode)
  (setq mode-name  "nml")
  (run-hooks 'nml-mode-hook))
