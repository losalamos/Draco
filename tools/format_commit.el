;; Lifted from draco-config.mode.el

(autoload 'c++-mode "cc-mode" "C++ Editing Mode" t)
(autoload 'c-mode   "cc-mode" "C Editing Mode" t)

(setq auto-mode-alist
      (append '(("\\.C$"      . c++-mode)
                ("\\.cc$"     . c++-mode)
                ("\\.cxx$"    . c++-mode)
                ("\\.pt$"     . c++-mode)
                ("\\.hh$"     . c++-mode)
                ("\\.hpp$"    . c++-mode)
                ("\\.cpp$"    . c++-mode)
                ("\\.hh.in$"  . c++-mode)
                ("\\.h.in$"   . c-mode)
                ("\\.c$"      . c-mode)   ; to edit C code
                ("\\.cu$"     . c-mode)   ; to edit CUDA kernels
                ("\\.h$"      . c-mode)   ; to edit C code
                ("\\.dcc$"    . c-mode)   ; to edit C code
                ("\\.dcc.in$" . c-mode)   ; to edit C code
                ("\\.dot$"    . c-mode)  ; for dot files
                ) auto-mode-alist))

(c-add-style
 "draco" '
 (
  (c-basic-offset . 4)
  (c-electric-pound-behavior . 'alignleft)
  (c-offsets-alist . (
                      (access-label . -2 )
                      (block-close . 0)
                      (block-open  . 0)
                      (case-label  . +)
                      (class-close . 0)
                      (class-open  . 0)
                      (defun-close . 0)
                      (defun-open  . 0)
                      (do-while-closure  . 0)
                      (else-clause       . 0)
                      ;;(extern-lang-close . +)
                      ;;(extern-lang-open  . +)
                      (extern-lang-close . 0)
                      (extern-lang-open  . 0)
                      (inline-close      . 0)
                      (inline-open       . 0)
                      (innamespace       . 0)
                      (statement-case-intro . +)
                      (statement-cont    . c-lineup-math)
                      (substatement-open . 0)
                      ))))

;; ================================================================================
;; Actions on the buffer
(c-set-style "draco")
(set-fill-column 80)
(untabify (point-min) (point-max))
(indent-region (point-min) (point-max))
(add-hook 'before-save-hook 'delete-trailing-whitespace)
(set-buffer-file-coding-system 'iso-latin-1-unix t)
