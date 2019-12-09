;; Copy this file to $HOME/.emacs and edit as needed.

;; Ref: www.emacswiki.org
;;      www.linuxpowered.com/html/tutorials/emacs.html


;; FAQ
;;
;; - increase font size: C-x X-+

;; Use:
;; 1. byte compile draco elisp: (cd ~/draco/environment/elisp; ./compile-elisp)
;; 2. byte compile this file.  From emacs: M-x byte-compile-file
;; 3. Set draco environment directory:
(defvar my-draco-env-dir "~/draco/environment/")
;; 4. Update customizations found in this file

;;------------------------------------------------------------------------------

;;
;; System information
;;
(defvar emacs-load-start-time (current-time))
(defconst linuxp
    (or (eq system-type 'gnu/linux)
        (eq system-type 'linux)) "Are we running on a GNU/Linux system?")
(defconst linux-x-p
    (and window-system 'linuxp) "Are we running under X on a GNU/Linux system?")
(defconst emacs>=24p (> emacs-major-version 23)
  "Are we running GNU Emacs 24 or above?")
(defvar emacs-debug-loading nil)

;;
;; Draco Environment
;;
(load-library (concat my-draco-env-dir "elisp/draco-setup"))

;;---------------------------------------------------------------------------------------
;; Personal Settings below this line
;;---------------------------------------------------------------------------------------

;; this is the same title bar format all Gnome apps seem to use
(setq frame-title-format (list "[" (system-name) "] %b" ))

; Turn off warning:
(setq large-file-warning-threshold nil)

;; Automatically keep buffers synchronized with file system
;; (i.e. reload the buffer automatically if the file changes outside
;; of emacs).
(global-auto-revert-mode t)

;; CEDET
;; Generate tag files for Fortran: 'gtags --gtagslabel=ctags'
(defvar file-at-root-level-draco "~/.draco_ede")
(defvar file-at-root-level-eap "~/cassio/.eap_ede")
; (draco-start-ecb)

;; Allow 'emacsclient' to connect to running emacs.
(server-start)
(add-hook 'server-switch-hook
          (lambda ()
            (when (current-local-map)
              (use-local-map (copy-keymap (current-local-map))))
            (when server-buffer-clients
              (local-set-key (kbd "C-c C-c") 'server-edit))))

;; ========================================
;; GNU Emacs Custom settings
;; ========================================

(custom-set-variables
 ;; custom-set-variables was added by Custom.
 ;; If you edit it by hand, you could mess it up, so be careful.
 ;; Your init file should contain only one such instance.
 ;; If there is more than one, they won't work right.
 '(column-number-mode t)
 '(compilation-auto-jump-to-first-error nil)
 '(compilation-scroll-output (quote first-error))
 '(cua-mode t nil (cua-base))
 '(draco-code-comment-width 80)
 '(draco-elisp-dir (concat my-draco-env-dir "elisp/"))
 '(draco-env-dir my-draco-env-dir)
; '(ecb-layout-window-sizes (quote (("ecb-layout-kt" (0.25 . 0.25) (0.15 . 0.25) (0.4 . 0.45) (0.4 . 0.3)))))
; '(ecb-options-version "2.40")
; '(ecb-prescan-directories-for-emptyness nil)
; '(ecb-primary-secondary-mouse-buttons (quote mouse-1--mouse-2))
; '(ecb-vc-supported-backends (quote ((ecb-vc-dir-managed-by-SVN . ecb-vc-state) (ecb-vc-dir-managed-by-CVS . ecb-vc-state) (ecb-vc-dir-managed-by-GIT . ecb-vc-state) (ecb-vc-dir-managed-by-MTN . ecb-vc-state))))
 '(font-lock-maximum-decoration t)
 '(global-font-lock-mode t nil (font-lock))
 '(inhibit-startup-screen t)
 '(scroll-bar-mode (quote right))
 '(show-paren-mode t nil (paren))
 '(text-mode-hook (quote (turn-on-auto-fill text-mode-hook-identify)))
 '(tool-bar-mode nil)
 '(uniquify-buffer-name-style (quote forward) nil (uniquify)))

(custom-set-faces
 ;; custom-set-faces was added by Custom.
 ;; If you edit it by hand, you could mess it up, so be careful.
 ;; Your init file should contain only one such instance.
 ;; If there is more than one, they won't work right.
 '(default ((t (:family "Inconsolata" :foundry "unknown" :slant normal
                        :weight normal :height 90 :width normal
                        :background "mint cream" ))))
 '(font-lock-comment-face ((((class color) (min-colors 88) (background light)) (:foreground "purple"))))
 '(font-lock-constant-face ((((class color) (min-colors 88) (background light)) (:foreground "CadetBlue"))))
 '(font-lock-doc-face ((t (:inherit font-lock-string-face :foreground "bisque3"))))
 '(font-lock-function-name-face ((((class color) (min-colors 88) (background light)) (:foreground "Blue3"))))
 '(font-lock-keyword-face ((((class color) (min-colors 88) (background light)) (:foreground "firebrick2"))))
 '(font-lock-preprocessor-face ((t (:inherit font-lock-builtin-face :foreground "hotpink"))))
 '(font-lock-string-face ((((class color) (min-colors 88) (background light)) (:foreground "orange3"))))
 '(font-lock-variable-name-face ((((class color) (min-colors 88) (background light)) (:foreground "Royalblue"))))
 '(menu ((((type x-toolkit)) (:weight bold :height 0.8 :family "helvetica"))))
 '(mode-line ((t (:background "wheat" :foreground "black" :box (:line-width -1 :style released-button) :weight bold :height 0.9 :width normal :foundry "inconsolata" :family "unknown"))))
 '(mode-line-emphasis ((t (:foreground "red" :weight bold))))
 '(mode-line-inactive ((t (:background "wheat" :family "Inconsolata" :foundry "unknown" :slant normal :weight normal :height 90 :width normal))))
 '(modeline ((t (:family "Inconsolata" :foundry "unknown" :slant normal :weight normal :height 90 :width normal))))
 '(modeline-face ((((class color) (min-colors 88) (background light)) (:background "wheat")))))
(setq x-alt-keysym 'meta)
