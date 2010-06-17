;; Copy this file to $HOME/.emacs and edit as needed.

;;
;; Initial personal stuff
;;
(defvar emacs-load-start-time (current-time))

(defconst win32p
    (eq system-type 'windows-nt)
  "Are we running on a WinTel system?")

(defconst cygwinp
    (eq system-type 'cygwin)
  "Are we running on a WinTel cygwin system?")

(defconst linuxp
    (or (eq system-type 'gnu/linux)
        (eq system-type 'linux))
  "Are we running on a GNU/Linux system?")

(defconst unixp
  (or linuxp
      (eq system-type 'usg-unix-v)
      (eq system-type 'berkeley-unix))
  "Are we running unix")

(defconst linux-x-p
    (and window-system 'linuxp)
  "Are we running under X on a GNU/Linux system?")

(defconst xemacsp (featurep 'xemacs)
  "Are we running XEmacs?")

(defconst emacs>=21p 
  (and (not xemacsp) (or (= emacs-major-version 21)
                         (= emacs-major-version 22)))
  "Are we running GNU Emacs 21 or above?")

(defvar emacs-debug-loading t)

;; ================================================================================
;; CCS-4 default settings
;; ================================================================================

(defun run-command (command args)
  "Return the result of running a system command."
  (replace-regexp-in-string
   "[ \n]+$" "" (shell-command-to-string (concat command " " args))))

(defconst os-name (run-command "uname" nil))
(defconst machine-name (system-name))

(if (string-match "Linux" os-name)
    (progn
      (setq draco-env-dirs (list
                            "~/draco/environment/"
                            "~/.xemacs" ))
      (load-library "~/draco/environment/elisp/draco-setup")
      ))

(when emacs-debug-loading
  (defadvice load (before debug-log activate)
    (message "Now loading: %s" (ad-get-arg 0))))

(if emacs>=21p
    ;; If this is GNU Emacs, start the server
    (server-start)
  (progn
    ;; If this is XEmacs, start the server
    (load "gnuserv")
    (gnuserv-start))
  )

; Display function name in buffer
(if emacs>=21p
    (which-function-mode t))

;;---------------------------------------------------------------------------------------
;; Personal Settings below this line
;;---------------------------------------------------------------------------------------


(custom-set-variables
  ;; custom-set-variables was added by Custom.
  ;; If you edit it by hand, you could mess it up, so be careful.
  ;; Your init file should contain only one such instance.
  ;; If there is more than one, they won't work right.
 '(inhibit-startup-screen t)
 '(scroll-bar-mode (quote right))
 '(show-paren-mode t)
 '(tool-bar-mode nil)
 '(transient-mark-mode t))

(custom-set-faces
  ;; custom-set-faces was added by Custom.
  ;; If you edit it by hand, you could mess it up, so be careful.
  ;; Your init file should contain only one such instance.
  ;; If there is more than one, they won't work right.
 '(default ((t (:stipple nil :background "white" :foreground "black" :inverse-video nil :box nil :strike-through nil :overline nil :underline nil :slant normal :weight normal :height 98 :width semi-condensed :family "misc-fixed")))))


(unless xemacsp (set-language-environment 'Latin-1))
;(set-keyboard-coding-system 'iso-8859-1) ??

;;; Default Settings
(setq next-line-add-newlines nil)
(setq track-eol nil)
(setq scroll-step 1)
;;(setq scroll-conservatively 10000)

(setq hscroll-step 1)
(setq make-backup-files nil)
(line-number-mode 1)     ; line-numbers
(column-number-mode 1)
(setq visible-bell t) ; no beeping
(when (fboundp 'blink-cursor-mode) (blink-cursor-mode -1)) ; no blinking cursor
(setq default-tab-width 4)
(setq-default indent-tabs-mode nil)
(setq imenu-max-items 40)
(setq message-log-max 3000)
;(setq echo-keystrokes 0.1)
(setq history-length 100)
(setq line-number-display-limit 10000000)
(setq sentence-end-double-space nil)
(setq read-quoted-char-radix 10) ; accept decimal input when using ^q, e.g.: ^q 13 [RET] -> ^M
(setq yank-excluded-properties t) ; do not paste any properties
(setq confirm-nonexistent-file-or-buffer t)
;;(setq directory-sep-char ?\\)
(setq completion-ignored-extensions (remove ".pdf" completion-ignored-extensions))
(setq completion-ignored-extensions (remove ".dvi" completion-ignored-extensions))

(setq max-specpdl-size 32000)
(setq max-lisp-eval-depth 32000)

(setq apropos-sort-by-scores t)

(setq-default case-fold-search t)

(require 'uniquify)
(setq uniquify-non-file-buffer-names t)
(setq uniquify-after-kill-buffer-p t)
(setq uniquify-buffer-name-style 'post-forward-angle-brackets)
(setq uniquify-ignore-buffers-re "\\`\\*")
;(toggle-uniquify-buffer-names)


;;; Window System specific code
(cond (window-system
       ;; use some nicer colors for font-lock mode
       ;=(set-face-foreground 'font-lock-comment-face "gray50")
       ;=(set-face-foreground 'font-lock-string-face "green4")

       (unless xemacsp (global-font-lock-mode t))
       (setq font-lock-maximum-decoration t)
       (setq lazy-lock-defer-on-scrolling t)
       (if emacs>=21p
           (progn
             (setq font-lock-support-mode 'jit-lock-mode)
             ;; (setq jit-lock-stealth-time 16)
             ;; (setq jit-lock-stealth-nice 0.5)
             (setq font-lock-multiline t))
         (setq font-lock-support-mode 'lazy-lock-mode))
       ;(setq lazy-lock-defer-contextually t)
       ;(setq lazy-lock-defer-time 0)
))

;; http://www.emacswiki.org/emacs/AddKeywords
(unless (fboundp 'font-lock-add-keywords)
  (defalias 'font-lock-add-keywords 'ignore))

