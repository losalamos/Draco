;; ======================================================================
;; draco-global-keys.el
;; Kelly Thompson
;; 8 Dec 2004
;;
;; Define some global key bindings
;;
;; Usage: (require 'draco-global-keys)
;; ======================================================================

;; Clang-format
(if (file-exists-p (concat (getenv "VENDOR_DIR") "/bin/clang-format" ))
    (if (> emacs-major-version 23)
        (progn
          (require 'clang-format)
          (global-set-key [(f12)] 'clang-format-region) ;; Windows/Linux
          (global-set-key [(C-M-tab)] 'clang-format-region) ;; Mac/Linux
          )))

;; Kill default XEmacs key binding annoyances:

(define-key global-map "\C-x\C-k" 'kill-buffer)
(define-key global-map "\e?"      'help-for-help)
(define-key global-map "\C-x?"    'describe-key-briefly)
(define-key global-map "\C-z"     'undo)

;; Force these onto the global map as well as the draco-mode-map
(define-key global-map [(button3)]          'kill-region)
(define-key global-map [(meta button3)]     'delete-rectangle)
(define-key global-map [(control button1)]  'popup-buffer-menu)
(define-key global-map [(control button2)]  'function-menu)
(define-key global-map [(control button3)]  'popup-mode-menu)

;; Comments
(define-key global-map [(f2)]               'comment-region)

;; Search
(define-key global-map [(f3)]               'isearch-forward)
(define-key global-map [(shift f3)]         'isearch-backward)

; Refresh
(define-key global-map [(f5)]               'font-lock-fontify-buffer)

;; Buffer management
(define-key global-map [(f7)]               'draco-save-and-kill-current-buffer)
(define-key global-map [(shift f7)]         'delete-frame)  ;; delete-window might work better
(define-key global-map [(control f7)]       'kill-this-buffer)
(define-key global-map [(f8)]               'draco-toggle-previous-buffer)
(define-key global-map [(control f8)]       'draco-find-companion-file)

;; Newer keyboards (Microsoft 4000 via Reflection X)
(define-key global-map [(f17)]       'kill-this-buffer)
(define-key global-map [(f18)]       'draco-find-companion-file)

;; Mouse wheel for Lambda, Flash and QSC
(if (getenv "HOST")
    (setq machine-name (getenv "HOST"))
(setq machine-name "none"))

;; (if (string-match "lambda" machine-name)
;;     (progn
;;       (define-key global-map [(button4)]
;;         '(lambda () (interactive) (scroll-down 5)))
;;       (define-key global-map [(button5)]
;;         '(lambda () (interactive) (scroll-up 5)))))

;(if (or (string-match "ffe1" machine-name)
;	(string-match "qscfe1" machine-name)) (mwheel-install))

(mwheel-install)

;; Ack/Grep
;; http://www.emacswiki.org/Ack

(defvar ack-history nil
  "History for the `ack' command.")
(defun ack (command-args)
  (interactive
   (let ((ack-command
          ;; "ack --nocolor --nogroup --with-filename --all "))
          "ack --nocolor --nogroup --with-filename "))
     (list (read-shell-command "Run ack (like this): "
                               ack-command
                               'ack-history))))
  (let ((compilation-disable-input t))
    (compilation-start (concat command-args " < " null-device)
                       'grep-mode)))

(define-key draco-mode-map [(shift f11)] 'grep)
(define-key draco-mode-map [(f11)] 'ack)

(provide 'draco-global-keys)
