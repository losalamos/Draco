;; ======================================================================
;; draco-global-keys.el
;; Kelly Thompson
;; 8 Dec 2004
;;
;; $Id$
;;
;; Define some global key bindings
;; 
;; Usage: (require 'draco-global-keys)
;; ======================================================================

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

(if (string-match "lambda" machine-name)
    (progn
      (define-key global-map [(button4)] 
        '(lambda () (interactive) (scroll-down 5)))
      (define-key global-map [(button5)] 
        '(lambda () (interactive) (scroll-up 5)))))

;(if (or (string-match "ffe1" machine-name)
;	(string-match "qscfe1" machine-name)) (mwheel-install))

(mwheel-install)

(provide 'draco-global-keys)
