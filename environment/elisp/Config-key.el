;;==============================================================================
;; Config-key.el
;; Geoffrey Furnish
;; 23 June 1997
;;
;; $Id$
;;
;; Customized keymaps
;;
;;==============================================================================

;;
;; Additional key commands for global map
;;

;; XEmacs key binding annoyances:

; Disable C-x f command key (C-x C-f is open file and may be confused
; with C-x f)
(define-key global-map "\C-xf" nil) 

; Assign C-x C-k to (kill-buffer).  C-x k also calls this function.
; When killing lots of buffers it is easier to just hold down the C
; key.
(define-key global-map "\C-x\C-k" 'kill-buffer)

;
; Fix the delete key on some PC keyboards - option 1
;
;(define-key global-map "\C-h" 'delete-backward-char)

;
; Fix the delete key on some PC keyboards - option 2
;
;(defun forward-delete-char ()
;  "Delete character to the right of the point."
;  (interactive)
;  (backward-delete-char -1))

;(define-key global-map '(delete) 'forward-delete-char)      ;; delete char to the right
;(define-key global-map '(backspace) 'delete-backward-char)  ;; delete char to the left
;(global-set-key [?\C-h] 'delete-backward-char)
;(global-set-key [?\C-x ?h] 'help-command)

;; RTT Customizations to the global keymap

;; misc

(define-key global-map "\e?"     'help-for-help)
(define-key global-map "\C-x?"   'describe-key-briefly)
(define-key global-map [(f11)]   'grep)
(define-key global-map [(control f11)] 'speedbar)
(define-key global-map "\C-csb"  'speedbar)

;; CVS

(define-key global-map [(f9)]         'cvs-examine)
(define-key global-map [(meta f9)]    'cvs-status)
(define-key global-map [(control f9)] 'cvs-checkout)

;; Mouse commands

(define-key global-map 'button3            'kill-region)
(define-key global-map [(meta button3)]    'delete-rectangle)
(define-key global-map [(control button1)] 'popup-buffer-menu)
(define-key global-map [(control button2)] 'popup-mode-menu)
(define-key global-map [(control button3)] 'popup-mode-menu)

;; Use the keypad commands:

;(define-key global-map [(kp-divide)] 	     'byte-compile-file)
(define-key global-map [(kp-multiply)]       'start-kbd-macro)
(define-key global-map [(kp-subtract)]       'end-kbd-macro)
(define-key global-map [(kp-add)] 	     'call-last-kbd-macro)
;(define-key global-map [(kp-enter)] 	     'other-window)
;(define-key global-map [(control kp-enter)] 'gmf-top-other-window)
;(define-key global-map [(shift kp-enter)]   'gmf-bot-other-window)
;(define-key global-map [(kp-7)] 	     'font-lock-fontify-buffer)
;(define-key global-map [(delete)] 	     'delete-char)
;(define-key global-map [(hpDeleteChar)]     'delete-char)
;(define-key global-map [(shift   delete)]   'delete-char)
;(define-key global-map [(control delete)]   'delete-char)
;(define-key global-map [(delete)]           'delete-backward-char)
;(define-key global-map [(prior)]            'scroll-down)
;(define-key global-map [(next)]             'scroll-up)

;; Compiling

(define-key global-map [(f1)]              'compile)
(define-key global-map [(control meta g)]  'what-line)
(define-key global-map [(f3)]              'previous-error)
(define-key global-map [(f4)]              'next-error)

;; Inserting comment blocks

(define-key global-map [(f5)]               'rtt-insert-function-doc)
(define-key global-map [(meta control f5)]  'rtt-insert-class-doc)
(define-key global-map [(f6)]               'rtt-insert-comment-divider)

;; Insert ChangeLog comment.
(define-key global-map [(shift control f5)]   'add-change-log-entry)

;; Buffer management

(define-key global-map [(f7)]               'rtt-save-and-kill)
(define-key global-map [(shift f7)]         'delete-window)
(define-key global-map [(control f7)]       'kill-this-buffer)
(define-key global-map [(f8)]               'rtt-previous-buffer)
(define-key global-map [(control f8)]       'rtt-find-companion-file)

; Some bindings to drop Brief style bookmarks.
; Shift Fx sets bookmakr x, Control Fx returns to bookmark x.

(define-key global-map [(meta f1)]     [(control x) r (space) a])
(define-key global-map [(control f1)]  [(control x) r j a])
(define-key global-map [(meta f2)]     [(control x) r (space) b])
(define-key global-map [(control f2)]  [(control x) r j b])
(define-key global-map [(meta f3)]     [(control x) r (space) c])
(define-key global-map [(control f3)]  [(control x) r j c])
(define-key global-map [(meta f4)]     [(control x) r (space) d])
(define-key global-map [(control f4)]  [(control x) r j d])

;; Dired 

(require 'dired)
(define-key dired-mode-map 'f5 'dired-redisplay-subdir)
