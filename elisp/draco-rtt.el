;;; draco-rtt.el
;;; Tom Evans
;;; Feb 24, 1999

;;; $Id$

;;---------------------------------------------------------------------------;;
;; Load the DRACO-RTT Elisp packages
;; ---------------------------------
;;
;; This file should be called from the user's .emacs file in the form:
;; 
;; (setq my-home-dir (getenv "HOME"))
;; (setq my-elisp-dir (concat my-home-dir "/lib/elisp"))
;; (setq load-path (cons my-elisp-dir load-path))
;; (load "draco-rtt")
;;---------------------------------------------------------------------------;;

;;
;; get the Elisp files
;;

(if want-rtt-font-lock-keywords
    (require 'fl-keywords))
(require 'rtt-hacks)
(require 'tme-hacks)
(require 'draco-hacks)

;;
;; do any setups
;;

(rtt-setup)
(tme-setup)
(draco-setup)

;; RTT Menu

(if want-rtt-menu
    (progn
      (require 'rtt-menu)
      (rtt-menu-setup)
      ))

;;
;; Final Message
;;

(provide 'draco-rtt)

(message "DRACO-RTT Environment has been configured.")

;;---------------------------------------------------------------------------;;
;; end of draco-rtt.el
;;---------------------------------------------------------------------------;;
