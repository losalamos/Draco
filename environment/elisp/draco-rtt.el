;;; draco-rtt.el
;;; Tom Evans
;;; Feb 24, 1999

;;; $Id$

;;----------------------------------
;; Load the DRACO-RTT Elisp packages
;; ---------------------------------

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
