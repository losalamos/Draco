;;; pooma-hacks.el
;;; Geoffrey Furnish
;;; 17 June 1997
;;;
;;; The purpose of this file is to provide some elisp support for
;;; editing Pooma files in a way which conforms to the evidential norm
;;; for C++ coding by Pooma team members.  Specifically, a Pooma
;;; "style" for cc-mode is introduced, which is intended to document
;;; the specific settings that people use when editing Pooma sources.
;;; The reason this deserves to be promoted to the level of a cc-mode
;;; style, is because people who work on more projects than just
;;; Pooma, may well be using significantly different styles, and
;;; cannot be reasonably expected to "convert" to the Pooma editing
;;; style in to-to.  However, by introducing a formal "Pooma editing
;;; style", it becomes possible for Pooma sources to use the Emacs
;;; "Local variables" support, to enforce the Pooma style for Pooma
;;; source files, even when hacked on by poeple using other styles by
;;; default.  In other words, this will enable Pooma files to retain
;;; their normal appearance even when hacked upon by people with
;;; different personal preferences, without imposing burden upon said
;;; individuals. 

;;; USAGE:
;;; The intended usage paradigm is that everyone puts:
;;;         (require 'pooma-hacks)
;;;         (pooma-setup)
;;; somewhere in their .emacs.  This will make the Pooma editing style
;;; known to emacs, so that cc-mode will know what to do when it
;;; visits a file which contains a specification of the Pooma style in
;;; the Local variables section at the end.  Additionally, a variable
;;; "want-pooma-style-by-default" is available for customization.  If
;;; this variable is set to true, then code is added to the cc-mode
;;; hooks mechanism which causes Pooma style to be on by default for
;;; C++ files.  People who do not want the Pooma style by default
;;; simply do not set "want-pooma-style-by-default", but still retain
;;; the ability to nondestructively edit Pooma sources if they are
;;; adorned with the requisite support in the Local variables section.

(require 'cc-mode)

;;; A standard style for Pooma.

;;; NOTE: I (gmf) do not advocate this style.  This is merely an
;;; attempt to codify a specific quantification of exactly what it is
;;; that people are doing.  I do not actually know specifically what
;;; people are doing--these settings are inferences on my part based
;;; on what I see in the Pooma sources.  I hope other Pooma team
;;; members who actually know what they are doing, will update this
;;; record to approximate as closely as possible, the reality of
;;; coding in the Pooma framework.  

(defconst pooma-c-style
  `((c-basic-offset	. 2)
    (c-offsets-alist	. ((inline-open . 0)
			   (substatement-open . 0)))))

;;; The following is not being promoted to the actual level of a style
;;; in this file, since "Geoff's one true editing style" is actually
;;; being provided through rtt-hacks.el.  I include this here just so
;;; Pooma people can have a point of reference for what I would
;;; actually prefer to see adopted in the Pooma sources.  The most
;;; important features are that c-basic-offset is 4, not 2, and that
;;; comments are outdented one level.  The 4 space indenting is what I
;;; consider to be the net-wide emerging multi-language standard.  And
;;; the one-stop-outdenting of comments makes them stand out a little
;;; from the code without being stuck left.

(defconst geoffs-advocacy-c-style
  `((c-basic-offset	. 4)
    (c-offsets-alist	. ((access-label . -2)
			   (comment-intro . -)
			   (inline-open . 0)
			   (substatement-open . 0)))))

;;; Install the Pooma style so that it is known to Emacs cc-mode.

(c-add-style "POOMA" pooma-c-style nil)

;;; Provide a way for individuals to indicate if they want to use the
;;; Pooma style by default.  Presumably people who spend most of their
;;; time hacking the Pooma sources would put:
;;;        (setq want-pooma-style-by-default t)
;;; in their .emacs, and people who hack lots of stuff other than
;;; Pooma, would not.

(defvar want-pooma-style-by-default nil
  "Does the user want to use the Pooma editing style by default?")

;;; Define a standard hook for Pooma customization.  It would be
;;; highly desirable if this could be limited to installing the Pooma
;;; style, so that Pooma customizations can be fully recovered by
;;; people who do not choose to use this by default, simply by setting
;;; the style in the Local variables section.  If other customizations
;;; are made here, those would not be recoverable (as far as I know).
;;; So, anything beyond applying the style, should have good
;;; justification. 

(defun pooma-c-mode-common-hook ()
  "Perform Pooma customization of the cc-mode engine."
  (c-set-style "POOMA")			; This should do it, no?
  )

;;; Entry point for configuring Pooma customizations.  Call this from
;;; your .emacs somewhere.  

(defun pooma-setup ()
  (message "Configure POOMA editing style support")

  (if want-pooma-style-by-default
      (progn
	(add-hook 'c-mode-common-hook 'pooma-c-mode-common-hook)
	(message "Configuring to use Pooma C++ style by default.")))
  )

(provide 'pooma-hacks)
