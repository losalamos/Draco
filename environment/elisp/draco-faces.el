;; Set face colors
;; ===============

;; See definition of draco-faces in draco-setup.el

(defun draco-set-faces ()
  "Examine the value of variable \"draco-faces\" and execute the appropriate
function to set face colors and fonts."
  (interactive)
  (if (string= draco-faces "tme") (draco-faces-tme))
  (if (string= draco-faces "kgt") (draco-faces-kgt)))

(defcustom draco-faces "none"
  "*Set the default face/font values for draco-mode?  The 
function draco-set-faces is run by draco-mode-hook.
Known values are: 

tme \t -- runs draco-faces-tme()
kgt \t -- runs draco-faces-kgt()"
  :group 'draco-mode
  :type 'string)

(add-hook 'draco-mode-hook 'draco-set-faces)

;;---------------------------------------------------------------------------;;

(defun draco-faces-kgt ()
  "Set face defaults to KT's scheme."

  (interactive)
  (make-face 'cvs-unknown-face)

  (custom-set-faces
   '(bold-italic ((t (:foreground "yellow2" :italic nil ))) t )
   '(bold        ((t (:foreground "red"  :italic      nil      ))) t)
   '(custom-comment-face ((t (:foreground "purple" ))) t )
   '(cvs-msg-face                 ((t (:foreground "purple"             ))) t )
   '(default                      ((t (:foreground "black" :background "ivory"  ))) t )
   '(dired-date-time-face ((t (:foreground "midnightblue" :inverse-video nil))) t)
   '(dired-face-boring ((t (:foreground "lightblue" ))) t)
   '(dired-face-directory         ((t (:foreground "blue" ))) t )
   '(dired-face-executable        ((t (:foreground "red"             ))) t )
   '(dired-face-permissions       ((t (:foreground "midnightblue" :background "lightblue" ))) t )
   '(dired-face-setuid            ((t (:foreground "light sea green" ))) t )
   '(dired-face-symlink           ((t (:foreground "steel blue"      ))) t )
   '(dired-group-name-face ((t (:foreground "purple4" :bold t))) t)
   '(dired-owner-name-face ((t (:foreground "purple4" :bold t))) t)
   '(ediff-current-diff-face-A    ((t (:foreground "firebrick" :background "yellow"          ))) t )
   '(ediff-current-diff-face-B    ((t (:foreground "firebrick"  :background "yellow"          ))) t )
   '(ediff-even-diff-face-A       ((t (:foreground "forest green"   :background "light grey"      ))) t )
   '(ediff-even-diff-face-B       ((t (:foreground "forest green"   :background "light grey"      ))) t )
   '(ediff-fine-diff-face-B       ((t (:foreground "Navy"      :background "sky blue"        ))) t )
   '(ediff-fine-diff-face-C       ((t (:foreground "black"    :background "pale green"      ))) t )
   '(ediff-odd-diff-face-A        ((t (:foreground "purple"  :background "Turquoise"       ))) t )
   '(ediff-odd-diff-face-B        ((t (:foreground "purple"  :background "Turquoise"       ))) t )
   '(font-latex-sedate-face       ((t (:foreground "firebrick3"      ))) t )
   '(font-lock-builtin-face       ((t (:foreground "white"
						   :background "dark violet"  ))) t )
   '(font-lock-comment-face       ((t (:foreground "purple"            ))) t )
   '(font-lock-doc-string-face    ((t (:foreground "orange3" ))) t )
   '(font-lock-draco-dbc-face     ((t (:foreground "magenta" ))) t )
   '(font-lock-function-name-face ((t (:foreground "blue3"        ))) t )
   '(font-lock-keyword-face       ((t (:foreground "firebrick2"         ))) t )
   '(font-lock-preprocessor-face  ((t (:foreground "hotpink"       ))) t )
   '(font-lock-reference-face     ((t (:foreground "orange3"    ))) t )
   '(font-lock-string-face        ((t (:foreground "brown"      ))) t )
   '(font-lock-type-face          ((t (:foreground "forestgreen"      ))) t )
   '(font-lock-variable-name-face ((t (:foreground "RoyalBlue"       ))) t )
   '(highlight                    ((t (:foreground "red"             ))) t )
   '(hyper-apropos-documentation  ((t (:foreground "bisque3"         ))) t )
   '(hyper-apropos-hyperlink      ((t (:foreground "royalblue"       ))) t )
   '(info-node                    ((t (:foreground "mediumseagreen"  ))) t )
   '(info-xref                    ((t (:foreground "pink1"  :bold ))) t )
   '(isearch                      ((t (:foreground "white"    :background "dark violet"     ))) t )
   '(italic      ((t (:foreground "yellow2" :italic      nil :bold        nil      ))) t)
   '(man-bold                     ((t (:foreground "bisque3"         ))) t )
   '(man-heading                  ((t (:foreground "mediumseagreen"  ))) t )
   '(man-italic                   ((t (:foreground "bisque3"         ))) t )
   '(man-xref                     ((t (:foreground "lightskyblue"    ))) t )
   '(modeline                     ((t (:foreground "black"   :background "wheat"          ))) t )
   '(paren-blink-off              ((t (:foreground "red"  :background "ivory"           ))) t )
   '(paren-match                  ((t (:foreground "red"  :background "ivory"           ))) t )
   '(paren-mismatch               ((t (:foreground "black" :background "deeppink"        ))) t )
   '(primary-selection            ((t (:foreground "black" :background "MediumSeaGreen"))) t)
   '(shell-option-face            ((t (:foreground "royal blue"      ))) t )
   '(shell-output-2-face          ((t (:foreground "red"             ))) t )
   '(shell-output-3-face          ((t (:foreground "light sea green" ))) t )
   '(shell-output-face            ((t (:foreground "forestgreen"    :italic      nil              ))) t )
   '(shell-prompt-face            ((t (:foreground "blue" 			   ))) t )
   '(text-cursor ((t (:foreground "ivory1" :background "Blue4"))) t)
   '(zmacs-region                 ((t (:background "pale turquoise"  ))) t )
   )
 
;; Set faces fonts
;; ===============
  (set-face-font 'default     "6x13")
  (set-face-font 'bold        "6x13bold")
  (set-face-font 'dired-face-directory        "6x13bold")
  (set-face-font 'font-lock-draco-dbc-face    "6x13bold")
;(set-face-font 'default "-dt-interface user-medium-r-normal-xs sans-14-80-100-100-m-*-iso8859-1")
  (set-face-font 'paren-match "6x13bold")
  (set-face-font 'paren-mismatch "6x13bold")
  (set-face-font 'modeline    "-dt-application-bold-r-normal-sans-12-120-75-75-p-70-iso8859-1")
  (set-face-font 'modeline-buffer-id "-dt-application-bold-r-normal-sans-12-120-75-75-p-70-iso8859-1")
  (set-face-font 'modeline-mousable  "-dt-application-bold-r-normal-sans-12-120-75-75-p-70-iso8859-1")
  (set-face-font 'modeline-mousable-minor-mode "-dt-application-bold-r-normal-sans-12-120-75-75-p-70-iso8859-1")
  (set-face-font 'man-heading "helvetica-bold")
  (set-face-font 'man-bold    "-*-courier-bold-r-*-*-*-120-*-*-*-*-iso8859-*")
  (set-face-font 'bold        "6x13bold")
  (set-face-font 'bold-italic "6x13bold")
  (set-face-font 'shell-output-face "6x13")
  (set-face-font 'ediff-odd-diff-face-A "6x13")
  (set-face-font 'ediff-odd-diff-face-B "6x13")
  (set-face-font 'ediff-even-diff-face-A "6x13")
  (set-face-font 'ediff-even-diff-face-B "6x13")
  (set-face-font 'info-node "6x13")
  (set-face-font 'info-xref "6x13")
  (set-face-font 'cvs-msg-face "6x13")
)

;;---------------------------------------------------------------------------;;

(defun draco-faces-tme ()
  "Set face defaults to Tom's scheme."

  (custom-set-faces
   '(default ((t (:foreground "black" :background "ivory1" :size "14pt"))) t)
   '(cvs-msg-face ((t (:foreground "purple" :background "" :dim t))))
   '(dired-date-time-face ((t (:foreground "midnightblue" :inverse-video nil))) t)
   '(dired-face-boring ((((class color)) (:foreground "lightblue"))))
   '(dired-face-directory ((t (:foreground "blue" :bold t))))
   '(dired-face-flagged ((((class color)) (:foreground "yellow" :background "red" :size "12pt" :italic t))))
   '(dired-face-marked ((((class color)) (:foreground "black"   :background "gold"))))
   '(dired-face-permissions ((t (:foreground "midnightblue" :background "lightblue"))))
   '(dired-group-name-face ((t (:foreground "purple4" :bold t))) t)
   '(dired-owner-name-face ((t (:foreground "purple4" :bold t))) t)
   '(font-latex-sedate-face ((((class color) (background light)) (:foreground "firebrick3"))))
   '(font-lock-comment-face ((((class color) (background light)) (:foreground "purple"))))
   '(font-lock-doc-string-face ((((class color) (background light)) (:foreground "orange3"))))
   '(font-lock-function-name-face ((((class color) (background light)) (:foreground "blue3"))))
   '(font-lock-keyword-face ((((class color) (background light)) (:foreground "firebrick2"))))
   '(font-lock-kull-macros-face ((t (:foreground "magenta"))) t)
   '(font-lock-preprocessor-face ((((class color) (background light)) (:foreground "hotpink"))))
   '(font-lock-reference-face ((((class color) (background light)) (:foreground "orange3"))))
   '(font-lock-string-face ((((class color) (background light)) (:foreground "brown"))))
   '(font-lock-type-face ((((class color) (background light)) (:foreground "forestgreen"))))
   '(font-lock-variable-name-face ((((class color) (background light)) (:foreground "RoyalBlue"))))
   '(info-xref ((t (:foreground "pink1" :bold t))))
   '(modeline ((t (:foreground "black" :background "wheat" :family "times" :bold t))) t)
   '(primary-selection ((t (:foreground "black" :background "MediumSeaGreen"))) t)
   '(text-cursor ((t (:foreground "ivory1" :background "Blue4"))) t)
   '(zmacs-region ((t (:foreground "black" :background "LightBlue"))) t)
   ))


;;---------------------------------------------------------------------------;;

(provide 'draco-faces)
