;;---------------------------------------------------------------------------;;
;; draco-menu.el
;; Kelly (K.T.) Thompson
;; May 25, 2001
;;---------------------------------------------------------------------------;;
;; DRACO MENU
;;---------------------------------------------------------------------------;;

(require 'easymenu)

(defun draco-menu-new-files () "Submenu to use for new files."
  (list "New file..."
   ["New C++ package"              draco-package      t]
   ["New C++ package/test"         draco-package-test t]
   ["New C++ package/autodoc"      draco-package-doc  t]
   ["New C++ class"     draco-class     t]
   ["New C++ header"    draco-cc-head   t]
   ["New C++ header.in" draco-cc-headin t]
   ["New C header"      draco-c-head    t]
   ["New C header.in"   draco-c-headin  t]
   ["New C++ template impl. file (.t.hh)" draco-cc-imp  t]
   ["New C++ instantiation file (_pt.cc)" draco-cc-pt   t]
   ["New C++ unit test executable"        draco-cc-test t]
   ["New C++ class parser"        draco-parser t]
   "----"
   ["Insert Name and Time" draco-name-and-time t]
   ["New Python file"                     draco-python  t]
   ["New specialized makefile"            draco-make    t]))

(defun draco-menu-new-latex-files () "Submenu to use for new LaTeX files."
  (list "New LaTeX document..."
	["New memo"                 draco-memo     t]
	["New research note"        draco-note     t]
	["New article"              draco-article  t]
        ["New report"               draco-report   t]
        ["New bibliography"         draco-bib      t]
        ["New vision and scope"     draco-viscope  t]
        ["New bug post-mortem memo" draco-bug-pm   t]))

(defun draco-menu-insert-comments-default ()
  "Submenu for inserting comments (context sensitive)."
  (list "Insert comment..."
	["Insert C++ comment divider" draco-insert-comment-divider t]))

(defun draco-menu-insert-extras() "Submenu with extra stuff."
  (list "Extras..."
	["Find companion file"  draco-find-companion-file t]
	["Customize Draco-mode" (customize-group 'draco-mode) t]
	["Clang-format region"  clang-format-region t]
	["Clang-format buffer"  clang-format-buffer t]))

(defvar draco-menu nil
  "The Draco menu for XEmacs.")

(defun draco-mode-update-menu ( submenu1 )
  "Update the Draco-mode menu based on current context."
  (interactive)
  (easy-menu-remove draco-menu) ;; kill old menu
  (easy-menu-define draco-menu global-map "Draco"
	(list "Draco"
	(append (draco-menu-new-files))
	(append (draco-menu-new-latex-files))
	(append submenu1)
	(append (draco-menu-insert-extras))))
  (easy-menu-add draco-menu) ;; put the menu back
)

(if want-draco-menu
    (add-hook 'draco-mode-hook
	      '(lambda ()
		 (draco-mode-update-menu (draco-menu-insert-comments-default)))))

;;---------------------------------------------------------------------------;;
;; Draco toolbar
;;---------------------------------------------------------------------------;;

(defcustom draco-default-toolbar-spec
  '([toolbar-file-icon     find-file                  t "Open a file"           ]
    [toolbar-folder-icon   toolbar-dired              t "View a directory"      ]
    [toolbar-disk-icon     toolbar-save               t "Save buffer"           ]
    [toolbar-printer-icon  toolbar-print              t "Print region/buffer"   ]
    [toolbar-cut-icon      toolbar-cut                t "Kill region"           ]
    [toolbar-copy-icon     toolbar-copy               t "Copy region"           ]
    [toolbar-paste-icon    toolbar-paste              t "Paste from clipboard"  ]
    [toolbar-undo-icon     toolbar-undo               t "Undo edit"             ]
    [toolbar-spell-icon    toolbar-ispell             t "Spellcheck"            ]
    [toolbar-replace-icon  toolbar-replace            t "Replace text"          ]
    [toolbar-mail-icon     toolbar-mail               t "Mail"                  ]
    [toolbar-news-icon     toolbar-news               t "Gnu News"              ]
    [toolbar-info-icon     toolbar-info               t "Information"           ]
    [toolbar-compile-icon  toolbar-compile            t "Compile"               ]
    [toolbar-debug-icon    toolbar-debug              t "Debug"                 ]
    )
  "The initial settings for draco-toolbar-spec."
  :group 'draco-mode
  :type 'list)

(defcustom draco-toolbar-spec
  '([toolbar-file-icon     find-file                  t "Open a file"           ]
    [toolbar-folder-icon   toolbar-dired              t "View a directory"      ]
    [toolbar-disk-icon     toolbar-save               t "Save buffer"           ]
    [toolbar-printer-icon  toolbar-print              t "Print region/buffer"   ]
    [toolbar-cut-icon      toolbar-cut                t "Kill region"           ]
    [toolbar-copy-icon     toolbar-copy               t "Copy region"           ]
    [toolbar-paste-icon    toolbar-paste              t "Paste from clipboard"  ]
    [toolbar-undo-icon     toolbar-undo               t "Undo edit"             ]
    [toolbar-spell-icon    toolbar-ispell             t "Spellcheck"            ]
    [toolbar-replace-icon  toolbar-replace            t "Replace text"          ]
    [toolbar-info-icon     toolbar-info               t "Information"           ]
    [toolbar-compile-icon  toolbar-compile            t "Compile"               ]
    )
  "A specialized version of draco-default-toolbar-spec."
  :group 'draco-mode
  :type 'list)

;; use my toolbar as default toolbar
(interactive (set-specifier default-toolbar draco-toolbar-spec))

;;---------------------------------------------------------------------------;;
;; XEmacs frame title
;;---------------------------------------------------------------------------;;

(defun short-system-name ()
  "Return the first word of a machine's name."
  (substring (system-name) 0 (string-match "\\." (system-name))))

(if (featurep 'xemacs)
    (setq frame-title-format (concat "[" (short-system-name) "] %f")))

; default value is "%S: %b"
;
;(setq frame-title-format
;      '("%S: " (buffer-file-name "%f" (dired-directory
;				       dired-directory "%b"))))
; Formats:
;
; %b - buffer name
; %c - column number
; %l - line number
; %f - file name
; %s - status
; %p - precent of buffer
; %S - selected X-Windows Frame

(provide 'draco-menu)

;;---------------------------------------------------------------------------;;
;; end of draco-menu.el
;;---------------------------------------------------------------------------;;


;; Notes

; See
; /usr/share/xemacs/xemacs-packages/lisp/view-process/view-process-xemacs.el
; for more examples on using radio buttons, etc in menus.
