;;; infer-mode.el
;;; This file is Emacs LIsp!!!

(defvar infer-mode-path-alist
  '(
    ("^/usr/include/CC/.*\\.h$" . (lambda ()
				   (if (not (string= mode-name "C++"))
				       (progn
					 (c++-mode)
					 (font-lock-fontify-buffer)))))
    )
  "Alist of file path patterns vs corresponding major mode functions.")

(defvar infer-mode-content-alist 
  '(
    ("^;;; This file is Emacs LIsp!!!$" . emacs-lisp-mode )
    ("^\\.TH[ \t]+" . nroff-mode)
    ("^#![ \t]*.*/perl\\>" . perl-mode)
    ("^#![ \t]*.*python\\>" . python-mode)
    ("^#![ \t]*.*wish -f" . tcl-mode)
    ("^#![ \t]*.*/sh" . sh-mode)
    ("The POOMA Framework" . (lambda ()
			       (if (not (string= mode-name "C++"))
				   (progn
				     (c++-mode)
				     (font-lock-fontify-buffer)))))
    ("KAI C\\+\\+ Compiler" . (lambda ()
			       (if (not (string= mode-name "C++"))
				   (progn
				     (c++-mode)
				     (font-lock-fontify-buffer)))))
    )
  "Alist of file content patterns vs corresponding major mode functions.
Each element looks like (REGEXP . FUNCTION).
Visiting a file whose which contains REGEXP within the first infer-mode-limit
characters causes FUNCTION to be called.")

(defvar infer-mode-content-limit 1000
  "Number of characters to search when inferring mode based on file contents.
Nil means search entire buffer.  See infer-mode-alist.")


(defun infer-file-mode ()
  "Infer mode of file based on its contents.
To use, add this function to find-file-hooks.
To customize, see the variable infer-mode-alist."

  (catch 'found-by-path
    ;; Infer file mode based on its path.
    (let ((bufnam (buffer-file-name))
	  (path-alist infer-mode-path-alist))
      (catch 'found
	(while path-alist
	(if (string-match (car (car path-alist)) bufnam)
	    (throw 'found t)
	  (setq path-alist (cdr path-alist)))))
      (if path-alist
	  (progn
	  (funcall (cdr (car path-alist)))
	  (throw 'found-by-path t))))
    
    ;; Okay, mode not infered on basis of path, now inspect buffer
    ;; contents. 
    (let ((mode-alist infer-mode-content-alist))
      (catch 'found
	(while mode-alist
	  (if (save-excursion
		(goto-char (point-min))
		(re-search-forward (car (car mode-alist)) 
				   infer-mode-content-limit t))
	      (throw 'found t)
	    (setq mode-alist (cdr mode-alist)))))
      (if mode-alist
	  (funcall (cdr (car mode-alist)))))))

(add-hook 'find-file-hooks 'infer-file-mode)

(provide 'infer-mode)
