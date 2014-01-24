%% LaTeX2e letter class based on the standard letter.cls and work done
%% by ryxm of T-2 13 Jan 95
%% ---------------------------------------------------------------------
%% The use of \begin{letter}...\end{letter} has been eliminated in 
%% favor of giving the address in a command.  The salutation in
%% \opening has been removed, and the complimentary close in \closing
%% has been removed.  These three changes are consistant with the
%% old LANL letter style.  The \opening command writes a LANL logo
%% and adds an equal opportunity message to the first page footer.
%% Additional new commands have been added to imitate laletter.sty.
%% Use the template /usr/local/teTeX/texmf/tex/latex/lettest.tex
%% on the IBM Cluster to have easy access to the options.
%% ---------------------------------------------------------------------
\typeout{^^J   Document Class 'laletter'. October 1996^^J
 *****  This is a DIFFERENT version of the LANL letter style  *****^^J
    See /usr/local/teTeX/texmf/tex/latex/local/lettest.tex^^J
    for a template to start from.^^J}
%%
\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{laletter}[1996/10/01 LANL Letter Class]
%% Load the normal letter class.
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{letter}}
\ProcessOptions
\LoadClass{letter}
%% Modify the normal text area on the page.
\setlength{\textwidth}{6.5in}
\setlength{\textheight}{8.5in}
\setlength{\oddsidemargin}{0in}
\setlength{\topmargin}{0pt}
\setlength{\headsep}{12pt}
\pagestyle{headings}
%% Define new commands for the LANL letter class.
\newcommand*{\complimentaryclose}[1]{\def\cclose{#1}}
\newcommand*{\originator}[1]{\def\originatr{\uppercase{#1}}}
\newcommand*{\typist}[1]{\def\typst{#1}}
\newcommand*{\towhom}[1]{\def\toaddr{#1}\@processto{#1}}
\newcommand*{\referto}[1]{\def\refno{#1}}
\newcommand*{\salutation}[1]{\def\salutatn{#1}}
\newcommand*{\subject}[1]{\def\subjectline{#1}}
\newcommand{\distribution}[1]{\def\distrib{#1}}
\newcommand{\cy}[1]{\def\cyy{#1}}
\newcommand{\enc}[1]{\def\encl{#1}}
\newcommand{\encas}{\def\encass{Enc. a/s}}
\newcommand{\attachments}[1]{\def\attach{#1}}
\newcommand{\attachmentas}{\def\attachas{Attachment a/s}}
\newcommand{\attachmentsas}{\def\attachsas{Attachments a/s}}
%% Modify some commands for the T-2 letter class.  The \opening
%% command is modified to add a LANL logo and standard address lines.
%% The header and footer lines are modified to LANL standards.
%% The letter environment still exists, but we do not use it.
\renewcommand*{\opening}{
  \thispagestyle{firstpage}%
  \parbox[t]{2.0in}{
     {\sffamily\bfseries\fontsize{25}{30}\selectfont Los Alamos}\\
     {\sffamily\fontsize{10}{11}\selectfont N\kern.07em A\kern.07em T\kern.07em I\kern.07em O\kern.07em N\kern.07emA\kern.07em L\kern.07em ~\kern.07em L\kern.07em A\kern.07em B\kern.07em O\kern.07em R\kern.07em A\kern.07em T\kern.07em O\kern.07em R\kern.07em Y}}
  \par
  \begin{minipage}{4.35in}
     \scriptsize\fromaddress
  \end{minipage}
  \begin{minipage}[t]{2.1in}
     \begin{tabular}{rl}
        {\scriptsize Date:}&\@date\\
        {\scriptsize Refer to:}&\refno
     \end{tabular}
  \end{minipage}
  \par\vspace{2\parskip}%
  {\raggedright \toaddr \par}%
  \vspace{2\parskip}%
  \salutatn\par\nobreak\par%
  {\bf SUBJECT: }\subjectline}
\renewcommand{\closing}{
  \stopbreaks
  \noindent
  \hspace*{\longindentation}
  \parbox{\indentedwidth}{\raggedright
       \vspace{11pt}\cclose\\[6\medskipamount]%
       \ifx\@empty\fromsig
           \fromname
       \else \fromsig \fi\strut}%
   \par
   \noindent
   \originatr:\typst
   \ifx\attach\@empty\else
     \par 
     \noindent
     Attachments:\\
     \attach
   \fi
   \ifx\attachas\@empty\else
     \par 
     \noindent
     \attachas
   \fi
   \ifx\attachsas\@empty\else
     \par 
     \noindent
     \attachsas
   \fi
   \ifx\encl\@empty\else
     \par 
     \noindent
     Enclosures:\\
     \encl
   \fi
   \ifx\encass\@empty\else
     \par 
     \noindent
     \encass
   \fi
   \ifx\distrib\@empty\else
     \par 
     \noindent
     Distribution:\\
     \distrib
   \fi
   \ifx\cyy\@empty\else
     \par 
     \noindent
     Cy:\\
     \cyy
   \fi}
\renewcommand{\ps@headings}{%
      \let\@oddfoot\@empty\let\@evenfoot\@empty
      \def\@oddhead{\slshape\headtoname{} \ignorespaces\toname
                    \hfil \pagename{} \thepage
                    \hfil \@date}%
      \let\@evenhead\@oddhead}
\renewcommand{\ps@firstpage}{%
     \let\@oddhead\@empty
     \def\@oddfoot{
        \footnotesize
        \hfill 
        An Equal Opportunity Employer/Operated 
        by the Los Alamos National Security, LLC
        \hfill}}
% Initializations
\distribution{}
\cy{}
\enc{}
\attachments{}
\def\attachas{}
\def\attachsas{}
\def\encass{}
\def\toname{}
\def\tojunk{}
\def\through{}
\def\originatr{}
\def\typst{}
\pagestyle{headings}
