%  ========================================================================  %
% 
%       Author: Mark G. Gray
%               Los Alamos National Laboratory
%       Date:   Wed Jan 13 09:10:59 MST 1999
% 
%       Copyright (c) 1999 U. S. Department of Energy. All rights reserved.
% 
%       $Id$
%
% memo class based on standard t2memo.cls   mgg 13 Jan 99
% ---------------------------------------------------------------
% This is based on the standard latex article class with changes
% to make it act more like the t2memo class which is in turn based
% on the Sydoriak LANL memo style.
% Use the file larmemo.tex as a memo template.  It also contains
% some comments on special features of the class.
%
% ========================================================================  %

			  %% identification
\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{larmemo}[1999/01/13 My Memo Class]

			      %% options
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{article}}
\ProcessOptions

			       %% main
\LoadClass{article}
\RequirePackage{color}

%% Modify the normal text area on the page.
\setlength{\textwidth}{6.5in}
\setlength{\textheight}{8.5in}
\setlength{\oddsidemargin}{0in}
\setlength{\topmargin}{0pt}
\setlength{\headheight}{2.5\baselineskip}
\setlength{\headsep}{12pt}
\setlength{\parskip}{0.7em}	% Block paragraph with
\setlength{\parindent}{0pt}	% no indent

%% Add new commands for the memo class.
\newcommand*{\divisionname}[1]{\def\divname{#1}}
\newcommand*{\groupname}[1]{\def\grpname{#1}}
\newcommand*{\fromms}[1]{\def\fromwho{#1}}
\newcommand*{\phone}[1]{\def\phonenum{#1}}
\newcommand*{\originator}[1]{\def\originatr{\uppercase{#1}}}
\newcommand*{\typist}[1]{\def\typst{#1}}
\newcommand*{\toms}[1]{\def\towho{#1}\gettoname{#1}}
\newcommand*{\thru}[1]{\def\through{#1}}
\newcommand*{\refno}[1]{\def\refnum{#1}}
\newcommand*{\subject}[1]{\def\subjectline{#1}}
\newcommand{\distribution}[1]{\def\distrib{#1}}
\newcommand{\cy}[1]{\def\cyy{#1}}
\newcommand{\enc}[1]{\def\encl{#1}}
\newcommand{\encas}{\def\encass{Enc. a/s}}
\newcommand{\attachments}[1]{\def\attach{#1}}
\newcommand{\attachmentas}{\def\attachas{Attachment a/s}}
\newcommand{\attachmentsas}{\def\attachsas{Attachments a/s}}
\newcommand{\HRule}{\rule{\linewidth}{1pt}}
\newcommand{\Secret}{\def\classification{SECRET}}

% various pagestyle headers
% Unclassified

\newcommand{\ps@larhead}{
    \renewcommand{\@oddhead}{\refnum \hfill -\thepage -\hfill \@date}
    \renewcommand{\@oddfoot}{}
    \renewcommand{\@evenhead}{\@oddhead}
    \renewcommand{\@evenfoot}{\@oddfoot}
}

% Secret Restricted Data

\newcommand{\SRD}{\hfill \large\scshape
                  \textcolor{red}{SECRET/RD} \hfill}

\newcommand{\ps@secret}{
    \renewcommand{\@oddhead}{
        \begin{minipage}{\textwidth}
	    \begin{tabular*}{\textwidth}{l@{\extracolsep\fill}cr}
                &\SRD& \\
                \refnum & -\thepage - & \@date
            \end{tabular*}
         \end{minipage}
    }
    \renewcommand{\@oddfoot}{\SRD}
    \renewcommand{\@evenhead}{\@oddhead}
    \renewcommand{\@evenfoot}{\@oddfoot}
}

\newcommand{\ps@secretfirst}{
    \renewcommand{\@oddhead}{\SRD}
    \renewcommand{\@oddfoot}{\SRD}
    \renewcommand{\@evenhead}{\@oddhead}
    \renewcommand{\@evenfoot}{\@oddfoot}
}


% redefine section and subsections
\renewcommand{\thesection}{}
\renewcommand{\section}{\@startsection {section}{1}%
  {-2ex}%
  {\baselineskip}%
  {0.1\baselineskip}%
  {\normalfont\bfseries}}

\renewcommand{\thesubsection}{}
\renewcommand{\subsection}{\@startsection {subsection}{2}%
  {6pt}%
  {\baselineskip}%
  {-6pt}%
  {\normalfont\bfseries}}

% maketitle command to add logo and info block
\renewcommand{\maketitle}{
  \newcommand{\nullmacro}{}
  \ifx\@oddfoot\nullmacro
      \thispagestyle{empty}
  \else
      \thispagestyle{secretfirst}
  \fi
  \noindent
  \begin{minipage}{3.25in}
     {\sffamily\bfseries\fontsize{25}{30}\selectfont Los\kern-.1em ~Alamos}\\
     {\sffamily\fontsize{10}{11}\selectfont N\kern.06em A\kern.03em T\kern.06em I\kern.06em O\kern.06em N\kern.06em A\kern.06em L\kern.06em ~\kern.03em L\kern.06em A\kern.06em B\kern.06em O\kern.06em R\kern.06em A\kern.03em T\kern.03em O\kern.06em R\kern.06em Y}\\
     {\sffamily\bfseries\fontsize{21}{21}\selectfont 
      \rule{0pt}{21pt}research\\memo\kern.02em randum}\\
     {\rmfamily\itshape\fontsize{10}{11}\selectfont \divname}\\
     {\rmfamily\itshape\fontsize{10}{11}\selectfont \grpname}
    \end{minipage}
   \ifx\through\@empty
    \begin{minipage}[t]{3.25in}
      \begin{tabular}{rl}
        {\scriptsize To/MS:}&\begin{minipage}[t]{2.75in}
                              \towho
                             \end{minipage}\\
        {\scriptsize From/MS:}&\begin{minipage}[t]{2.75in}
                              \fromwho
                             \end{minipage}\\
        {\scriptsize Phone/FAX:}&\phonenum\\
        {\scriptsize Symbol:}&\refnum\\
        {\scriptsize Date:}&\@date
      \end{tabular}
    \end{minipage}
   \else
    \begin{minipage}[t]{3.25in}
      \begin{tabular}{rl}
        {\scriptsize To/MS:}&\begin{minipage}[t]{2.75in}
                              \towho
                             \end{minipage}\\
        {\scriptsize Thru/MS:}&\begin{minipage}[t]{2.75in}
                                \through
                               \end{minipage}\\
        {\scriptsize From/MS:}&\begin{minipage}[t]{2.75in}
                                \fromwho
                               \end{minipage}\\
        {\scriptsize Phone/FAX:}&\phonenum\\
        {\scriptsize Symbol:}&\refnum\\
        {\scriptsize Date:}&\@date
      \end{tabular}
    \end{minipage}
   \fi
%   \vspace{12pt}
   \par
%   \noindent
   \par
   \HRule
%   \vspace{12pt}
   \par
   \textbf{\large\subjectline}
   \par
  }

% makeclosing command to add additional memo options.
\newcommand{\makeclosing}{
   \newpage
   \par
   \nopagebreak
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
   \fi
   }

% Procedure to strip the toname out of "toname, group, MS".
\def\Dist{Distribution}
\def\gettoname#1{\ifx\towho\Dist\let\toname=\towho\else\strip #1\par\fi}
\long\def\strip #1,#2\par{\def\toname{#1}\def\tojunk{#2}}

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
\pagestyle{larhead}
\AtEndDocument{\makeclosing}
