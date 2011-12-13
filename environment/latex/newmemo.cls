%%---------------------------------------------------------------------------%%
%% NEWMEMO.CLS:
%% ----------
%% My additions to the t2memo class for memorandum.
%% ADDITIONS:
%%           1) added additional sections
%%           2) added float figres and tables
%%           3) added bibtex (bibliography environment)
%%           4) AW: Added the footer; reformatted the header and fonts
%%           5) AW: Removed macro shortcuts \@xp (\expandafter) and \@nx (\noexpand)
%%---------------------------------------------------------------------------%%
%%
%% t-2 memo class based on standard letter.cls   ryxm   13 Jan 95
%% ---------------------------------------------------------------
%% This is based on the standard latex letter class with changes
%% to make it act more like the Sydoriak LANL memo style.
%% Use the file t2memo.tex as a memo template.  It also contains
%% some comments on special features of the class.
%% ---------------------------------------------------------------
\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{newmemo}
        [2010/07/22 LANL Memorandum Class, T.M. Evans, A. Wollaber]
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{letter}}
\ProcessOptions
\LoadClass{letter}
\RequirePackage{graphicx}
%% Apparantly "palatino" is now the official font..
\renewcommand{\rmdefault}{ppl}
\renewcommand{\sfdefault}{phv}    % Helvetica (close enough to Arial)
%% Modify the normal text area on the page.
\setlength{\textwidth}{7in}
\setlength{\textheight}{8.5in}
\setlength{\oddsidemargin}{0in}
\setlength{\topmargin}{0pt}
\setlength{\headheight}{3.5ex}
\setlength{\headsep}{18pt}
%% Add new commands for the LANL memo class.
\newcommand*{\groupname}[1]{\def\grpname{#1}}
\newcommand*{\groupmail}[1]{\def\grpmail{#1}}
\newcommand*{\groupphone}[1]{\def\grpphone{#1}}
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
% Revise opening command to add logo and info block
\renewcommand*{\opening}{
    \thispagestyle{firstpage}
    \begin{minipage}{2.95in}
     \vspace{-8em}
     \hspace{-4em}
     {\includegraphics[width=2.6in]{logo}}\\
     {\sffamily\fontsize{20}{20}\selectfont \rule{0pt}{21pt}memorandum}\\
     {\sffamily \itshape\fontsize{10}{11}\selectfont \grpname}\\
     {\itshape\fontsize{11}{12}\selectfont \grpmail}\\
     {\itshape\fontsize{11}{12}\selectfont \grpphone}
    \end{minipage}
   \ifx\through\@empty
    \begin{minipage}[t]{3.5in}
     \vspace{-8em}
      \begin{tabular}{rl}
        {\sffamily \scriptsize \itshape To/MS:}&\begin{minipage}[t]{2.75in}
               \sffamily \towho
                             \end{minipage}\\
        {\sffamily \scriptsize \itshape From/MS:}&\begin{minipage}[t]{2.75in}
             \sffamily \fromwho
                             \end{minipage}\\
        {\sffamily \scriptsize \itshape Phone/FAX:}&{\sffamily \phonenum}\\
        {\sffamily \scriptsize \itshape Symbol:}&{\sffamily \refnum }\\
        {\sffamily \scriptsize \itshape Date:}& {\sffamily \@date }
      \end{tabular}
    \end{minipage}
   \else
    \begin{minipage}[t]{3.5in}
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
    \par
    \vspace{24pt}
    \textbf{\large \subjectline}
    \par
  }
% Revise closing command to add additional memo options.
\renewcommand*{\closing}{
   \stopbreaks
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
   \fi
   }
% Procedure to strip the toname out of "toname, group, MS".
\def\Dist{Distribution}
\def\gettoname#1{\ifx\towho\Dist\let\toname=\towho\else\strip #1\par\fi}
%\long\def\strip #1,#2@@@{\def\toname{#1}\def\tojunk{#2}}
\long\def\strip #1,#2\par{\def\toname{#1}\def\tojunk{#2}}
% Revise Headings
% \refnum
\if@twoside
  \def\ps@headings{%
      \let\@oddfoot\@empty\let\@evenfoot\@empty
      \def\@oddhead{\slshape%
        \begin{tabular}{l}
          \headtoname{} \ignorespaces\toname \\
          {footnotesize\refnum}
        \end{tabular}
                    \hfil --\thepage--
                    \hfil \@date}%
      \let\@evenhead\@oddhead}
\else
  \def\ps@headings{%
      \def\@oddfoot{\sffamily \fontsize{8}{8}\selectfont 
       \hfil An Equal Opportunity Employer / Operated by 
        Los Alamos National Security, LLC for 
        the U.S.~Department of Energy's NNSA }
      \def\@oddhead{\slshape%
        \begin{tabular}{l}
          \headtoname{} \ignorespaces\toname \\
          {\footnotesize\refnum}
        \end{tabular}
                    \hfil --\thepage--
                    \hfil \@date}}
\fi
\def\ps@firstpage{%
   \let\@oddhead\@empty 
   \def\@oddfoot{\sffamily \fontsize{8}{8}\selectfont 
    \hfil An Equal Opportunity Employer / Operated by 
     Los Alamos National Security, LLC for 
     the U.S. Department of Energy's NNSA }
    }
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
\pagestyle{headings}
%--Bibliographystyle
\newdimen\bibindent
\setlength\bibindent{1.5em}
\newcommand\refname{References}
\newenvironment{thebibliography}[1]
%     {\section*{\refname
%        \@mkboth{\MakeUppercase\refname}{\MakeUppercase\refname}}%
      {\normalfont\bfseries\refname\normalfont
      \list{\@biblabel{\@arabic\c@enumiv}}%
           {\settowidth\labelwidth{\@biblabel{#1}}%
            \leftmargin\labelwidth
            \advance\leftmargin\labelsep
            \@openbib@code
            \usecounter{enumiv}%
            \let\p@enumiv\@empty
            \renewcommand\theenumiv{\@arabic\c@enumiv}}%
      \sloppy\clubpenalty4000\widowpenalty4000%
      \sfcode`\.\@m}
     {\def\@noitemerr
       {\@latex@warning{Empty `thebibliography' environment}}%
      \endlist}
\newcommand\newblock{\hskip .11em\@plus.33em\@minus.07em}
\let\@openbib@code\@empty
%--Figures and Tables
\newcommand\figurename{FIG.}
\newcommand\tablename{TABLE}
\newcounter{figure}
\renewcommand\thefigure{\@arabic\c@figure}
\def\fps@figure{tbp}
\def\ftype@figure{1}
\def\ext@figure{lof}
\def\fnum@figure{\figurename~\thefigure}
\newenvironment{figure}
               {\@float{figure}}
               {\end@float}
\newenvironment{figure*}
               {\@dblfloat{figure}}
               {\end@dblfloat}
\newcounter{table}
\renewcommand\thetable{\@arabic\c@table}
\def\fps@table{tbp}
\def\ftype@table{2}
\def\ext@table{lot}
\def\fnum@table{\tablename~\thetable}
\newenvironment{table}
               {\@float{table}}
               {\end@float}
\newenvironment{table*}
               {\@dblfloat{table}}
               {\end@dblfloat}
\newlength\abovecaptionskip
\newlength\belowcaptionskip
\setlength\abovecaptionskip{10\p@}
\setlength\belowcaptionskip{0\p@}
\long\def\@makecaption#1#2{%
  \vskip\abovecaptionskip
  \sbox\@tempboxa{#1: #2}%
  \ifdim \wd\@tempboxa >\hsize
    #1: #2\par
  \else
    \global \@minipagefalse
    \hb@xt@\hsize{\hfil\box\@tempboxa\hfil}%
  \fi
  \vskip\belowcaptionskip}
%--Sections
%-DEFINITION OF SECTIONS
\def\@startsection#1#2#3#4#5#6{%
 \if@noskipsec \leavevmode \fi
 \par \@tempskipa #4\relax
 \@afterindenttrue
 \ifdim \@tempskipa <\z@ \@tempskipa -\@tempskipa \@afterindentfalse\fi
 \if@nobreak \everypar{}\else
     \addpenalty\@secpenalty\addvspace\@tempskipa\fi
 \@ifstar{\@dblarg{\@sect{#1}{\@m}{#3}{#4}{#5}{#6}}}%
         {\@dblarg{\@sect{#1}{#2}{#3}{#4}{#5}{#6}}}%
}
\def\@sect#1#2#3#4#5#6[#7]#8{%
  \edef\@toclevel{\ifnum#2=\@m 0\else\number#2\fi}%
  \ifnum #2>\c@secnumdepth \let\@secnumber\@empty
  \else \expandafter \let \expandafter \@secnumber\csname the#1\endcsname\fi
 \ifnum #2>\c@secnumdepth
   \let\@svsec\@empty
 \else
    \refstepcounter{#1}%
    \edef\@svsec{\ifnum#2<\@m
       \@ifundefined{#1name}{}{%
         \ignorespaces\csname #1name\endcsname\space}\fi
       \noexpand \textup{%
      \noexpand \mdseries
         \csname the#1\endcsname.}\enspace
    }%
  \fi
  \@tempskipa #5\relax
  \ifdim \@tempskipa>\z@ % then this is not a run-in section heading
    \begingroup #6\relax
    \@hangfrom{\hskip #3\relax\@svsec}{\interlinepenalty\@M #8\par}%
    \endgroup
    %\ifnum#2>\@m \else \@tocwrite{#1}{#8}\fi
  \else
  \def\@svsechd{#6\hskip #3\@svsec
    \@ifnotempty{#8}{\ignorespaces#8\unskip
       \@addpunct.}%
    %\ifnum#2>\@m \else \@tocwrite{#1}{#8}\fi
  }%
  \fi
  \global\@nobreaktrue
  \@xsect{#5}}
\let\@ssect\relax
%-DEFINITION OF SECTIONS
\setcounter{secnumdepth}{3} \newcounter{section}
\newcommand{\section}{\@startsection {section} {1} {0em}
  {\baselineskip} {\baselineskip} {\center\normalfont\large\bfseries}}
\newcounter{subsection}[section]
\newcommand{\subsection}{\@startsection {subsection} {2} {0em}
  {.5\baselineskip} {.5\baselineskip}
  {\normalfont\bfseries\normalsize}}
\newcounter{subsubsection}[subsection]
\newcommand{\subsubsection}{\@startsection {subsubsection} {3} {0em}
  {-.5\baselineskip} {-2ex} {\normalfont\itshape\normalsize}}
\renewcommand\thesection {\@arabic\c@section}
\renewcommand\thesubsection {\thesection.\@arabic\c@subsection}
\renewcommand\thesubsubsection{\thesubsection.\@arabic\c@subsubsection}
\newcommand{\appendix}{
  \setcounter{section}{0}
  \renewcommand{\thesection}{\Alph{section}}}
\newcommand\abstractname{Executive Summary}
\newenvironment{abstract}{
  \small
  \begin{center}%
    {\bfseries \abstractname\vspace{-.5em}\vspace{\z@}}
  \end{center}\begin{quotation}}{\end{quotation}\normalsize}
