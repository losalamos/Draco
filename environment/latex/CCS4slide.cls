%%---------------------------------------------------------------------------%%
%% CCS4Slide.cls
%%
%% OPTIONS:
%%
%%   SRD: Puts SECRET/RD or UNCLASSIFIED markers at top and botom of
%%        each page
%%
%%   All other options for the slide class should work as
%%   well. (e.g. landscape)
%%
%% COMMANDS:
%%
%%  These are active only if the SRD option was passed:
%%
%%   \SRDtext{}    Override default label (SECRET/RD in red text) for
%%                 classified slides 
%%
%%   \notSRDtext{} Override default label (UNCLASSIFIED in black text)
%%                 for unclassified slides in an otherwise SRD presentation
%%
%%   \notSRD       Removes classification markers from a single slide
%%                 when used in slide environment, or for all
%%                 subsequent slides when used outside slide
%%                 envirnoment.
%%
%%   \isSRD        Adds (restores) classification markers to a single
%%                 slide if used in slide environment, or to all
%%                 subsequent slides if used outside slide
%%                 environment.
%%
%%
%%  These options are always active:
%%
%%   \lines          Draws a slide-wide double line
%%
%%   \slidetitle{}   Puts a centered large-font title on the slide
%%                   with two lines underneath
%%
%%   \theTitle       Contents of the \title{} command.
%%
%%   \theAuthor      Contents of the \author{} command.
%%
%%   \theDate        Contents of the \date{} command. Defaults to
%%                   \today
%%
%%  The following are always active and set the contents of the
%%  left or right part of the header or footer. Default values are as
%%  given
%%
%%   \setLeftFoot{}  Default = empty
%%
%%   \setRightFoot{} Default = slide number
%%
%%   \setLeftHead{}  Default = \theAuthor
%%
%%   \setRightHead{} Default = \theDate
%%
%%   
%%
%%---------------------------------------------------------------------------%%
%% $Id$
%%---------------------------------------------------------------------------%%

\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{CCS4slide}[2004/02/02 CCS-4 Slide Class, M.W. Buksas]

\def\classif{}

%% OPTIONS
\DeclareOption{SRD}{
  \def\srd{\textcolor{red}{SECRET/RD}}
  \def\unclass{UNCLASSIFIED}
  \def\classif{\srd}
  \newcommand*{\SRDtext}[1]{\def\srd{#1}}
  \newcommand*{\notSRDtext}[1]{\def\unclass{#1}}
  \newcommand*{\notSRD}{\def\classif{\unclass}}
  \newcommand*{\isSRD}{\def\classif{\srd}}
}
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{slides}}
\ProcessOptions

%% This is based on the slides class
\LoadClass{slides}
\RequirePackage{graphicx}
\RequirePackage[dvips]{color}
\RequirePackage{pstcol}
\RequirePackage{latexmake}

\newcommand{\lines}{
  \rule[\baselineskip]{\textwidth}{0.3mm} \\[-1\baselineskip] \vspace{0.2mm}
  \rule[\baselineskip]{\textwidth}{0.3mm} \\[-2.5\baselineskip] 
}

\newcommand{\slidetitle}[1]{
  {\Large \centerline {#1}}
  \vspace{-\baselineskip}\vspace{1mm}
  \lines \par
}

\newcommand{\pageNumber}{\makebox[3em][l]{\theslide}}

\newcommand{\theTitle}{\@title}
\newcommand{\theAuthor}{\@author}
\newcommand{\theDate}{\@date}

\newcommand{\leftFoot}{}
\newcommand{\rightFoot}{\pageNumber}
\newcommand{\leftHead}{\theAuthor}
\newcommand{\rightHead}{\theDate}

\newcommand{\setLeftFoot} [1]{\renewcommand{\leftFoot}{#1}}
\newcommand{\setLeftHead} [1]{\renewcommand{\leftHead}{#1}}
\newcommand{\setRightFoot}[1]{\renewcommand{\rightFoot}{#1}}
\newcommand{\setRightHead}[1]{\renewcommand{\rightHead}{#1}}


%% headings (regular pages)
\def\ps@headings{%
  \def\ps@slide{%
    \def\@oddhead{\@mainsize 
      \makebox[0em][l]{\leftHead} \hfill \classif \hfill 
      \makebox[0em][r]{\rightHead}
    }
    \def\@oddfoot{\@mainsize 
      \makebox[0em][l]{\leftFoot} \hfill \classif \hfill
      \makebox[0em][r]{\rightFoot} 
    }
  }
}

\newcommand{\lanlLogo}{\includegraphics[width=2.6in]{logo.eps}}
\newcommand{\ourLogo}{\includegraphics[width=1.0in]{ccs4logo.eps}}

%% title page -> opening
\renewcommand{\maketitle}{%
  \parbox[t]{0in}{\mbox{} \\[-1em] \makebox[0in][l]{\lanlLogo}}
  \hfill {\small \classif} \hfill 
  \parbox[t]{0in}{\mbox{} \\[-1em] \makebox[0in][r]{\ourLogo}}
  \vspace{-0.5in}
  \par
  {\centering {\Large \theTitle} \par \theAuthor \par \theDate \par}
  \vfill
  \parbox{\textwidth}{\mbox{} \hfill {\small \classif} \hfill \mbox{}}
}    



%%---------------------------------------------------------------------------%%
%% end of CCS4slide.cls
%%---------------------------------------------------------------------------%%
