%%---------------------------------------------------------------------------%%
%% $Id$
%%---------------------------------------------------------------------------%%

\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{CCS4seminar}[2004/02/02 CCS-4 Seminar Class, M.W. Buksas]

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
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{seminar}}
\ProcessOptions

%% This is based on the slides class
\LoadClass{seminar}
\RequirePackage{graphicx}
\RequirePackage[dvips]{color}
\RequirePackage{pstcol}
\RequirePackage{latexmake}
\RequirePackage{fancyhdr}
\RequirePackage{fancybox}
\RequirePackage{slidesec}

% Import bug-fixes for seminar class
\input{seminar.bug}
\input{seminar.bg2}

% Set the header width to the text width
\renewcommand{\headwidth}{\textwidth}

% Eliminate the decorative lines at the header and footer
\renewcommand{\headrulewidth}{0pt}
\renewcommand{\footrulewidth}{0pt}

% Direct access to the logos used in the title page
\newcommand{\lanlLogo}{\includegraphics[width=2.6\semin]{logo.eps}}
\newcommand{\ourLogo} {\includegraphics[width=1.0\semin]{ccs4logo.eps}}

% Make the title information available with ordinary LaTeX commands
% for use elsewhere:
\newcommand{\theTitle}{\@title}
\newcommand{\theAuthor}{\@author}
\newcommand{\theDate}{\@date}

% Define variables for the default headers/footers. This allows for
% per-slide definitions.
\newcommand{\leftFootText}{\small \theslidesection--\theslideheading}
\newcommand{\rightFootText}{\theSlideOfTotal}
\newcommand{\leftHeadText}{\theAuthor}
\newcommand{\rightHeadText}{\theDate}

% Create commands which change the header/footer variables. 
\newcommand{\leftFoot} [1]{\renewcommand{\leftFootText}{#1}}
\newcommand{\rightFoot}[1]{\renewcommand{\rightFootText}{#1}}
\newcommand{\leftHead} [1]{\renewcommand{\leftHeadText}{#1}}
\newcommand{\rightHead}[1]{\renewcommand{\rightHeadText}{#1}}

\AtEndDocument{%
  \immediate\write\@auxout{\string
    \newlabel{LastSlide}{{}{\theslide}{}{}{}}}}%

\AtEndDocument{%
  \immediate\write\@auxout{\string
    \newlabel{LastPage}{{}{\thepage}{}{}{}}}}

\newcommand{\lastSlide}{\pageref{LastSlide}}
\newcommand{\theSlideOfTotal}{\theslide\mbox{ of }\lastSlide}
\newcommand{\lastPage}{\pageref{LastPage}}
\newcommand{\thePageOfTotal}{\thepage\mbox{ of }\lastPage}

% Re-define the plain style to use the fancyhdr headings.
\fancypagestyle{plain}{
  \fancyhead[l]{\leftHeadText}
  \fancyhead[c]{\classif}
  \fancyhead[r]{\rightHeadText}
  \fancyfoot[l]{\leftFootText}
  \fancyfoot[c]{\classif}
  \fancyfoot[r]{\rightFootText}
}

% Make this the stlye for slides only.
\slidepagestyle{plain}

% Re-define the empty style to include the classification
% information. This style applies to the notes and the article pages
\fancypagestyle{empty}{
  \fancyhf{}
  \fancyhead[c]{\classif}
  \fancyfoot[c]{\classif}
  \fancyfoot[r]{\thepage}
}

% Make this the default sytle for notes and article pages
\pagestyle{empty}

% Define a style for the title page which includes classification but
% omits the page number.
\fancypagestyle{title}{
  \fancyhf{}                 % Clear all headers and footers
  \fancyhead[c]{\classif}
  \fancyfoot[c]{\classif}
}

% Create a seminar title page command which adds the logos at the
% bottom of the original slide title.
\newcommand{\seminarTitle}{
  \mbox{}  \vfill
  \slidemaketitle
  \vfill
  \lanlLogo \hfill \ourLogo
  \addtocounter{slide}{-1}  % Don't count this page in the total
  \pagestyle{title}
}







