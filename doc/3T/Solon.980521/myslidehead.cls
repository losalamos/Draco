\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{myslidehead}[1997/07/01 1.0 (RMR)]

\DeclareOption*{\PassOptionsToClass{\CurrentOption}{slides}}
\ProcessOptions
\LoadClass{slides}

\newsavebox{\@slideheaderbox}
\newlength{\@slideheaderheight}
\newlength{\@slideheaderdepth}

\newcommand{\setslideheader}[1]{%
  \sbox{\@slideheaderbox}{#1}%
  \settoheight{\@slideheaderheight}{\usebox{\@slideheaderbox}}%
  \settodepth{\@slideheaderdepth}{\usebox{\@slideheaderbox}}%
  %
  \addtolength{\textheight}{\headheight}%
  \setlength{\headheight}{\@slideheaderheight}%
  %\addtolength{\headheight}{\@slideheaderdepth}%
  \addtolength{\textheight}{-\headheight}%
  \setlength{\topmargin}{0in}%
  \setlength{\headsep}{0in}%
}

\newsavebox{\@slidefooterbox}
\newcommand{\setslidefooter}[1]{%
  \sbox{\@slidefooterbox}{#1}%
}

\def\ps@headings{%
  \def\ps@slide{%
    \def\@oddfoot{\@mainsize \usebox{\@slidefooterbox} %
      \hfil\hb@xt@3em{\theslide \hss}}%
    \def\@oddhead{\@mainsize \hfil \usebox{\@slideheaderbox} \hfil }%
    \def\@evenfoot{\@oddfoot}%
    \def\@evenhead{\@oddhead}}

  \def\ps@overlay{%
    \def\@oddfoot{\@mainsize \mbox{}\hfil\hb@xt@3em{\theoverlay\hss}}%
    \def\@oddhead{}%
    \def\@evenfoot{\@mainsize \mbox{}\hfil\hb@xt@3em{\theoverlay\hss}}%
    \def\@evenhead{}}

  \def\ps@note{%
    \def\@oddfoot{%
      \@mainsize
      \if@clock \fbox{\large \@arabic\c@minutes\space min}\else \hbox{}\fi
      \hfil\thenote}%
    \def\@oddhead{}%
    \def\@evenfoot{%
      \@mainsize
      \if@clock \fbox{\large \@arabic\c@minutes\space min}\else \hbox{}\fi
      \hfil\thenote}%
    \def\@evenhead{}}}
