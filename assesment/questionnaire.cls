%  ========================================================================  %
% 
% 	Author:	Mark G. Gray
% 		Los Alamos National Laboratory
% 	Date:	Thu Feb  4 13:41:46 MST 1999
% 
% 	Copyright (c) 1999 U. S. Department of Energy. All rights reserved.
% 
%  	$Id$
% 
%  ========================================================================  %
			  %% identification
\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{questionnaire}[1999/02/04 Questionnaire Class]

			      %% options
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{article}}
\ProcessOptions
			  %% Package loading
\RequirePackage{calc}

			       %% main
\LoadClass{article}

%% Modify the normal text area on the page.
\setlength{\hoffset}{0pt}
\setlength{\textwidth}{6.5in}
\setlength{\oddsidemargin}{0pt}
\setlength{\marginparsep}{0pt}
\setlength{\marginparwidth}{0pt}

% fix section numbering to conform to KPA.## template
\renewcommand{\thesection}
{
    KPA.
    \ifnum \value{section}<10 0\fi
    \arabic{section}
}

% Squrill away Key Process Area name for use in KPARate
\newcommand{\KPA@name}{}
\let\KPA@section=\section
\renewcommand{\section}[1]{\renewcommand{\KPA@name}{#1}\KPA@section{#1}}

% Key Process Activity environment
\newcounter{activity}[section]		% KPA counter
\newcounter{approach}[section]
\newcounter{deployment}[section]
\newcounter{result}[section]
\newcounter{NA}[section]
\newcounter{score}
\newcounter{total}[section]
\newenvironment{KPAActivity}
{
    {\bf Activity Indicators (rate your organization's capability
    on an integer scale from 0 to 10)} % header
    \begin{center}
    \begin{tabular}{|p{0.5in}|p{6.0in}|} \hline % table
}
{
    \end{tabular}
    \end{center}
    \setcounter{activity}{\value{activity}-\value{NA}}
    \setcounter{approach}{\value{approach}/\value{activity}}
    \setcounter{deployment}{\value{deployment}/\value{activity}}
    \setcounter{result}{\value{result}/\value{activity}}
    \setcounter{total}{\value{total}/\value{activity}}
    \setcounter{total}{\value{total}/3}
    The average score for \KPA@name\ is \arabic{total}.
}

% \Activity{Approach score}{Deployment score}{Results score}{Description}
\newcommand{\Activity}[4]
{
    \stepcounter{activity}
    \ifx N#1
        \stepcounter{NA} N/A
    \else
        \setcounter{score}{#1+#2+#3}
	\setcounter{total}{\value{total}+\value{score}}
        \setcounter{score}{\value{score}/3}
        \setcounter{approach}{\value{approach} + #1}
        \setcounter{deployment}{\value{deployment} + #2}
        \setcounter{result}{\value{result} + #3}
        \arabic{score}
    \fi &  \arabic{activity}. #4 \\ \hline 
} % step counter, print N/A or score, text

% Key Process Area Score
\newenvironment{KPARate}
{
    {\bf Rate your approach, deployment, and results for
    this process area on an integer scale from 0 to 10, and 
    identify evidence below that will be provided to support your
    ratings for this process area.}
    \begin{center}
    \begin{tabular}{|p{1.0in}|p{0.5in}|p{5.0in}|} \hline
    & Rating (0-10) & \multicolumn{1}{c|}{Evidence} \\ \hline
}
{
    \end{tabular}
    \end{center}
}

% \Approach{Evidence}
\newcommand{\Approach}[1]
{
    Approach to \KPA@name & \arabic{approach} & #1 \\ \hline
} 

% \Deployment{Evidence}
\newcommand{\Deployment}[1]
{
    Deployment of \KPA@name & \arabic{deployment} & #1 \\ \hline
}

% \Results{Evidence}
\newcommand{\Results}[1]
{
    Results Achieved in \KPA@name & \arabic{result} & #1 \\ \hline
}

\newenvironment{KPAComment}
{
    {\bf Provide comments on other approaches, plans,
    or applicability of this process area to your organization's
    environment.  Use a separate sheet if preferred.}
    \begin{center}
    \begin{tabular}{|p{6.50in}|} \hline
}
{
    \\ \hline
    \end{tabular}
    \end{center}
}


%  ========================================================================  %
