%% LIGODOC_v2.CLS
%%
%% LIGO Document Class
%%
%% All options are passed to the underlying article class implementation
%% 
%% The following macros are defined:
%% 
%%   \title{title}
%%      document title
%%   \author{authors}
%%      list of authors
%%   \institution{institution}
%%      list of author's institutions
%%   \date{date}
%%      publication date
%%   \ligodccnumber{T}{YY}{NNNN}{VV}{G}
%%      LIGO document control number
%%        T     document type identifier
%%        YY    two digit year
%%        NNNN  four digit serial number
%%        VV    two digit version number
%%        G     working group identifier
%%   \ligodistribution{group}
%%      includes ``Distribution of this document: <group>'' on the title page
%%   \ligonote{note}
%%      includes the specified note on the title page
%%   \ligointernal
%%      sets note to ``This is an internal working note of the LIGO project''
%%   \ligopublished{journal}
%%      sets note to ``To be published in <journal>''
%%   \ligosubmitted{journal}
%%      sets note to ``Submitted for publication in <journal>''
%%   \ligodraft
%%      includes the label ``Draft'' on every page
%%   \ligorestricted{group}
%%      includes the label ``Circulation resitricted to <group>'' on every page
%%   \ligosites
%%      include LIGO site information on the title page
%%   \cvsid{$Id: ligodoc.cls,v 1.1 2008/02/24 07:11:08 rana Exp $}
%%      include the file's CVS ID on the title page
%% 
%% Shourov K. Chatterji
%% shourov@ligo.mit.edu
%% 2003-Jun-18
%%
%% edited by Rana, July 2006

%---------- LATEX VERSION ----------
\NeedsTeXFormat{LaTeX2e}

%---------- CLASS NAME AND REVISION ----------
\ProvidesClass{ligodoc_v2}[2006/07/17 LIGO document class]

%---------- PROCESS OPTIONS ----------
\DeclareOption{letterpaper}{\PassOptionsToClass{\CurrentOption}{article}}
\DeclareOption{12pt}{\PassOptionsToClass{\CurrentOption}{article}}
\DeclareOption{oneside}{\PassOptionsToClass{\CurrentOption}{article}}
\DeclareOption{onecolumn}{\PassOptionsToClass{\CurrentOption}{article}}
\DeclareOption{final}{\PassOptionsToClass{\CurrentOption}{article}}
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{article}}
\ExecuteOptions{letterpaper,12pt,oneside,onecolumn,final}
\ProcessOptions\relax

%---------- LOAD PARENT CLASS ----------
\LoadClass{article}

%---------- LOAD REQUIRED PACKAGES ----------
\RequirePackage{latexsym}
\RequirePackage{amsmath}
\RequirePackage{amsfonts}
\RequirePackage{amssymb}
\RequirePackage{graphicx}
\RequirePackage[usenames,dvipsnames]{color}
\RequirePackage{hyperref}
\RequirePackage{ifthen}
\RequirePackage{array}
\RequirePackage{calc}
\RequirePackage{vmargin}

%---------- SET MARGINS, INDENTATION, AND SPACING ----------
\setpapersize[portrait]{USletter}
\setmarginsrb{1in}{0.5in}{1in}{0.5in}{0.3in}{0.2in}{0.3pt}{0.2in}
\setlength{\parindent}{0pt}
\setlength{\parskip}{0.5\baselineskip}



%---------- SET PAGE COLOR ----------
 \pagecolor{white}

%---------- DEFAULT LIGO DOCUMENT PROPERTIES ----------
\def\@ligodcctype{\relax}
\def\@ligodccyear{\relax}
\def\@ligodccserial{\relax}
\def\@ligodccversion{\relax}
\def\@ligodccgroup{\relax}
\def\@ligodccnumber{\relax}
\def\@ligodistribution{\relax}
\def\@ligodraft{\relax}
\def\@ligorestricted{\relax}
\def\@ligonote{\relax}
\def\@ligosites{\relax}
\def\@institution{\relax}
\def\@cvsid{\relax}

%---------- MACRO TO PARSE CVS ID TAG ----------
\def\setcvsid#1{\def\next##1#1{\def\@cvsid{\mbox{%
  CVS ##1 --- Processed with \LaTeX\ on \today}}}\next}

%---------- COMMANDS TO SET LIGO DOCUMENT PROPERTIES ----------
\newcommand{\ligodccnumber}[4]{% 
  \def\@ligodcctype{#1}%
  \def\@ligodccyear{#2}%
  \def\@ligodccserial{#3}%
  \def\@ligodccversion{#4}%
  %\def\@ligodccgroup{#5}% JM edit for new-style dcc numbers
  %\def\@ligodccnumber{LIGO-#1#2#3-#4-#5}% JM edit for new-style dcc numbers
  \def\@ligodccnumber{LIGO-#1#2#3-#4}%
  \ifthenelse{\equal\@ligodcctype{C}}{%
    \def\@ligodoctype{Contractual Note}%
  }{%
  \ifthenelse{\equal\@ligodcctype{D}}{%
    \def\@ligodoctype{Drawing}%
  }{%
  \ifthenelse{\equal\@ligodcctype{E}}{%
    \def\@ligodoctype{Engineering Note}%
  }{%
  \ifthenelse{\equal\@ligodcctype{G}}{%
    \def\@ligodoctype{Presentation}%
  }{%
  \ifthenelse{\equal\@ligodcctype{L}}{%
    \def\@ligodoctype{Correspondence}%
  }{%
  \ifthenelse{\equal\@ligodcctype{M}}{%
    \def\@ligodoctype{Administrative Note}%
  }{%
  \ifthenelse{\equal\@ligodcctype{P}}{%
    \def\@ligodoctype{Publication}%
  }{%
  \ifthenelse{\equal\@ligodcctype{T}}{%
    \def\@ligodoctype{Technical Note}%
  }{%
    \def\@ligodoctype{}%
  }}}}}}}}%
}
\newcommand{\ligodraft}{\def\@ligodraft{Draft}}
\newcommand{\ligorestricted}[1]{\def\@ligorestricted{%
  Circulation restricted to #1}}
\newcommand{\ligodistribution}[1]{\def\@ligodistribution{%
  \textit{Distribution of this document:}\vspace{1ex}\\#1}}
\newcommand{\ligonote}[1]{\def\@ligonote{#1}}
\newcommand{\ligointernal}{\def\@ligonote{%
  This is an internal working\\note of the LIGO project}}
\newcommand{\ligopublished}[1]{\def\@ligonote{%
  To be published in\\#1}}
\newcommand{\ligosubmitted}[1]{\def\@ligonote{%
  Submitted for publication in\\#1}}
\newcommand{\ligosites}{\def\@ligosites{true}}
\newcommand{\institution}[1]{\def\@institution{\\#1}}
\newcommand{\cvsid}[1]{\setcvsid#1}

% Colorized Section,subsSection,& subsub
\renewcommand\section{\@startsection {section}{1}{\z@}%
                                   {-3.5ex \@plus -1ex \@minus -.2ex}%
                                   {1.4ex \@plus.2ex}%
                                   {\color{red}\normalfont\Large\bfseries}}
\renewcommand{\subsection}{\@startsection{subsection}{2}{\z@}%
                                     {-3.25ex\@plus -1ex \@minus -.2ex}%
                                     {0.5ex \@plus .2ex}%
                                     {\color{blue}\normalfont\bfseries}}
\renewcommand{\subsubsection}{\@startsection{subsubsection}{3}{\z@}%
                                     {-3.25ex\@plus -1ex \@minus -.2ex}%
                                     {0.2ex \@plus .2ex}%
                                     {\color[named]{OliveGreen}\normalfont\bfseries}}


%---------- REDEFINE DATE FORMAT ----------
\renewcommand{\today}{%
  \number\year/\ifnum\month<10 0\fi\number\month/\ifnum\day<10 0\fi\number\day}

%---------- ? ----------
\newcommand{\PreserveBackslash}[1]{\let\temp=\\#1\let\\=\temp}

%---------- LABEL LAST PAGE ----------
\AtEndDocument{\label{ligodoc-lastpage}}

%---------- DEFINE FOOTERS AND HEADERS ----------
\def\ps@headings{
  \def\@evenfoot{\vbox to 0in{%
     \centering\mbox{page \thepage{}}%
     }}
  \def\@oddfoot{\vbox to 0in{%
     \centering\mbox{page \thepage{}}%
     }}
  \def\@oddhead{\vbox to 0in{%
     \vskip -0.3in%
     \centering\mbox{\@ligodccnumber}\rule{0pt}{0ex}\\
     \vspace{0.1in}
     \centering\mbox{\color{red} \@ligodraft%
     \ifthenelse{\equal\@ligodraft{\relax}}{}{%
       \ifthenelse{\equal\@ligorestricted{\relax}}{}{:\ }}%
     \@ligorestricted}%
     }}
  \def\@evenhead{\vbox to 0in{%
     \vskip -0.3in%
     \centering\mbox{\@ligodccnumber}\rule{0pt}{0ex}\\
     \vspace{0.1in}
     \centering\mbox{\color{red} \@ligodraft%
     \ifthenelse{\equal\@ligodraft{\relax}}{}{%
       \ifthenelse{\equal\@ligorestricted{\relax}}{}{:\ }}%
     \@ligorestricted}%
     }}
}

%---------- ENABLE HEADERS AND FOOTERS ----------
\pagestyle{headings}


%---------- DEFINE TITLEPAGE ----------
\renewcommand\maketitle{%
  \thispagestyle{titleheadings}%
  \begin{center}
    \fontsize{12pt}{14pt}\selectfont%
    \mbox{%
    LASER INTERFEROMETER GRAVITATIONAL WAVE OBSERVATORY}\\
    - LIGO -\\
    CALIFORNIA INSTITUTE OF TECHNOLOGY\\
    MASSACHUSETTS INSTITUTE OF TECHNOLOGY\\[0.75in]
    \begin{tabular*}{5.5in}%
      {|*{1}{>{\PreserveBackslash\centering\hspace{0pt}}%
             p{1.5in-\tabcolsep*2-\arrayrulewidth*2/3}}%
        *{1}{>{\PreserveBackslash\centering\hspace{0pt}}%
             p{2.5in-\tabcolsep*2-\arrayrulewidth*2/3}}%
        *{1}{>{\PreserveBackslash\centering\hspace{0pt}}%
             p{1.5in-\tabcolsep*2-\arrayrulewidth*2/3}}|}
      \hline
      \textbf{\@ligodoctype} &
      \textbf{\@ligodccnumber} &
      \@date\rule[-2ex]{0pt}{5ex}\\
      \hline
      \multicolumn{3}{|p{5.5in-\tabcolsep*2-\arrayrulewidth*2}|}%
        {\rule{0pt}{5ex}\centering\parbox[t]{5in}%
          {\centering\fontsize{22pt}{16pt}\selectfont%
           \bfseries\@title\rule[-2ex]{0pt}{0pt}}}\\
      \hline
      \multicolumn{3}{|p{5.5in-\tabcolsep*2-\arrayrulewidth*2}|}%
        {\rule{0pt}{3ex}\centering\@author\textit{\@institution}%
         \rule[-2ex]{0pt}{0ex}}\\
      \hline
    \end{tabular*}\\
    \vspace{0.5in}
    \@ligodistribution\rule{0pt}{0ex}\\
    \vspace{0.25in}
    {\color{red}\fontsize{14pt}{14pt}\selectfont\@ligodraft\rule{0pt}{0ex}}\\
    {\color{red}\fontsize{14pt}{14pt}\selectfont\@ligorestricted\rule{0pt}{0ex}}\\
    \vspace{0.25in}
    \@ligonote\rule{0pt}{0ex}\\
    \vfill
    \begin{tabular}{c@{\qquad}c}
      \textbf{California Institute of Technology} &
      \textbf{Massachusetts Institute of Technology} \\
      \textbf{LIGO Project, MS 18-34} &
      \textbf{LIGO Project, Room NW17-161} \\
      \textbf{Pasadena, CA 91125} &
      \textbf{Cambridge, MA 02139} \\
      Phone (626) 395-2129 & Phone (617) 253-4824 \\
      Fax (626) 304-9834 & Fax (617) 253-7014 \\
      E-mail: info@ligo.caltech.edu & E-mail: info@ligo.mit.edu \\
        \vspace{0.3in} \\
        \textbf{LIGO Hanford Observatory} &
        \textbf{LIGO Livingston Observatory} \\
        \textbf{Route 10, Mile Marker 2} &
        \textbf{19100 LIGO Lane} \\
        \textbf{Richland, WA 99352} &
        \textbf{Livingston, LA 70754} \\
        Phone (509) 372-8106 & Phone (225) 686-3100 \\
        Fax (509) 372-8137 & Fax (225) 686-7189 \\
        E-mail: info@ligo.caltech.edu & E-mail: info@ligo.caltech.edu \\
      \vspace{0.3in}
    \end{tabular}\\[1ex]
    \href{http://www.ligo.caltech.edu/}{\tt http://www.ligo.caltech.edu/}\\
  \end{center}
  \normalsize
  \newpage
}
\def\ps@titleheadings{%
  \let\@evenhead\@empty\let\@oddhead\@empty
  \def\@oddfoot{\hfil\mbox{\fontsize{9pt}{11pt}\selectfont\@cvsid}\hfil}%
  \def\@evenfoot{\hfil\mbox{\fontsize{9pt}{11pt}\selectfont\@cvsid}\hfil}%
}

%---------- AUTOMATICALLY GENERATE TITLE PAGE ----------

% Didn't work with colortbl. Need to insert manually.

% \AtBeginDocument{%
%  \maketitle
%   \setcounter{page}{1}
% }
