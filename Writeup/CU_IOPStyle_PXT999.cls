
\NeedsTeXFormat{LaTeX2e}[1995/12/01]
\ProvidesClass{CU_IOPStyle_PXT999}[Cardiff University IOP style]

\DeclareOption*{\PassOptionsToClass{\CurrentOption}{article}} 
\ProcessOptions\relax
\LoadClass[10pt,twocolumn]{scrartcl}


% colors
\RequirePackage{xcolor}

% font settings
%\RequirePackage[adobe-utopia]{mathdesign}
%\RequirePackage[T1]{fontenc}
%\RequirePackage[latin1]{inputenc}


% page settings
\RequirePackage[a4paper,margin=1.6cm,top=2.6cm, bottom=2.6cm]{geometry}


% header and footer settings
\RequirePackage{fancyhdr}
\renewcommand{\headrulewidth}{1pt}
\setlength{\headheight}{30pt}
\setlength{\footskip}{40pt}
\pagestyle{fancy}
\lhead{\color{gray} Cardiff University Publishing}
\rhead{\color{gray} School of Physics and Astronomy}
\lfoot{}
\cfoot{\vskip1.5em \color{gray}\thepage}
\rfoot{}


% language settings
\RequirePackage[american]{babel}


% captions settings
\RequirePackage{caption}
\captionsetup{labelfont={small,bf}}
\setlength\belowcaptionskip{3pt}


% graphics settings
\RequirePackage{graphicx}


% tweaking the footnotes
\RequirePackage[hang,flushmargin]{footmisc} 


% optimized typography
\RequirePackage{ifpdf}
\ifpdf
  \RequirePackage[protrusion=true,expansion=true]{microtype}
\fi


% new commands to pass authors and affiliations
\def\addtok#1#2{\global#1\expandafter{\the#1#2}}
\newcounter{n@authors}
\newtoks\@authors
\newtoks\@pdfauthors
\newtoks\@shortauthors

\renewcommand{\author}[3]{%
	\stepcounter{n@authors}%
  	\ifnum \then@authors > 1 \addtok\@authors{~$\cdot$~} \addtok\@pdfauthors{, } \addtok\@shortauthors{~$\cdot$~} \fi%
  	\addtok\@authors{#2~#1$^{\,#3}$}%
  	\addtok\@pdfauthors{#2~#1}%
  	\addtok\@shortauthors{#1}%
}

\newcounter{n@affiliations}
\newtoks\@affiliations
\newcommand{\affiliation}[2]{%
	\stepcounter{n@affiliations}%
  	\ifnum \then@affiliations > 1 \addtok\@affiliations{\newline} \fi%
  	\addtok\@affiliations{$^{#1}$#2}%
}


% title page
\renewcommand{\maketitle}{
  	\vskip1em%
  	{\flushleft\huge\bfseries\@title\par}%
  	\vskip1em%
  	{\flushleft\large\the\@authors\par}%
  	{\flushleft\the\@affiliations}%
  	\rhead{}%
  	\hypersetup{pdftitle=\@title}%
  	\hypersetup{pdfauthor=\the\@pdfauthors}%
}


% abstract 
\renewenvironment{abstract}{
	\centering
  	\vskip1em
  	\begin{minipage}[t]{\textwidth}
  		\noindent%
  		{\flushleft\bfseries\abstractname\vskip0.25em}%
	}{
  	\end{minipage}\par%
  	\vskip2em
}


% bibliography settings
%\RequirePackage[square,numbers]{natbib}
\bibliographystyle{ieeetr}


% font sizes and vertical spacings
\setlength\parindent{1.5em}
\setlength\partopsep{2\p@ \@plus 1\p@ \@minus 1\p@}
\def\@listi{\leftmargin\leftmargini
            \parsep 2\p@ \@plus\p@ \@minus\p@
            \topsep 4\p@ \@plus\p@ \@minus2\p@
            \itemsep2\p@ \@plus\p@ \@minus\p@}
\let\@listI\@listi
\@listi
\def\@listii {\leftmargin\leftmarginii
              \labelwidth\leftmarginii
              \advance\labelwidth-\labelsep
              \topsep    2\p@ \@plus\p@\@minus\p@
              \parsep    \z@
              \partopsep \p@ \@plus\z@ \@minus\p@
              \itemsep   \parsep}
\def\@listiiii {\leftmargin\leftmarginiii
              \labelwidth\leftmarginiii
              \advance\labelwidth-\labelsep}
\def\@listiv {\leftmargin\leftmarginiv
              \labelwidth\leftmarginiv
              \advance\labelwidth-\labelsep}
\def\@listv  {\leftmargin\leftmarginv
              \labelwidth\leftmarginv
              \advance\labelwidth-\labelsep}
\def\@listvi {\leftmargin\leftmarginvi
              \labelwidth\leftmarginvi
              \advance\labelwidth-\labelsep}

\renewcommand\thesection{\arabic{section}}  
\setkomafont{section}{\normalfont\bfseries\color{black}}
%\renewcommand\thesubsection{\thesection.\Roman{subsection}}
\setkomafont{subsection}{\normalfont\normalsize\itshape\color{black}}
\setkomafont{subsubsection}{\normalfont\normalsize\itshape\color{black}}
\setkomafont{paragraph}{\normalfont\normalsize\bfseries\color{black}}
\setkomafont{subparagraph}{\normalfont\normalsize\bfseries\color{black}}
  

% hyperref settings
\RequirePackage{hyperref}
\hypersetup{pdfpagelayout=SinglePage,pdfdisplaydoctitle}


