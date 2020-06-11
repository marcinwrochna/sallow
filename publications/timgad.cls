%%% LaTeX class "timgad" for modern-looking preprints
%%% author: Marcin Wrochna
%%% version: 0.03
%%% licence: public domain
%%% http://mimuw.edu.pl/~mw290715/timgad
%%%
%%% Based on KOMA-script "scrartcl" class
%%% 

\NeedsTeXFormat{LaTeX2e}[1995/12/01]
\ProvidesClass{timgad}[2015/09/29 LaTeX class "timgad" for modern-looking preprints]

%----- options ----------------------------------------------------------------
% Fallback handler for all other options: pass them to scrartcl.
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{scrartcl}}
% Default options:
\ExecuteOptions{12pt,DIV=13,numbers=noendperiod} %,version=3.16
% End options (process them).
\ProcessOptions\relax

\LoadClass{scrartcl}

% Basic page geometry is ok with DIV=13 option, but we prefer to set it by hand.
% Setting fixed totalwidth and totalheight in pt will allow the breaks to be exactly the same
% on A4, letter and legal, and in particular on arXiv (which messes with it).
\RequirePackage[totalwidth=460pt, totalheight=680pt]{geometry}
%\RequirePackage[vmarginratio={4:5},vscale=0.79,hscale=0.77]{geometry}
% Add basic utf8 support
\RequirePackage[utf8]{inputenc}
% Make word-wrapping and micro-spaces like kerning a little nicer
\RequirePackage{microtype}
\frenchspacing

%A macro tool.
\usepackage{letltxmacro}


%----- colors -----------------------------------------------------------------
%Modern version of color package
\RequirePackage{xcolor}
%That touch of color used everywhere
\definecolor{timgadBgColor}{RGB}{185,222,152}
%Color of header text
\definecolor{dispositionColor}{RGB}{0,0,0}

%----- maketitle --------------------------------------------------------------
%A nice package for defining many authors with many affiliations.
\RequirePackage{authblk}

\renewcommand\Affilfont{\small}
\setlength{\affilsep}{0em}

% Define \title with optional argument for thanks, handle pdf.
\LetLtxMacro{\old@title}{\title}
\def\title{\@ifnextchar[{\@dotitlewiththanks}{\@dotitle}}
\gdef\@dotitle#1{\old@title{#1}}
\gdef\@dotitlewiththanks[#1]#2{\old@title{#2\texorpdfstring{\thanks{#1}}{}}}

%Define \author to additionally add to a comma-separated list \string@author.
\def\string@author{\@gobbletwo}
\LetLtxMacro{\old@author}{\author}
\renewcommand\author[2][]{%
	\g@addto@macro\string@author{, \space #2}%
	\old@author[#1]{#2}%
}

%Hide \thanks in pdf title/author meta-data.
\LetLtxMacro{\old@thanks}{\thanks}
\renewcommand\thanks[1]{\texorpdfstring{\old@thanks{#1} }{}}

%Declare \keywords,\pdftitle,\pdfauthor variables to be optionally set by author.
\def\keywords#1{\gdef\@keywords{#1}}
\let\@pdftitle\@empty
\let\@pdfauthor\@empty
\def\pdftitle#1{\gdef\@pdftitle{#1}}
\def\pdfauthor#1{\gdef\@pdfauthor{#1}}

%Define the 'title' part (excluding abstract).
\def\@maketitle{%
	\newpage
	\noindent\colorbox{timgadBgColor}{\hspace*{0.6\textwidth}}\\
	{\raggedright%
		{\LARGE\bfseries\@title\par\vspace{0.4em}}%
		\vskip 0.6em%
		{\@author\par\vspace{0.4em}}%
		{\small\@date}%
	}%
	\vskip 1em}


%----- section style ----------------------------------------------------------
%Use the color for header text.
\addtokomafont{disposition}{\color{dispositionColor}}
%Define the colored rectangle for section number.
\renewcommand*{\othersectionlevelsformat}[3]{%
		\colorbox{timgadBgColor}{%
			\color{dispositionColor}\kern0.1em#3\kern0.1em%
		}\enskip%
}
%With v>3.17 of KOMA-script you'd renew \sectionformat
%	for all levels instead (maybe  \dispositionformat?)
%Without KOMA-script you'd redefine \@seccntformat, probably.

% That's for unnamed \section* variant
\LetLtxMacro\old@section{\section}
\renewcommand\section{\@ifstar{%
	\@startsection{section}{\sectionnumdepth}{\z@}%
	{-3.5ex \@plus -1ex \@minus -.2ex}%
	{2.3ex \@plus.2ex}%
	{\ifnum \scr@compatibility>\@nameuse{scr@v@2.96}\relax
		\setlength{\parfillskip}{\z@ plus 1fil}\fi
		\raggedsection\normalfont\sectfont\nobreak\size@section {\color{timgadBgColor}\rule[-1pt]{0.5em}{2ex}}\kern0.1em}*%
}{\old@section}}

%----- pdf and hyperref -------------------------------------------------------
% Has to be before hyperref to avoid warnings.
\RequirePackage{float}
%Make pdf include nice meta-data
\RequirePackage{hyperref}
%Color links
\hypersetup{
    colorlinks,
    linkcolor={green!40!black},
    citecolor={blue!50!black},
    urlcolor={blue!70!black},
%	linkbordercolor={red!50!black},
%	pdfborderstyle={/S/U/W 1}% border style will be underline of width 1pt
	pdfdisplaydoctitle=true,
}
\AtBeginDocument{\hypersetup{
    pdftitle={\ifx\@pdftitle\empty\@title\else\@pdftitle\fi},
    pdfauthor={\ifx\@pdfauthor\empty\string@author\else\@pdfauthor\fi},
    pdfkeywords={\@keywords},
    %pdfcreator, pdfproducer, pdfsubject,
}}
%Define sectioning meta-data for pdf-readers
\RequirePackage[numbered]{bookmark}










\RequirePackage{amssymb}
\RequirePackage{amsmath}
\RequirePackage{amsfonts}
\RequirePackage{amsthm}

%----- theorem style ----------------------------------------------------------
%Our super-nice thm number rectangle
\newtheoremstyle{bgnumberbefore}
{\topsep}{\topsep}{\itshape}{}{\sffamily\bfseries\color{dispositionColor}}{}{5pt plus 1pt minus 1pt}
{\colorbox{timgadBgColor}{\color{dispositionColor}\kern0.01em\thmnumber{#2}\kern0.01em} \thmname{#1.}\thmnote{ \normalfont(#3) }}
\theoremstyle{bgnumberbefore}

\newtheorem{theorem}{Theorem}[section]
\newtheorem{lemma}[theorem]{Lemma}
\newtheorem{proposition}[theorem]{Proposition}
\newtheorem{corollary}[theorem]{Corollary}
\newtheorem{observation}[theorem]{Observation}
\newtheorem{remark}[theorem]{Remark}
\newtheorem{statement}[theorem]{Statement}
\newtheorem{definition}[theorem]{Definition}
\newtheorem{claim}{Claim}[theorem]
\newtheorem{fact}[theorem]{Fact}
\newtheorem{conjecture}[theorem]{Conjecture}

%Allow restating theorems.
\newcommand{\newrep@theorem}[1]{%
	\newenvironment{rep#1}[2][\@gobbletwo]{%
		\begingroup%
		\def\thetheorem{##2}%
		\begin{#1}[##1, {}{\sffamily restated}]%
	}{%
		\end{#1}%
		\addtocounter{theorem}{-1}%
		\endgroup%
	}%
}
\newrep@theorem{theorem}
\newrep@theorem{lemma}
\newrep@theorem{proposition}
\newrep@theorem{corollary}
\newrep@theorem{fact}
\newrep@theorem{observation}
\newrep@theorem{remark}
\newrep@theorem{statement}
\newrep@theorem{definition}

%----- proof style ------------------------------------------------------------
% mdframed style with corners:
\usepackage[framemethod=tikz]{mdframed}
\mdfdefinestyle{proofframe}{
	skipabove         = 2pt minus 2pt plus 2pt,
	skipbelow         = 4pt plus 1pt minus 2pt,
	leftmargin        = 0pt ,
	rightmargin       = 0pt ,
	innermargin       = 0pt ,
	innertopmargin    = .3em ,
	innerleftmargin   = .5em ,
	innerrightmargin  = 0pt ,
	innerbottommargin = .3em ,
	hidealllines      = true ,
	singleextra       = {
		\draw (O) -- ++(0,1.4em) (O) -- ++(.675em,0) ;
		\draw (P-|O) -- ++(0,-1.4em) (P-|O) -- ++(.675em,0) ;
	},
	firstextra        = {
		\draw (P-|O) -- ++(0,-1.675em) (P-|O) -- ++(.675em,0) ;
	},
	secondextra       = {
		\draw (O) -- ++(0,1.675em) (O) -- ++(.675em,0) ;
	},
}

\renewenvironment{proof}[1][\proofname]{\par
	\vspace{-3pt}
	\begin{mdframed}[style=proofframe]%	
	\pushQED{\qed}%
	\normalfont \topsep6\p@\@plus6\p@\relax
	\trivlist\item[\hskip\labelsep\rmfamily\itshape#1\@addpunct{.}]\ignorespaces
}{%
	\popQED\endtrivlist%
	\end{mdframed}\@endpefalse%
}

%----- figure -----------------------------------------------------------------
\setkomafont{captionlabel}{\sffamily\bfseries\color{dispositionColor}}
\renewcommand*{\figureformat}{{\color{timgadBgColor}\rule[-1pt]{0.45em}{2	ex}}\kern 0.15em\figurename~\thefigure}
\renewcommand*{\captionformat}{\ } %removes trailing ":"
\setcapindent{0.8em}
\setcapmargin{1em}
\setlength{\belowcaptionskip}{-0.4\baselineskip}


%----- misc -------------------------------------------------------------------
\RequirePackage{enumitem}
\setlist{topsep=5pt,itemsep=3pt,parsep=3pt}



\endinput
