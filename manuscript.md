---
title: "astroBayes!"
author:
  - Robin B. Trayler
  - Stephen R. Meyers
  - Mark D. Schmitz
bibliography: /Users/robintrayler/Zotero/ref_library.bib
csl: /Users/robintrayler/Zotero/styles/earth-and-planetary-science-letters.csl
mainfont: "Baskerville"
fontsize: 12pt
geometry: margin=1.0in
tblPrefix: Table
figPrefix: Figure
secPrefix: Section
link-citations: true
indent: true
header-includes:
    - \usepackage{lineno}
    - \linenumbers
    - \usepackage{setspace}
    - \doublespacing
---

<!-- pandoc -s -o manuscript.pdf --pdf-engine=xelatex --filter pandoc-crossref --citeproc --number-sections manuscript.md --> 

# Introduction
Developing precise and accurate models that relate stratigraphic position to absolute age (e.g., age-depth models) is a crucial step in interpreting the rate and tempo of geologic and climatological processes [@blaauw2012; @parnell2011]. Constructing these chronologies relies on a variety of geochronologic information. Various radioisotopic techniques (e.g., ^40^Ar/^39^Ar, U-Pb) allow the age of discrete stratigraphic points to be determined. 


# Statistical Methods 

$$P(parameters~|data) = \frac{P(data~|parameters)}{P(data)} \times P(parameters)$$


In our case, our data takes two forms. First, our cyclostratigraphic proxy record, which consists measurements $[d_1, d_2, ... d_i]$ where $i$ is the stratigraphic position of each measurement. We assume that cyclic signals in $d$ is derived from orbital forcing and the record can therefore be tuned. 

\newpage

# References {.unnumbered}
:::{#refs}
:::

