---
title: "BKP: An R Package for Beta Kernel Process Modeling"
shorttitle: "BKP"
# date: "2026.04.11"
header-includes:
  - \AtBeginDocument{\author[Jiangyan Zhao (ECNU)]{Jiangyan Zhao}}
  - \AtBeginDocument{\date[2026.07.24 (Kunming)]{2026.07.24}}
  - \renewcommand*{\bibfont}{\footnotesize}
  - \usepackage{graphicx}
  - \graphicspath{{../code/figure/}}
institute: 
  - School of Statistics, ECNU
  - jyzhao@sfs.ecnu.edu.cn
  - joint work with Kunhai Qing and Jin Xu
latex_engine: pdflatex
output:
  beamer_presentation:
    slide_level: 3
    keep_tex: true
    keep_md: true
    citation_package: natbib
    includes:
      in_header: preamble.tex
bibliography: ../paper/refs.bib
biblio-style: "chicago"
link-citations: true
classoption: "aspectratio=169"
# colorlinks: true    # 启用链接颜色
# citecolor: blue     # 设置引文链接颜色（如 blue, red, green, magenta）
# urlcolor: blue      # 设置 URL 链接颜色
# linkcolor: red      # 设置内部跳转链接颜色（如目录）
---



# Background

- Estimating a **continuous probability function** from **binary** or **binomial observations** is fundamental in statistics and machine learning.

- Applications include:

  - Binary classification \citep{MacKenzie2014, Wen2025KLR}
  
  - Probability calibration \citep{Sung2020binaryCalibration, Dimitriadis2023}
  
  - Relative abundance modeling \citep{Martin2020}
  
  - Longitudinal patient-reported outcomes \citep{Najera2018HRQoL, Najera2019PRO}
  
  - Infectious disease risk prediction \citep{Diggle2007Loa}
  
  - ...

### Classical methods & Challenges

- **Parametric method**: logistic regression is \textcolor{red}{computationally efficient} but \textcolor{blue}{limited in capturing complex nonlinear patterns}.

- **Nonparametric method**: Gaussian Process (GP) offers \textcolor{red}{flexibility and principled uncertainty quantification}, but \textcolor{blue}{approximate inference} for binary (*non-Gaussian*) likelihood increases computational cost.

### Beta Kernel Process

- **Beta Kernel Process (BKP)**: a scalable and interpretable nonparametric model for binomial probability surfaces \citep{Goetschalckx:2011, MacKenzie2014, Rolland2019BKP}

- **Key features**:
  - \textcolor{red}{Computational efficiency}: closed-form posterior inference 
  
  - \textcolor{red}{More transparency and straightforward}: beta-binomial conjugate pair
  
  - \textcolor{red}{Uncertainty quantification}: adaptive decision-marking
  
- **Dirichlet Kernel Process (DKP)**: a natural extension for categorical or multinomial data



# Summary

### Summary

- Developed the \textcolor{red}{\bf first publicly available \pkg{BKP} package} implementing BKP and DKP models

- \textcolor{red}{\bf Methodological contributions:}

  - Provided a **Bayesian interpretation** of the BKP model

  - Proposed a **flexible framework for constructing priors**, including data-adaptive options

  - Developed a **robust and efficient hyperparameter optimization strategy** using LOOCV and multi-start reparameterization.

### Future Development Directions

- **Model extensions**: Extend the BKP framework to multivariate, functional, and spatio-temporal data, as well as mixed-type covariates

- **Alternative likelihoods**: Explore likelihoods beyond the binomial family, e.g., negative binomial for over-dispersed counts, geometric for waiting-time modeling

- **Community contributions**: We invite developers to contribute via GitHub pull requests \href{https://github.com/Jiangyan-Zhao/BKP/pulls}{\texttt{https://github.com/Jiangyan-Zhao/BKP/pulls}}

### Resources

- \textcolor{blue}{\faRProject}\  \href{https://cran.r-project.org/web/packages/BKP/index.html}{\texttt{https://cran.r-project.org/web/packages/BKP/index.html}}


- \text{\hspace{0.1cm}}\textcolor{black}{\faGithub}\hspace{0.1cm} \href{https://github.com/Jiangyan-Zhao/BKP}{\texttt{https://github.com/Jiangyan-Zhao/BKP}}


- \text{\hspace{0.15cm}}\textcolor{red}{\faFilePdf}\hspace{0.15cm} \href{https://arxiv.org/abs/2508.10447}{\texttt{https://arxiv.org/abs/2508.10447}}


---


\begin{center}\includegraphics[width=0.6\linewidth]{Thanks} \end{center}
