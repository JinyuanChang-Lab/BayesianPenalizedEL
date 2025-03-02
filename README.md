## Official Implementations of `Bayesian penalized empirical likelihood and Markov Chain Monte Carlo sampling`

## Introduction

In this study, we introduce a novel methodological framework called Bayesian penalized empirical likelihood (BPEL), designed to address the computational challenges inherent in empirical likelihood (EL) approaches. Our approach has two primary objectives: (i) to enhance the inherent flexibility of EL in accommodating diverse model conditions, and (ii) to facilitate the use of well-established Markov Chain Monte Carlo sampling schemes as a convenient alternative to the complex optimization typically required for statistical inference using EL. To achieve the first objective, we propose a penalized approach that regularizes the Lagrange multipliers, significantly reducing the dimensionality of the problem while accommodating a comprehensive set of model conditions. For the second objective, our study designs and thoroughly investigates two popular sampling schemes within the BPEL context. We demonstrate that the BPEL framework is highly flexible and efficient, enhancing the adaptability and practicality of EL methods. Our study highlights the practical advantages of using sampling techniques over traditional optimization methods for EL problems, showing rapid convergence to the global optima of posterior distributions and ensuring the effective resolution of complex statistical inference challenges.

## Data and Files

* **Data Generation Process**: The folder "Data Generation Process" includes all the Matlab codes for the data generation process in the numerical studies, and the associated simulated data can be accessed via the [Google Drive](https://drive.google.com/file/d/1OmH2cOvSc7XZSP7VywsrEzI9-x0B9-61/view).

* **data_SES.mat**: The "data_SES.mat" in folder "Sampling Efficiency and Stability" is the data file utilized for Section 3.2 in the main paper.

* **BInit.mat**: The "BInit.mat" in folder "Real Data Analysis" is the initial points used for the real data analysis in this study.

* **trade_data.mat**: The "trade_data.mat" in [release 1.0.0](https://github.com/JinyuanChang-Lab/BayesianPenalizedEL/releases/tag/1.0.0) is the data file utilized for the real data analysis in this study, as sourced from [Shi (2016, JoE)](https://github.com/zhentaoshi/REL).

## Codes

We have adhered to the section titles and provided the main codes used for numerical studies and real data analysis in the main paper. Furthermore, the main codes for the numerical studies in the supplementary material are included in folder "Additional Numerical Results".

