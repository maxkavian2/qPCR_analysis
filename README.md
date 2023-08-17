# qPCR_analysis
qPCR analysis pipelines developed at A. Sanz laboratory, School of Molecular Biosciences, University of Glasgow (July 2023).

For any questions regarding this code, please send a message to: maxsanara@gmail.com.

# Basic requirements
* Quantum Studio 3 datasets of target gene and reference gene (Tab-separated text-quoted CSV).
* Bash based terminal.
* R (including the xtable package).

# Introduction
The following pipelines have been designed to work with a CSV export of the Excel file format of the results provided by QuantumStudio 3 software for real time polymerase chain reaction (qPCR) analysis. It is mandatory that the source CSVs are tab separated and their text fields quoted.

Two different analysis are included that can be executed independently. In both cases the procedure depends on the existence of a two datasets, one for the a target gene whose expression is analyzed and one for the reference gene (as a housekeeping gene) whose expression does not change presumably over the course of the experiment. Each data set contains meaningful experimental categories, such as different types of controls or mutants under study. Each category must be represented by different samples (i.e. _biological replicas_) which in turn are split into more replicas (i.e. _technical replicas_). Technical replicas try to cope with instrumental errors, e.g. pipetting, fluorescence noise, etc. and could be treated separatedly.

The two aformentioned analysis can be defined as follows: 

* **Through standard curve of dilution factors.** This method relies on the estimated dilution factor of the target mRNA through a standard curve of known *dilution factors*. Technical replicas of each biological sample are first assessed, kept or discarded according to their standard deviation[^1], before they are reduced to a central estimation such as their averages. These averages are regarded as representative biological measurements gathered under one category (e.g. one type of control). These averages are later normalised over the partnered average from the reference gene dataset. Two normalizations are provided: dilution factors ratios and differences of C_t (i.e. Delta C_t). Finally every pair of categories is compared by both Kolmogorov-Smirnov and t-Student tests.
* **DeltaDetal Ct method with Montecarlo correction** In this case no distinction is made between technical replicas and biological samples as long as the number of technical replicas per sample is the same. Measurements associated to one category are put together and each category mean distribution estimated, corresponding to the parametrization of a t-Student distribution. Then distribution of Delta C_t values are determined between the target gene and the reference gene, and the distribution of differentials between Delta C_t ascertained (i.e. DeltaDelta C_t) afterwards. Since these distrubutions have no analytic expression, a Montecarlo approach is undertaken over a large number of points (e.g. 10^4). Finally the relative level of expression is established by the distribution of 2^{-DeltaDelta C_t}, where 1 means no change, >1 means upregulation and <1 downregulation.

The results are given in PDF, SVG (charts) and LaTex compatible codes. Blank charts with no label are also included for TikZ annotation in LaTeX code. Tables results are also rendered for further analysis in R.

[^1] : Quantum Studio 3 includes two different informations, i.e. *flags*, to this regard: high _standard deviation_ and _outlier_. When t
