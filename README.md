# qPCR analysis for Quantum Studio 3 Data Formats.
qPCR analysis pipelines developed at A. Sanz laboratory, School of Molecular Biosciences, University of Glasgow (July 2023).

For any questions regarding this code, please send a message to: maxsanara@gmail.com.

# Basic requirements
* Quantum Studio 3 datasets of target gene and reference gene (Tab-separated text-quoted CSV).
* Bash based terminal.
* R (including the ```xtable``` package).

# Introduction
The following pipelines have been designed to work with a CSV export of the Excel file format of the results provided by QuantumStudio 3 software for real time polymerase chain reaction (qPCR) analysis. It is mandatory that the source CSVs are tab separated and their text fields quoted.

Two different analysis are included that can be executed independently. In both cases the procedure depends on the existence of a two datasets, one for the a target gene whose expression is analyzed and one for the reference gene (as a housekeeping gene) whose expression does not change presumably over the course of the experiment. Each data set contains meaningful experimental categories, such as different types of controls or mutants under study. Each category must be represented by different samples (i.e. _biological replicas_) which in turn are split into more replicas (i.e. _technical replicas_). Technical replicas try to cope with instrumental errors, e.g. pipetting, fluorescence noise, etc. and could be treated separatedly.

The two aformentioned analysis can be defined as follows: 

* **Through standard curve of dilution factors (1).** This method relies on the estimated dilution factor of the target mRNA through a standard curve of known *dilution factors*. Technical replicas of each biological sample are first assessed, kept or discarded according to their standard deviation[^1], before they are reduced to a central estimation such as their averages. These averages are regarded as representative biological measurements gathered under one category (e.g. one type of control). These averages are later normalised over the partnered average from the reference gene dataset. Two normalizations are provided: dilution factors ratios and differences of C_t (i.e. Delta C_t). Finally every pair of categories is compared by both Kolmogorov-Smirnov (KS) and t-Student tests.
* **DeltaDetal Ct estimation with Montecarlo approximation (2).** In this case no distinction is made between technical replicas and biological samples as long as the number of technical replicas per sample is the same. Measurements associated to one category are put together and each category mean distribution estimated, corresponding to the parametrization of a t-Student distribution. Then the distribution of C_t differentials (i.e. Delta C_t) are determined between the target gene and the reference gene, and the distribution of differentials between Delta C_t (i.e. DeltaDelta C_t) ascertained afterwards. Since these distributions have no analytic expression, a Montecarlo approach is undertaken over a large number of points (e.g. 10^4). Finally the relative level of expression is established by the distribution of 2^{-DeltaDelta C_t}, where 1 means no change, >1 means upregulation and <1 downregulation. KS and t-Student tests for 

The results are given in PDF, SVG (charts) and LaTex compatible codes. Blank charts with no label are also included for TikZ annotation in LaTeX code. Tables results are also rendered for further analysis in R.

# Execution

## Through a standard curve of dilution factors.
Execution can be performed by ```bash make_qpcr_analysis_1c.sh``` including any essential information according to the following flags:
```
bash make_qpcr_analysis.sh --wd ${HOME}
 	--analysis.folder qPCR_analysis
	--data.folder qPCR_data
	--exp.folder default
	--exp.file default.csv
	--plate.ranks a1:b4,h8:k9,...
	--categs standard,a1b4,h8k9,...
	--categs.names "Standard,control, experimental ...
	--alpha 0.05
	--wgraph 0.75
	--wh 2,6
	--delete.outliers TRUE
	--delete.highsd FALSE
```
* ```--wd <string>```: This is the top working directory (absolute path) were global folders for data and results are located. 
* ```--analysis.folder <string>```: Name of the analysis folder where results will be stored. Avoid spaces for compatibility!
* ```--data.folder <string>```: Name of the data folder. Compatible files are in CSV format, with tab as separator and quoted text. Decimal character is a dot. Avoid spaces for compatibility!
* ```--exp.folder <string>```: Name of the experimental folder. It will be created automatically if missing. Avoid spaces for compatibility.
* ```--exp.file <string>```: Name of the files from QS3. Avoid spaces for compatibility.
* ```--plate.ranks <string>```: String containing the ranks of different plate categories separated by commas, in excel format.
* ```--categs <string>```: Categories associated to the plate ranks for internal data structure (defaults are give as excel format ranks without the colon character).
* ```--categs.names <string>```: Actual name of the categories, expressing conditions, genotypes, temperature, etc. meant to be displayed in the graphical output (chart). It can contain spaces.
* ```--alpha <float>```: Signification level for the CI and tests.
* ```--wgraph <float>```: Relative width of each bar plot.
* ```--wh <float,float>```: Width and height of the PDF and SVG file, in inches.
* ```--delete.outliers <boolean>```: Should outliers from QS3 files be removed?
* ```--delete.highsd <boolean>```: Should points with high SD be removed from QS3 files?

Default values to some of these flags are provided when missing through the file ```config_1c```, which currently contains the following structure:
```
ANFLD=qPCR_analysis
DATAFLD=qPCR_data
DATAEXPFLD=default
DATAEXPFILE=default.cvs
DATAREFFILE=default_ref.csv
STDCURVE_POS=1
SKIPLINES=49
ALPHA=0.05
WGRAPH=0.75
WH=2,6
DELOUT=TRUE
DELHIGHSD=FALSE
EXCLUDING_FIELDS=1,3,7,8,10,11,13,14,25,30,31,32,33,29
FACTORING_FIELDS=Target.Name,Task,Sample.Name,Well.Position
```
```STDCURVE```, ```EXCLUDING_FIELDS``` and ```FACTORING_FIELDS``` are fixed parameters to be updated on newer versions of Quantum Studio 3 and should not be altered. ```SKIPLINES``` is currently deprecated. 

## DeltaDelta Ct estimation with Montecarlo approximation. 
Execution can be undertaken on a terminal by ```bash make_qpcr_analysis1.sh``` for a preliminary data extraction and later by ```bash make_qpcr_analysis2.sh``` for the Montecarlo assay and the statistical tests. Their syntaxes go as follows:
```
bash make_qpcr_analysis.sh --wd ${HOME}
	--analysis.folder qPCR_analysis
	--data.folder qPCR_data
	--exp.folder default
	--exp.file default.csv
	--plate.ranks a1:b4,h8:k9,...
	--categs standard,a1b4,h8k9,...
	--categs.names "Standard,control, experimental ..."
	--alpha 0.05
	--wgraph 0.75
	--wh 2,6
	--delete.outliers TRUE
	--delete.highsd FALSE 

```
* ```--wd <string>```: This is the top working directory (absolute path) were global folders for data and results can be found. 
* ```--analysis.folder <string>```: Name of the analysis folder where results will be stored. Avoid spaces for compatibility!
* ```--data.folder <string>```: Name of the data folder. Compatible files are in CSV format, with tab as separator and quoted text. Decimal character is a dot. Avoid spacies for compatibility!
* ```--exp.folder <string>```: Name of the experimental folder. Avoid spaces for compatibility!
* ```--exp.file <string>```: Name of the files from QS3. Avoid spaces for compatibility!
* ```--plate.ranks <string>```: String containing the ranks of different plate categories separated by commas, in excel format.
* ```--categs <string>```: Categories associated to the plate ranks for internal data structure (defaults are give as excel format ranks without the colon character).
* ```--categs.names <string>```: Actual name of the categories, expressing conditions, genotypes, temperature, etc. meant to be displayed in the graphical output (chart).
* ```--alpha <float>```: Signification level for the CI.
* ```--wgraph <float>```: Width of each bar plot.
* ```--wh <float,float>```: Width of the PDF file.
* ```--delete.outliers <boolean>```: Should outliers from QS3 files be removed?
* ```--delete.highsd <boolean>```: Should points with high SD be removed from QS3 files?


Default values are provided by the ```config``` file as follows:

```
ANFLD=qPCR_analysis
DATAFLD=qPCR_data
DATAEXPFLD=default
DATAEXPFILE=default.cvs
STDCURVE_POS=1
SKIPLINES=49
ALPHA=0.05
WGRAPH=0.75
WH=2,6
DELOUT=TRUE
DELHIGHSD=FALSE
EXCLUDING_FIELDS=1,3,7,8,10,11,13,14,25,30,31,32,33,29
FACTORING_FIELDS=Target.Name,Task,Sample.Name,Well.Position
```

Meaning the same as in the previous section.

As for ```make_qpcr_analysis.sh``` the syntaxis is defined as follows:

```
bash make_qpcr_analysis2.sh --wd ${HOME}
     --analysis.folder qPCR_analysis
     --data.folder qPCR_data
     --exp.folder default
     --exp.file default.csv
     --ref.file default_ref.csv
     --exp.levels A1-4,B5-8, ...
     --ctr.levels C1-4,D5-8, ...
     --montecarlo.n 1e+4
     --qranks 0.01,0.05,0.25,...
     --wgraph 0.75
     --wh 2,8
     --lve A1-4xC1-4,B5-8xD5-8,...
     --cve D1-4xE1-4,F5-8xH1-4,...
     --dig 2
```

* ```--wd <string>```: This is the top working directory (absolute path) were global folders for data and results can be found. 
* ```--analysis.folder <string>```: Name of the analysis folder where results will be stored. Avoid spaces for compatibility!
* ```--data.folder <string>```: Name of the data folder. Compatible files are in CSV format, with tab as separator and quoted text. Decimal character is a dot. Avoid spacies for compatibility!
* ```--exp.folder <string>```: Name of the experimental folder where experimental file and reference files are located. Avoid spaces ...
* ```--exp.file <string>```: Name of the experimental file from QS3. Avoid spaces ...
* ```--ref.file <string>```: Name of the reference gene file from QS3. Avoid spaces ...
* ```--exp.levels <string>```: Name of the categories that will be normalised to the control in the same position (check ```--ctr.levels```)
* ```--ctr.levels <string>```: Control names that will be compared to the specified exp.levels.
* ```--montecarlo.n <integer>```: Number of points in the montecarlo simulation.
* ```--qranks <integer,integer>```: Probability levels used to calculate the quantiles.
* ```--wgraph <float>```: Width of the graphical output
* ```--wh <float,float>```: Width and height in inches of the PDF output
* ```--lve <string>```: Composed levels used to make comparisons (tests) with the level at the same position in --cve
* ```--cve <string>```: Composed levels used to make comparisons (tests) with the level at the same position in --lve
* ```--dig <integer>```: Number of significant digits in the latex table

Deafults are found at ```config_ref```:
```
ANFLD=qPCR_analysis
DATAFLD=qPCR_data
DATAEXPFLD=default
DATAEXPFILE=default.cvs
DATAREFFILE=default_ref.csv
MONTECARLO_NPOINTS=1e+4
QRANKS=.01,.05,.25,.5,.75,.95,.99
WGRAP=0.75
WH=2,8
DIG=2
```

# Tasks
- [ ] So far method **(2)** does not include an error reduction due to the paired structure between the target gene and the reference gene datasets.

[^1]: Quantum Studio 3 includes two different informations, i.e. *flags*, on this regard: high _standard deviation_ and _outlier_. 

