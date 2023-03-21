Code used to generate results, figures and tables for:

Published article: Diversification and change in the R programming language, Royal Society Open Science (URL pending)
Pre-print [Expansion and evolution of the R programming language](https://arxiv.org/abs/2208.12382).

The published version is revised and best version to refer to. The latest update to this code is for this version.

The data used to generate this publication were obtained using the GitHub API.

galaxy\_processing.R uses the API to access and process publically accessible data. 

galaxy\_analysis.R imports the cleaned data and runs all analyses and produces all figures. 

Both scripts makes extensive use of the code folding functionality in RStudio. Alt + O is the default shortcut on Windows and Linux versions of RStudio to collapse all folds.

**galaxy\_processing.R**

The gitHub aquisition process in this script requires several discrete steps. Due to API request limits, some steps can only be performed at a maximum of 30 per minute, regardless of computer processing power. If you attempt to re-create the entire dataset from the publication, it will take a **LONG** time. Like weeks long.  Maybe months. Seriously. It's a huge investment! I'm pretty sure you don't care *that much* about my work. If all you want to do is understand how the process works, I strongly suggest you set the number of calendar days in step #3 to something manageable (e.g., 100 days).

This script has several requirements in order to run correctly.

1. You must set up a [personal access token](https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/creating-a-personal-access-token) in GitHub to use in the script. This token is used to query the GitHub API. This needs to be added at line 19.

2. You must download this entire repository, including all subfolders, and run galaxy\_processing.R from within this structure.

3. The number of calendar days queried by the script is set is the lapply() on line 27 (in the *Identify repositories with the "R" language flag* subsection). I set this to **1:4382**, to capture all data from Jan 1 2010 to Dec 31 2021 (though subsequently only 2014-2021 were used for analyses).

4. The latter part of the script associates R functions with packages. It will only work on packages that are locally installed on your machine. There are several blocks of code wrapped in if statements set to "FALSE". These will uninstall all packages in your R installation, and then install sets of packages. *I recommend NOT running these blocks*. They are included for completeness only.

5. The script executes the following process (with accompanying folders). None of the files from these intermediate steps are included in the repository, given size and number of files.
	
	1. Identifies all GitHub repositories created on a given calendar day that contain the "R" language flag. These are stored in the *gitCalls* folder.
	2. Each of these repositories is then queried individually to access a tree structure containing all files. The path to all ".R" files is accessed, which are stored in the *gitURLs* folder. Repositories with the R flag but no .R files are stored in *gitBlanks*. On the off chance that a repository has been deleted between steps 1 and 2, these missing repositories are stored in the *gitGone* folder.
	3. The script then downloads each .R file in its entirety, storing them in *gitScripts*. 
	4. Each script is run through an algorithm to extract function calls. Details are in the manuscript, but I search for text strings to the left of each open parenthesis ("("). A vector of functions for each script is then stored in *gitFuns*.
	5. Function calls are then converted into long-form count data, including repository level metadata, which are output as a large R list in */outputs/funTableList.rds*.
	6. Binding this list together is conducted in blocks of 10,000 scripts, each of which are stored in *gitFunTables*.
	7. The raw processed data are then stored in */outputs/gitFunctionLong.csv*.
	8. After this, functions are associated with R packages, and processed as per the manuscript, including applying a sampling filter to remove extremely rare functions.
	9. The final processed data is stored as *./outputs/commonFunctionLong.csv*.

**galaxy\_analysis.R**

This script imports the processed data from galaxy\_processing.R and conducts all statistical analyses in the main text and supplementary material. 

Given the time investment to recreate this file, I have housed a full version (including metadata) for download as a [Dryad dataset](https://doi.org/10.5061/dryad.h18931zrg). It is housed outside of the repository due to size (c. 7GB). To run galaxy\_analysis.R, this file must be downloaded and placed into the *outputs* subfolder of this repository.

galaxy\_analysis.R also relies on a number of functions in the *functions* subfolder, so please make sure these are included.

Some of the analysis steps take a long time to run, as this is a very large dataset. If you are interested in quickly reproducing intermediate steps, I suggest subsetting *funTableSub* when it is created on line 23. 100,000 rows should be sufficient to capture enough data to run analyses.

All plots are generated within the R script, with only labels being altered in post-R software. These are saved in the *plots* sub-folder.

Please let me know if you encounter any issues with these scripts or associated data.