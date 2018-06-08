In this study we examined whether electrical brain stimulation of the frontal eye field can influence reaction time and accuracy of human eye movements.

* __Project title__: No evidence that frontal eye field tDCS affects latency or accuracy of prosaccades.
* __Project code__: sacc-tDCS
* __Authors__: Reteig, L.C., Knapen T., Roelofs, F.J.F.W., Ridderinkhof, K.R., & Slagter, H.A.
* __Affiliation__: Department of Psychology, University of Amsterdam
* __Year__: 2018

# Resources

Please see this project's [page on the Open Science Framework](https://osf.io/8jpv9/). All the other files from this project can be obtained from there, along with much more (and more detailed) information then in this README.

__Print__: bioRxiv - [bioRxiv / self-archive link]

__Data__: figshare - [figshare link]

# Contents of this repository

* `sacc-tDCS.Rmd`: R notebook with analysis code that will reproduce all of the results, figures and statistics contained in the paper. See the `notebooks` folder for additional notebooks with a more in-depth exploration of the dataset.
* `sacc-tDCS.md`: Markdown output of the notebook, optimized for rendering on GitHub. Plots referenced in the `.md` file are contained in `sacc-tDCS_files/`.
* `analysis.m`: MATLAB script that loads the raw data and creates the `.csv` file that's read in by the R notebook.

* `notebooks/`: Additional notebooks with more in-depth analyses on each different aspect of the project. More information can be found in the notebooks themselves; they are also referred to in the main `sacc-tDCS` notebook. Like the main notebook, each has three types of files (`.Rmd`, `nb.html` and `.md`) and a folder with plots.
* `src/`: All other analysis code used in the project.
* `task/`: Code for the experimental task and control of the eye tracker; used during data collection.

N.B. Each folder has its own `README.md` file with more detailed information on its contents.
