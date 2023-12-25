# Internship code

This repository contains code which is used in a summer research project at the University of Warwick, under the supervision of Dr Paul Strøm.

The outline of the project is the following: taking Kepler data, or PLATO simulated data for Sun-like stars, randomly injecting transits of different periods, planet radii and impact parameters to them, then looking at the percentages of recovery. 

The repository contains the following folders:

- "inject_retrieve", which has most of the code. The only file that should be changed by the user (unless debugging is needed) is "current_parameters.py". "ir_main.py" is the one that must be run.

- "permanent", which has a programme to download Kepler curves, and one to bin PLATO lighcurves to any desired cadence (multiple of 25 seconds).

- "inject-retrieve_Thesis.pdf", as its name suggests, is my undergraduate thesis.


__Be aware:__ This code is for testing only.
