# SpatialSampling
Code to support Mastin et al. (2019) Optimising risk-based surveillance for early detection of invasive plant pathogens. https://doi.org/10.1101/834705

There are two steps:
	1) Generating a set of runs from the landscape scale simulation
	2) Using the results to generate a sampling pattern
[note, the performance of this pattern must of course be tested on an independent set of simulation runs!]

The program landscapeScaleSimulation does the first step; the program simulatedAnnealing does the second.

The Makefiles show how each can be compiled.

The .bat (Windows) or .sh (Linux) scripts shows how to run them.

Each program has a corresponding .cfg file, which controls configurable options (e.g. epidemic model parameterisation).

The .R script (which has only been tested on Windows) will process the outputs and georeference it, such that coordinates of individal sites can be identified. It will output plots of maps of the optimised distribution of sampling sites over a map of the mean end prevalence in each cell, the probability of infection of each cell, and another gridded dataset (currently selected as the citrus density).

The spatial maps of the host landscape and of primary infection rates are not freely sharable, so here we use different date (with corresponding changes to the parameterisation of the epidemic model).

The gridded dataset of citrus density ("propFull.txt") is taken from the county-level numbers of commerical citrus trees in Florida, taken from the Florida Citrus Statistics 2018-2019. 5,000 additional trees were added to each county to represent residential trees. Tree numbers were then adjusted into densities according to the total area of the county, and then scaled such that the total density of citrus was equal to a high-resolution dataset.

The gridded dataset of relative primary infection rates ("relPrimaryInf.txt") is taken from predictions from Gottwald et al 2019 ("A probabilistic census travel model to predict introduction sites of exotic plant, animal and human pathogens"; Philos Trans R Soc B Biol Sci;374: 20180260. doi:10.1098/rstb.2018.0260), based on predictions for the year 2010.
