::
:: Run the simulations. 
:: 
:: Note that parameters are specified in the configuration file landscapeScaleSimulation.cfg
:: although this can be overriding by passing a <key>=<value> pair on the command line
:: (illustrated below for maxIncidence)
::
landscapeScaleSimulation maxIncidence=0.02

::
:: Produce an optimised sampling pattern based on these simulation results
::
simulatedAnnealing
