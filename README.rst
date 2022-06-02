Work in progress


Usage
-------
::

	$ python3 calibrator.py cube-slice

where cube-slice can be replaced by another method::

	{cube,region}-{slice,harm,ortho-harm} region-seq-slice de-harm de1 depart-harm

This produces results/cube-slice-config.log

which contains the calibrated Nsteps and cost for each considered problem.

Code structure
---------------

`problems.py` implements the problem likelihood and volume functions.
This is used by `evaluate_sampling.py` sets up UltraNest to run for a few regions, records shrinkage samples and other diagnostics

`calibrator.py` iterates through the problems and applies the shrinkage test to vet configurations. It also builds the LRPS method from the first input argument.
