Work in progress
----------------

This repository accompanies the paper "Comparison of Step Samplers for Nested Sampling"
(published at https://www.mdpi.com/2673-9984/5/1/46).

Code for running and evaluating the simulations is provided here 
in the spirit of reproducible science.

Usage
-------
::

	$ python3 calibrator.py cube-slice

where cube-slice can be replaced by another method::

	{cube,region}-{slice,harm,ortho-harm} region-seq-slice de-harm de1 depart-harm

This produces `results/cube-slice-config.log`

which contains the calibrated `Nsteps` and cost for each considered problem.

Once that is done, the results are analysed::

	$ cd results/
	$ python3 plotcomparison.py


Code structure
---------------

`calibrator.py` iterates through the problems and applies the shrinkage test to vet configurations. 

It builds the LRPS method from the first input argument.
The LRPS methods are implemented within UltraNest here: 
`API documentation <https://johannesbuchner.github.io/UltraNest/ultranest.html#module-ultranest.stepsampler>`_,
`code <https://johannesbuchner.github.io/UltraNest/_modules/ultranest/stepsampler.html>`_.

Problem functions are taken from `problems.py`, which implements likelihood and volume functions.

Then, `calibrator.py` calls `evaluate_warmed_sampler` implemented in 
`evaluate_sampling.py`, which simulates a nested sampling run for a few thousand iteration, records shrinkage samples and other diagnostics.
