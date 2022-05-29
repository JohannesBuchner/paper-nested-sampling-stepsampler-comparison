import scipy.stats
import numpy as np
import matplotlib.pyplot as plt
#from ultranest.mlfriends import ScalingLayer, AffineLayer, MLFriends
#from ultranest.stepsampler import RegionMHSampler, CubeMHSampler
from ultranest.stepsampler import SliceSampler
from ultranest.stepsampler import generate_random_direction, generate_cube_oriented_direction, generate_cube_oriented_differential_direction
from ultranest.stepsampler import generate_differential_direction, generate_partial_differential_direction, generate_region_oriented_direction
from ultranest.stepsampler import generate_region_random_direction, generate_mixture_random_direction
from ultranest.stepsampler import OrthogonalDirectionGenerator, SequentialRegionDirectionGenerator
from ultranest.dychmc import DynamicCHMCSampler
from problems import get_problem
from evaluate_sampling import evaluate_warmed_sampler

def slicesampler_with(generate_direction):
    return lambda nsteps, ndim: SliceSampler(nsteps=nsteps, generate_direction=generate_direction)
def slicesampler_withclass(generate_direction_class):
    return lambda nsteps, ndim: SliceSampler(nsteps=nsteps, generate_direction=generate_direction_class())

def ortho_slicesampler_with(generate_direction):
    return lambda nsteps, ndim: SliceSampler(nsteps=nsteps, generate_direction=OrthogonalDirectionGenerator(generate_direction))

def chmc(adaptive_nsteps=False):
    return lambda nsteps, ndim: DynamicCHMCSampler(0.1 * 0.5**ndim, nsteps=nsteps, adaptive_nsteps=adaptive_nsteps)


samplers = {
    'cube-slice': slicesampler_with(generate_cube_oriented_direction),
    'region-slice': slicesampler_with(generate_region_oriented_direction),
    'region-seq-slice': slicesampler_withclass(SequentialRegionDirectionGenerator),
    'cube-harm': slicesampler_with(generate_random_direction),
    'region-harm': slicesampler_with(generate_region_random_direction),
    'cube-ortho-harm': ortho_slicesampler_with(generate_random_direction),
    'region-ortho-harm': ortho_slicesampler_with(generate_region_random_direction),
    'de-harm': slicesampler_with(generate_differential_direction),
    # this is the same as cube-slice, except with better adapt guess
    # so the needed number of steps will be the same!
    'de1': slicesampler_with(generate_cube_oriented_differential_direction),
    'depart-harm': slicesampler_with(generate_partial_differential_direction),
    'de-mix': slicesampler_with(generate_mixture_random_direction),
    'chmc': chmc(),
}


"""
from ultranest.popstepsampler import PopulationSliceSampler, generate_cube_oriented_direction, \
   generate_random_direction, generate_region_oriented_direction, generate_region_random_direction

def make_popsampler(generate_direction, popsize=40):
    def make_sampler(nsteps):
        return PopulationSliceSampler(
            popsize=popsize, nsteps=nsteps,
            generate_direction=generate_direction)
    return make_sampler

samplers = {
    'cube-slice': make_popsampler(generate_cube_oriented_direction),
    'region-slice': make_popsampler(generate_region_oriented_direction),
    'cube-harm': make_popsampler(generate_random_direction),
    'region-harm': make_popsampler(generate_region_random_direction),
}
"""

ordered_problems = [
    ('shell', 2, 3000),
    ('pyramid', 4, 10000),
    ('shell', 8, 6000),
    ('pyramid', 16, 10000),
    ('corrgauss', 16, 10000),
    ('corrgauss', 100, 10000),
]

# calibrate it
def calibrate(samplername, nlive=400, nshrinkages=25):
    nsteps_current = 1
    nevalsteps_total = nshrinkages * nlive - 2
    
    ordered_problems_remaining = list(ordered_problems)
    
    # select first problem
    problemname, ndim, nevalsteps = ordered_problems_remaining.pop(0)
    fout = open('calibration/%s-config.log' % samplername, 'w')
    
    while True:
        sampler = samplers[samplername](nsteps_current, ndim)
        
        print("evaluating sampler: %s on problem %s-%d" % (sampler, problemname, ndim))

        loglike, grad, volume, warmup = get_problem(problemname, ndim=ndim)

        all_shrinkages = []
        all_ncalls = 0
        seed = 0
        fishy = False
        while len(all_shrinkages) < nevalsteps_total:
            seed += 1
            sampler = samplers[samplername](nsteps_current, ndim)
            Lsequence, ncalls, steps = evaluate_warmed_sampler(problemname, ndim, nlive, nevalsteps, sampler, seed=seed)
            stepsize, angular_step, radial_step = steps.transpose()
            del sampler
            deltaLmin = np.nanmin(np.diff(Lsequence))
            print("    smallest steps in L:", deltaLmin, "euclidean:", stepsize.min(), "angle:", angular_step.min(), "radial:", radial_step.min())
            if deltaLmin == 0 or not angular_step.min() > 0 or not stepsize.min() > 0 or not radial_step.min() > 0:
                print("Got stuck with nsteps=%d, doubling" % nsteps_current)
                fishy = True
                break
            
            assert np.isfinite(Lsequence).all(), Lsequence
            vol = np.asarray([volume(Li, ndim) for Li in Lsequence])
            assert np.isfinite(vol).any(), ("Sampler has not reached interesting likelihoods", vol, Lsequence)
            shrinkage = 1 - (vol[np.isfinite(vol)][1:] / vol[np.isfinite(vol)][:-1])**(1. / ndim)
            
            all_shrinkages = np.concatenate((shrinkage, all_shrinkages))
            all_ncalls += ncalls
            print("   have %d/%d shrinkages" % (len(all_shrinkages), nevalsteps_total))

        if fishy:
            nsteps_current *= 2
            continue

        efficiency = len(all_shrinkages) / all_ncalls
        all_shrinkages = all_shrinkages[:nevalsteps_total]
        
        # convert to a uniformly distributed variable, according to expectations
        cdf_expected = 1 - (1 - shrinkage)**(ndim * nlive)
        
        result = scipy.stats.kstest(cdf_expected, 'uniform')
        print("    ", result)
        prefix = 'calibration/%s-%s-%d' % (samplername, problemname, ndim)
        np.savetxt('%s-shrinkage.txt.gz' % prefix, shrinkage)
        np.savetxt('%s-shrinkage-cdf.txt.gz' % prefix, cdf_expected)
        
        if result.pvalue < 0.01:
            print("Deviation detected in sampler with nsteps=%d, doubling" % nsteps_current)
            nsteps_current *= 2
        else:
            fout.write('%d\t%d\t%.10f\t%s\n' % (nsteps_current, ndim, efficiency, problemname))
            fout.flush()
            print("Sampler looks ok")
            print("Calibration with %s-%d: %d steps needed" % (problemname, ndim, nsteps_current))
            if ordered_problems_remaining:
                print()
                problemname, ndim, nevalsteps = ordered_problems_remaining.pop(0)
                print("Next problem: %s-%d" % (problemname, ndim))
            else:
                break
    fout.close()
    nsteps, ndims = np.loadtxt('calibration/%s-config.log' % samplername, usecols=(0,1), unpack=True)
    # i = np.argmin(nsteps)  # it is always the first, because monotonous
    i = 0
    k = int(np.ceil(np.nanmax((nsteps - nsteps[i]) / (ndims - ndims[i]))))
    N0 = max(0, nsteps[i] - k * ndims[i])
    k = int(np.ceil(np.nanmax((nsteps - N0) / ndims)))
    plt.plot(ndims, nsteps, drawstyle='steps-post', marker='*', label=samplername)
    plt.plot([1, ndims.max()], [N0 + k, ndims.max() * k + N0], ':', label=r'%d + %d $\times$ d' % (N0, k))
    plt.xlabel('Dimensionality')
    plt.ylabel('Number of steps')
    plt.legend(loc='upper left')
    plt.savefig('calibration/%s-config.pdf' % samplername)
    plt.close()

def main():
    import sys
    # pick a sampler
    return calibrate(sys.argv[1])

if __name__ == '__main__':
    main()
