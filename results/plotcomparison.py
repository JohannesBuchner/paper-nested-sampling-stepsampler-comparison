import sys
import matplotlib.pyplot as plt
import numpy as np

samplers = [
    'cube-slice',
    'cube-harm',
    'cube-ortho-harm',
    'region-slice',
    'region-seq-slice',
    'region-harm',
    'region-ortho-harm',
    'de1', 
    'de-harm',
    #'depart-harm',
    'de-mix',
]

cax = plt.figure('calib').gca()
eax = plt.figure('efficiency').gca()
x = np.arange(8, 30)

cax.plot(x, x / 3, ls='--', color='k')
cax.text(16, 16 / 3, 'linear', color='k', rotation=22, va='center', ha='right')

eax.plot(x, 4e-1 / x, ls='--', color='k') #, label='$1 / d$ scaling')
eax.text(8, 4e-1 / 8, '$1 / d$ scaling', color='k', rotation=-17, va='center', ha='left')

latextable = open('samplertable.tex', 'w')

latextable.write(r"""
\begin{tabular*}{1\textwidth}{@{\extracolsep{\fill}}>{\centering}p{0.33\textwidth}>{\centering}p{0.33\textwidth}>{\centering}p{0.33\textwidth}}
\toprule 
\textbf{Sampler} & \textbf{Factor $k$} & \textbf{Efficiency \%}\tabularnewline
\midrule 
""")

for samplername in samplers:
	point_kwargs = dict(
		marker = 'x' if 'region' in samplername else 'o' if 'de' in samplername else 's',
		ls=' ',
	)
	line_kwargs = dict(
		ls = '--' if 'region' in samplername else ':' if 'de' in samplername else '-',
		marker = point_kwargs['marker'],
		mfc='none', mew=2, ms=8,
		#mfc = (0, 0, 0, 0.0),
		#drawstyle='steps-post', 
	)
	try:
		nsteps, ndims, efficiencies, problemnames = np.loadtxt('%s-config.log' % samplername, unpack=True,
			#converters={0:int,1:int,2:float,3:str}
			dtype=dict(
				names='nsteps ndims efficiency problemname'.split(),
                formats=('i4', 'i4', 'f4', 'S10'))
			)
		#print(sampler, ndims, nsteps)
		if len(ndims) == 0: continue
		print("loaded %s" % samplername)
	except (Exception, IOError) as e:
		print("No data for %s" % samplername, e)
		continue
	bad_calibration = ndims.max() < 100
	if bad_calibration:
		ndims = np.concatenate((ndims, [100]))
		nsteps = np.concatenate((nsteps, [2048]))
		efficiencies = np.concatenate((efficiencies, [1e-10]))
	cax.plot(ndims, nsteps, label=samplername, **line_kwargs)
	eax.plot(ndims, efficiencies, label=samplername, **line_kwargs)
	k = int(np.ceil(np.nanmax(nsteps * 1. / ndims)))
	if bad_calibration:
		latextable.write("%s & >16 & --  \\tabularnewline\n" % (samplername)) # , 100 * min(efficiencies[:-1] * ndims[:-1])))
	else:
		latextable.write("%s & %d & %.2f \\tabularnewline\n" % (samplername, k, 100 * min(efficiencies * ndims)))

	if samplername == 'cube-slice':
		ax1 = plt.figure().gca()
		lastd = 0
		for d, nstep, problemname in zip(ndims, nsteps, problemnames):
			if (ndims == d).sum() == 1:
				ha = 'center'
				y = nstep * 1.2
			if (ndims == d).sum() == 2:
				if lastd != d:
					y = nstep * 1.2
					ha = 'right' 
				else:
					y = nstep * 1.6
					ha = 'left'
			
			#cax.text(d, 1000, problemname.decode('utf-8'), va='top', ha=ha, rotation=90)
			#if d == 16:
			#	ax1.text(d, 1000, problemname.decode('utf-8'), va='top', ha=ha, rotation=90)
			#else:
			ax1.text(d, y, problemname.decode('utf-8'), va='bottom', ha='center')
			cax.text(d, y, problemname.decode('utf-8'), va='bottom', ha='center')
			eax.text(d, 0.4, problemname.decode('utf-8'), va='top', ha=ha, rotation=90)
			lastd = d
		
		plt.plot([2, ndims.max()], [2 * k, ndims.max() * k], ':', label=r'%d $\times$ d scaling' % (k), color='k')
		plt.plot(ndims, nsteps, drawstyle='steps-post', **line_kwargs, label=samplername)
		plt.fill_between(ndims, nsteps * 0, nsteps, hatch='X', facecolor='none', edgecolor='lightgray', step='post', alpha=0.5)
		plt.text(10, 8, 'Rejected', size=14, bbox=dict(facecolor='white', edgecolor='none'))
		plt.text(3, 400, 'Allowed', size=14)
		plt.xscale('log')
		plt.yscale('log')
		plt.ylim(1, 1024)
		plt.xticks([2, 4, 8, 16, 32, 100], [2, 4, 8, 16, 32, 100])
		#plt.xlim(2, 100)
		plt.xlabel('Dimensionality')
		plt.ylabel('Calibrated number of steps')
		plt.legend(loc='upper right')
		plt.savefig('plotcomparison1.pdf')
		plt.close()

latextable.write(r"""
\bottomrule
\end{tabular*}
""")
latextable.close()

plt.figure('calib')
plt.ylabel('Calibrated number of steps')
plt.xlabel('Dimension')
plt.xscale('log')
plt.yscale('log')
plt.xticks([2, 4, 8, 16, 32, 100], [2, 4, 8, 16, 32, 100])
plt.ylim(1, 1024)
plt.legend(loc='lower right', numpoints=1)
plt.savefig('plotcomparison.pdf', bbox_inches='tight')
plt.close()

plt.figure('efficiency')
plt.ylabel('Efficiency (model eval. / iteration)')
plt.xlabel('Dimension')
plt.xscale('log')
plt.yscale('log')
plt.xticks([2, 4, 8, 16, 32, 100], [2, 4, 8, 16, 32, 100])
plt.ylim(3e-5, 0.42)
plt.legend(loc='upper right', numpoints=1)
plt.savefig('plotcomparison_efficiency.pdf', bbox_inches='tight')
plt.close()
