#!/usr/bin/env python3
import numpy as np

try:
    from pandas import read_csv

    pandas_avail = True
except ImportError:
    pandas_avail = False
    print("No pandas found, will fall back to the slower numpy-based version")
from scipy.interpolate import interp1d

try:
    import matplotlib.pyplot as plt

    mpl_enabled = True
except ImportError:
    mpl_enabled = False
from optparse import OptionParser
import sys


class Wham:
    def __init__(self, opt):
        """
        Implements the self-consistent WHAM algorithm to calculate
        the free energy profile from Umbrella Sampling simulations
        :param opt:
        """
        # options
        self.opt = opt
        # general settings
        self.nbins = None
        self.p_grid = None
        self.bin_edges = None
        self.nwindows = None
        # actual data and results
        self.pot_interp = None
        self.data = None
        self.combined = None
        self.weights = None
        self.num = None
        self.p = None
        self.g = None
        self.f = None
        # error analysis
        self.autocorr = None
        self.boot_framenums = None
        self.gs_boot = None
        self.gs_error = None
        self.precision = np.sqrt(np.finfo(float).eps)

    def __call__(self, meta, multi=False):
        """
        the main routine for running WHAM with a single metafile,
        wrapped together to enable usage from outside the script
        :param meta: str, metafile to read from
        :param multi: bool, whether multiple calculations were requested
        :return: self.p_grid: grid values for the collective variable
                 self.g: free energy profile
                 self.f: F values
        """
        self.read_data_and_setup_grids(meta)
        self.f, self.g = self.iterate_to_consistency(self.num, self.weights)
        if self.autocorr:
            self.bootstrap()
        self.log_results(multi, meta)
        return self.p_grid, self.g, self.f

    def read_data_and_setup_grids(self, meta):
        """
        reads the metafile to determine all initial parameters,
        including colvar bounds, grids, bootstrap flags; also
        reads data into memory
        :param meta: str, metafile to read from
        :return: None
        """
        files = [x.strip().split() for x in open(meta)]
        self.check_if_bootstrap(files)
        self.nwindows = len(files)
        if not all([len(q) == len(files[0]) for q in files]):
            raise ValueError('not all lines in metafile have the same num of cols, likely an error')
        if pandas_avail:
            self.data = [read_csv(q[0], header=None, delim_whitespace=True).as_matrix()[:, 1] for q in files]
        else:
            self.data = [np.loadtxt(q[0])[:, 1] for q in files]
        self.combined = np.concatenate(self.data)
        # determine min and max CV vals to specify bounds automatically
        vmin, vmax = np.min(self.combined), np.max(self.combined)
        self.nbins = self.opt.n
        self.weights = self.get_weights(files)
        potentials = self.get_pots(files, vmin, vmax, self.nbins)
        # setup grids, interpolate potentials and calc histograms
        self.p_grid = np.linspace(vmin, vmax, self.nbins)
        max_v = np.max(np.concatenate([x[:, 1] for x in potentials]))
        self.pot_interp = [interp1d(q[:, 0], q[:, 1], bounds_error=False, fill_value=max_v)(self.p_grid)
                           for q in potentials]
        spacing = self.p_grid[1] - self.p_grid[0]
        self.bin_edges = np.linspace(vmin - spacing / 2, vmax + spacing / 2, self.nbins + 1)
        self.num = np.histogram(self.combined, self.bin_edges, weights=np.concatenate(self.weights))[0]

    def get_weights(self, files):
        """
        if you want custom weights for individual frames, just add 3rd col
        to your data files; will be recognized and read by this fn
        :param files: list, paths to data files, as included in the metafile
        :return: list of arrays that contain weights corresponding to individual frames
        """
        contains_weights = np.loadtxt(files[0][0]).shape[1] == 3
        if contains_weights:
            if pandas_avail:
                weights = [read_csv(q[0], header=None, delim_whitespace=True).as_matrix()[:, 2] for q in files]
            else:
                weights = [np.loadtxt(q[0])[:, 2] for q in files]
            weights_cat = np.concatenate(weights)
            weights = [w * weights_cat.shape[0] / np.sum(weights_cat) for w in weights]
        else:
            weights = [np.ones(len(window), dtype=int) for window in self.data]
        return weights

    @staticmethod
    def get_pots(files, vmin, vmax, nbins):
        """
        for compatibility with Grossfield's version - 3 cols in the metafile
        imply that data are read as data_file, pot_center, frc_constant
        :param files: list, paths to files with tabulated potentials
        :param vmin: float, beg of interval to consider
        :param vmax: float, end of interval to consider
        :param nbins: int, number of bins
        :return: list, contains np.arrays with tabulated potentials
        """
        try:
            _ = float(files[0][1])
            pots_are_harmonic = True
        except ValueError:
            pots_are_harmonic = False
        if pots_are_harmonic:
            centers = [float(q[1]) for q in files]
            ks = [float(q[2]) for q in files]
            windows = len(files)
            potentials = [np.vstack([np.linspace(vmin, vmax, nbins),
                                     0.5 * ks[w] * (np.linspace(vmin, vmax, nbins) - centers[w]) ** 2]).T
                          for w in range(windows)]
        else:
            potentials = [np.loadtxt(q[1]) for q in files]
        return potentials

    def check_if_bootstrap(self, files):
        """
        reads metafile, returns False if no bootstrapping required
        or a list of autocorr times (in dataframe units) otherwise
        :return:
        """
        if len(files[0]) == 4:
            acorr = [float(line[3]) for line in files]
        elif len(files[0]) == 3:
            try:
                _ = float(files[0][1])
                acorr = False
            except ValueError:
                acorr = [float(line[2]) for line in files]
        else:
            acorr = False
        self.autocorr = acorr

    def get_boot_framenums(self):
        """
        calculates the number of frames to draw from each window while
        bootstrapping, based on autocorrelation time passed in the metafile
        :return: None
        """
        self.boot_framenums = [int(len(self.data[x]) / self.autocorr[x]) for x in range(len(self.data))]

    def iterate_to_consistency(self, num, weights):
        """
        iteratively recalculates G and F until self-consistency
        is reached within given tolerance
        :param num: float, numerator in the WHAM eqn
        :param weights: list, data on which to operate, list of np.arrays
        :return: f: np.array holding the F values for each US window
                 g: np.array with free energy values along the RC
        """
        self.p = np.ones(self.nbins) / self.nbins
        f = np.zeros(self.nwindows)
        kb_t = 0.001986 * self.opt.t  # kcal/mol
        tol = self.opt.tol
        diff = np.infty
        while diff > tol:
            den = np.zeros(len(self.p_grid))
            for j in range(len(self.pot_interp)):
                den += np.sum(weights[j]) * np.exp((f[j] - self.pot_interp[j]) / kb_t)
                # den += len(data[j]) * np.exp((f[j] - self.pot_interp[j]) / kb_t)
            temp = num / den
            diff = np.max(np.abs(temp - self.p))
            self.p = temp
            for j in range(len(f)):
                f[j] = -kb_t * np.log(np.sum(self.p * np.exp(-self.pot_interp[j] / kb_t)))
        np.place(self.p, self.p <= 0, self.precision)
        g = -kb_t * np.log(self.p) - np.min(-kb_t * np.log(self.p))
        return f, g

    def bootstrap(self):
        """
        subsamples the original data multiple times based on
        autocorrelation time, then redoes the iterative calculation
        and reports the distribution of results
        :return: None
        """
        n_tries = self.opt.boot_iter
        self.get_boot_framenums()
        gs_boot_temp = []
        for i in range(n_tries):
            if i % 10 == 0:
                print("Bootstrap iteration {}\n".format(i))
            selection = [np.random.choice(range(len(self.data[x])), self.boot_framenums[x], replace=True)
                         for x in range(len(self.data))]
            sel_frames = [self.data[x][selection[x]] for x in range(len(self.data))]
            sel_weights = [self.weights[x][selection[x]] for x in range(len(self.weights))]
            num = np.histogram(np.concatenate(sel_frames), self.bin_edges, weights=np.concatenate(sel_weights))[0]
            _, g = self.iterate_to_consistency(num, sel_weights)
            gs_boot_temp.append(g)
        self.gs_boot = np.vstack(gs_boot_temp)
        self.calc_error()

    def calc_error(self):
        """
        calculates the free energy error as the standard deviation
        of free energy profiles obtained from bootstrapping
        :return: None
        """
        mean = np.mean(self.gs_boot, axis=0)
        self.gs_error = np.std([g - mean for g in self.gs_boot], axis=0)

    def log_results(self, multi, meta):
        """
        logs the results to 3 files (F values, free energy profile
        and free energy errors, the latter only if desired)
        :param multi: bool, whether multiple calculations were requested
        :param meta: str, metafile name to distinguish between multiple calculations
        :return: None
        """
        if multi:
            np.savetxt(meta + '_' + self.opt.outf, self.f - self.f[0], fmt="%8.3f")
            np.savetxt(meta + '_' + self.opt.out, np.vstack([self.p_grid, self.g]).T, fmt="%8.3f")
            if self.autocorr:
                np.savetxt(meta + '_' + self.opt.outerr, self.gs_error, fmt="%8.3f")
        else:
            np.savetxt(self.opt.outf, self.f - self.f[0], fmt="%8.3f")
            np.savetxt(self.opt.out, np.vstack([self.p_grid, self.g]).T, fmt="%8.3f")
            if self.autocorr:
                np.savetxt(self.opt.outerr, self.gs_error, fmt="%8.3f")

    def do_graphics(self, multi, metas, m, color):
        """
        does the plotting part, plots profiles and error tunnels
        :param multi: bool, whether multiple calculations were requested
        :param metas: list, contains metafile names to be used as titles
        :param m: int, number of axis to use for plotting
        :param color: str, name of the color to be used for plotting
        :return: None
        """
        if self.opt.plotme or self.opt.graph:
            if multi:
                axis = axarr[m]
                title = metas[m]
            else:
                axis = axarr
                title = metas[0]
            axis.plot(self.p_grid[2:-2], self.g[2:-2], color=color)
            axis.set_title(title)
            if self.autocorr:
                axis.fill_between(self.p_grid[2:-2], self.g[2:-2] + self.gs_error[2:-2],
                                  self.g[2:-2] - self.gs_error[2:-2], color=color, alpha=0.3)

    def do_histogram(self, multi, metas, m):
        """
        does the plotting part, plots histograms
        :param multi: bool, whether multiple calculations were requested
        :param metas: list, contains metafile names to be used as titles
        :param m: int, number of axis to use for plotting
        :return: None
        """
        if self.opt.plothist:
            if multi:
                axis = axarr2[m]
                title = metas[m]
            else:
                axis = axarr2
                title = metas[0]
            axis.hist(self.combined, bins=self.p_grid, histtype='stepfilled', color='#dddddd')
            axis.hist(self.combined, bins=self.p_grid, histtype='step', color='#000000')
            axis.set_title(title)
            for ggg in range(self.nwindows):
                if (ggg + 1) % 2 == 0:
                    colcol='#0000ee'
                else:
                    colcol='#ee0000'
                axis.hist(self.data[ggg], bins=self.p_grid, histtype='step', color=colcol)                             

    def show_graphics(self):
        """
        either shows the graphics (with -p) or saves
        to a graphics file (with -g)
        :return: None
        """
        if self.opt.plotme:
            plt.tight_layout()
            plt.show()
        if self.opt.graph:
            plt.tight_layout()
            plt.savefig(self.opt.graph)


def parse_options():
    """
    option parser to run script from terminal
    :return: opt: an options object with run parameters
             metas: list, contains metafile names to be used as titles
             multi: bool, whether multiple calculations were requested
    """
    parser = OptionParser(usage="%prog -m metafile [-o FEP_outfile -f F_outfile -e stderr_outfile"
                                "-T temperature -t tolerance -n num_bins -b bootstrap_iters -g graph_outfile -p]")
    parser.add_option("-m", dest="meta", action="store", type="string",
                      help="metafile: col 1 lists time:CV files (or time:CV:weight) from individual US windows, "
                           "col 2 lists corresponding files with tabulated potentials; "
                           "alternatively col 2 contains the center of the harmonic restraining potential "
                           "and col 3 the force constant, consistently with Alan Grossfield's implementation")
    parser.add_option("-o", dest="out", action="store", type="string", default="free.out",
                      help="outfile for the free energy profile; default is free.out")
    parser.add_option("-f", dest="outf", action="store", type="string", default="F.out",
                      help="outfile for the F values; default is F.out")
    parser.add_option("-e", dest="outerr", action="store", type="string", default="err.out",
                      help="outfile for the standard error values; default is err.out")
    parser.add_option("-T", dest="t", action="store", type="float", default=300.0,
                      help="simulation temperature; default is 300")
    parser.add_option("-n", dest="n", action="store", type="int", default=100,
                      help="number of bins; default is 100")
    parser.add_option("-t", dest="tol", action="store", type="float", default=0.00001,
                      help="tolerance (convergence criterion); default is .00001")
    parser.add_option("-p", dest="plotme", action="store_true", default=False,
                      help="use as flag to plot the resulting FEP")
    parser.add_option("--hist", dest="plothist", action="store_true", default=False,
                      help="use as flag to plot the histogram")
    parser.add_option("-g", dest="graph", action="store", type="string",
                      help="optional filename to save the plotted FEP")
    parser.add_option("-b", dest="boot_iter", action="store", type="int", default=50,
                      help="number of bootstrap iterations to run; default is 50")
    (opt, args) = parser.parse_args()
    if opt.meta is None:
        print("########\n\nmeta file required; use -h to see all options\n\n########")
        sys.exit(1)
    if ' ' in opt.meta:
        metas = opt.meta.split()
        multi = True
    else:
        metas = [opt.meta]
        multi = False
    return opt, metas, multi


if __name__ == "__main__":
    opts, metafiles, multiple = parse_options()
    if opts.plotme or opts.graph:
        if not mpl_enabled:
            raise ImportError('matplotlib could not be loaded, consider installing it via'
                              'pip install matplotlib if visualization is required')
        fig, axarr = plt.subplots(1, len(metafiles), sharey=True)
        if opts.plotme and opts.graph:
            raise ValueError('Simultaneous showing and saving of the free energy profiles currently not supported, '
                             'choose either one')
    if opts.plothist:
        fig2, axarr2 = plt.subplots(1, len(metafiles))
    wh = Wham(opts)
    for mf in range(len(metafiles)):
        pot_grid, gg, _ = wh(metafiles[mf], multiple)
        col = 'C' + str(mf % 10)
        wh.do_graphics(multiple, metafiles, mf, col)
        wh.do_histogram(multiple, metafiles, mf)
    wh.show_graphics()
