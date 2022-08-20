
from pathlib import Path

import numpy as np
import h5py
from astropy.io import ascii
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import Normalize
import cmasher as cmr
from scipy.stats import binned_statistic_2d

from flags.pz import eazy
import flare.plt as fplt
# from flare.filters import add_filters




# --- open h5py catalogue


parameter_labels = {}
parameter_labels['z'] = 'z'
parameter_labels['log10tau_V'] = r'\log_{10}(\tau_V)'
parameter_labels['fesc'] = 'f_{esc}'
parameter_labels['log10duration'] = '\log_{10}(t_{SF}/yr)'
parameter_labels['log10Z'] = '\log_{10}(Z)'

parameter_range = {}
parameter_range['z'] = [5.01, 14.99]
parameter_range['log10duration'] = [6.01, 8.99]
parameter_range['log10Z'] = [-4.99, -1.39]
parameter_range['fesc'] = [0.01, 0.99]
parameter_range['log10tau_V'] = [-1.99, 0.99]


def simple_plt():

    left  = 0.15
    bottom = 0.15
    width = 0.8
    height = 0.8

    fig = plt.figure(figsize = (3.5, 3.5))
    ax = fig.add_axes((left, bottom, width, height))

    return fig, ax






class Analyser:


    def __init__(self, model, template_set = 'tweak_fsps_QSF_12_v3'):


        self.cat = h5py.File(f'{model}.hf', 'r')
        self.ids = self.cat[f'pz/eazy/{template_set}/id'][:]
        n = len(self.ids)

        self.z = self.cat['parameters/z'][()]
        self.z_pz = self.cat[f'pz/eazy/{template_set}/z_a'][()]
        self.dz = (self.z - self.z_pz)/(1+self.z)

        self.filters = [filter.decode("utf-8")  for filter in self.cat['obs/filters'][:]]
        self.pivwv = {filter: wv for filter, wv in zip(self.filters, self.cat['obs/pivwv'][:])}

        self.template_set =template_set
        self.templates = eazy.get_templates(f'{template_set}.spectra.param')
        self.template_lam, self.template_fnu = eazy.get_template_grid(self.templates)

        Path(model).mkdir(parents=True, exist_ok=True) # create output directory


    def explore(self):
        self.cat.visit(print) # explore datasets

    def print_info(self, i):

        def visitor_func(name, node):
            if isinstance(node, h5py.Dataset):
                 print(name, node[i])

        self.icat.visititems(visitor_func)

    def dz_hist_plot(self):

        fig, ax = simple_plt()



        bin_edges = np.arange(-1.5, 1.5, 0.02)

        ax.hist(self.dz, bins = bin_edges)

        ax.set_xlabel(r'$\rm (z-z_{pz})/z$')
        ax.set_ylabel(r'$\rm N$')

        fn = f'{model}/dz.pdf'
        print(fn)
        fig.savefig(fn)


    def beta_dz_plot(self):

        fig, ax = simple_plt()

        ax.scatter(self.cat['diagnostics/beta'][()], self.dz, color = 'k', s=5, lw=0, alpha = 0.5)

        ax.set_xlabel(r'$\rm \beta$')
        ax.set_ylabel(r'$\rm (z-z_{pz})/(1+z)$')

        ax.set_xlim([-3.5, 2.0])

        fn = f'{model}/beta_dz.pdf'
        print(fn)
        fig.savefig(fn)

    def corner(self, parameters, p, scatter = False, bins = 10):

        n = len(parameters)


        fig, axes = plt.subplots(n,n, figsize = (7,7))

        left  = 0.125  # the left side of the subplots of the figure
        right = 0.9    # the right side of the subplots of the figure
        bottom = 0.1   # the bottom of the subplots of the figure
        top = 0.9      # the top of the subplots of the figure
        wspace = 0.0   # the amount of width reserved for blank space between subplots
        hspace = 0.0   # the amount of height reserved for white space between subplots

        fig.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)

        if p == 'dz':
            z = self.dz
            cmap = cmr.get_sub_cmap('cmr.pride', 0.1, 0.9)
            vmin, vmax = -1.5, 1.5
            zlabel = '(z-z_{pz})/z'

        if p == 'beta':
            z = self.cat['diagnostics/beta']
            cmap = cmr.get_sub_cmap('cmr.pride', 0.1, 1.0)
            vmin, vmax = -3, 1
            zlabel = '\beta'

        norm = Normalize(vmin=vmin, vmax=vmax)

        cax = fig.add_axes((right, bottom, 0.02, top-bottom))
        cbar = fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax, orientation = 'vertical') # add the colourbar

        cax.yaxis.tick_right()
        cax.yaxis.set_label_position('right')
        cax.set_ylabel(rf'$\rm {zlabel}$')


        # ---- loop over parameters

        for i in np.arange(n):
            for j in np.arange(n):

                ax = axes[i,j]

                pi = parameters[i]
                pj = parameters[j]

                if i != n-1: ax.set_xticklabels([])

                if j != 0: ax.set_yticklabels([])

                if i == n-1:
                    ax.set_xlabel(r'${\rm'+parameter_labels[pj]+'}$')

                if j == 0 and i!=0:
                    ax.set_ylabel(r'${\rm'+parameter_labels[pi]+'}$')
                    # ax.get_yaxis().set_label_coords(-0.1*n,0.5)


                if j < i:

                    x = self.cat[f'parameters/{pj}'][()]
                    y = self.cat[f'parameters/{pi}'][()]

                    if scatter:
                        ax.scatter(x,y,s=5,alpha=0.5,color = cmap(norm(z)), lw=0)

                    else:

                        C, _,_,_ = binned_statistic_2d(x,y,z, statistic='median', bins=bins, range = [parameter_range[pj], parameter_range[pi]] )
                        # N, _,_ = np.histogram2d(z, log10FUV, bins=[z_bin_edges, log10L_bin_edges])
                        # s = np.array(N)<5
                        # C[s] = None
                        ax.imshow(C.T, aspect = 'auto', origin = 'lower', extent = (*parameter_range[pj], *parameter_range[pi]), vmin = vmin, vmax = vmax, cmap = cmap)


                    # xlims = [xe[0], xe[-1]]
                    # ylims = [ye[0], ye[-1]]

                    ax.set_xlim(parameter_range[pj])
                    ax.set_ylim(parameter_range[pi])

                elif j == i:

                    ax.set_axis_off()

                else:

                    ax.set_axis_off()

        fn = f'{model}/corner_{p}_{self.template_set}.pdf'
        print(fn)
        fig.savefig(fn)









    def make_pz_plots(self, id = False):

        if id:
            ids = [id]
        else:
            ids = self.ids

        for i in ids:

            z_a = self.cat['pz/eazy/z_a'][i]

            fig, ax = simple_plt()

            z = self.cat['pz/eazy/p_of_z/z'][:]
            pz = self.cat['pz/eazy/p_of_z/pz'][i, :]

            ax.plot(z, pz)
            ax.axvline(z_a, c='k', lw=1, alpha = 0.5)
            ax.set_xlabel(r'$\rm z$')
            ax.set_ylabel(r'$\rm P(z)dz$')

            fn = f'{model}/{i}_pz.pdf'
            print(fn)
            fig.savefig(fn)



    def make_sed_plots(self, id = False):

        if id:
            ids = [id]
        else:
            ids = self.ids

        for i in ids:

            z_a = self.cat['pz/eazy/z_a'][i]

            left  = 0.15
            bottom = 0.15
            width = 0.8
            height = 0.8

            fig = plt.figure(figsize = (3.5, 3.5))
            ax = fig.add_axes((left, bottom, width, height))

            template_norm = self.cat['pz/eazy/template_norms'][i]

            s = self.template_lam*(1+z_a)<50000

            for ii, t in enumerate(self.templates.values()):
                ax.plot(np.log10(t.lam[s]*(1+z_a)), template_norm[ii]*t.fnu[s], lw = 1, alpha = 0.2)

            total_fnu = np.sum(self.template_fnu*np.expand_dims(template_norm, 1), axis=0)
            ax.plot(np.log10(self.template_lam[s]*(1+z_a)), total_fnu[s], lw = 2, alpha = 0.3, c='k')

            for filter in self.filters:
                ax.scatter(np.log10(self.pivwv[filter]), self.cat[f'obs/{filter}/flux'][i], c='k', s=10)

            ax.set_xlim([3.5, 4.8])

            ax.set_xlabel(r'$\rm log_{10}(\lambda_{obs}/\AA)$')
            ax.set_ylabel(r'$\rm f_{\nu}/nJy$')

            fn = f'{model}/{i}_sed.pdf'
            print(fn)
            fig.savefig(fn)



if __name__ == "__main__":

    scenario = 'constant'
    z = 7.00
    template_set = 'tweak_fsps_QSF_12_v3'
    template_set = 'Larson22'
    model = f'data/{scenario}/{z:.2f}'

    for template_set in ['tweak_fsps_QSF_12_v3', 'Larson22', 'Wilkins22.bpass-2.2.1']:

        a = Analyser(model, template_set = template_set)
        # a.explore()
        # a.explore_icat()
        # a.dz_hist_plot()

        parameters = ['z','log10duration','fesc','log10tau_V','log10Z']

        # a.beta_dz_plot()

        a.corner(parameters, 'dz')
        # a.corner(parameters, 'beta')

        # id = 3

        # a.print_info(id)
        # a.make_pz_plots(id = id)
        # a.make_sed_plots(id = id)
