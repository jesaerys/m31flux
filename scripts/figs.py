import astropy.io.fits
import matplotlib.pyplot as plt
import numpy as np
import os


REPO_DIR = '/Users/Jake/Research/code/m31flux'


def log10(val):
    return np.where(val > 0, val, np.nan)


def plot_map(data, outfile, label, limits, stretch, cmap):
    fig_dx = 6.0
    ax2_dy = 0.15

    data = data[::-1].T  # rotate to landscape
    ny, nx = data.shape

    plt.close()
    h2 = float(ny)/nx * fig_dx
    s1, s2 = 0.5, 0.1
    h1 = s1 + ax2_dy + s2
    fig_dy = h1 + h2

    fig = plt.figure(figsize=(fig_dx, fig_dy))
    ax1 = fig.add_axes((0, h1/fig_dy, 1, h2/fig_dy))

    # Stretch
    data = np.clip(data, *limits)
    data = stretch(data)

    # Plot
    img = ax1.imshow(data, origin='lower', interpolation='nearest', cmap=cmap)
    ax1.xaxis.set_visible(False)
    ax1.yaxis.set_visible(False)
    for spine in ax1.spines.values():
        spine.set_visible(False)

    ax2 = fig.add_axes((0.02, s1/fig_dy, 0.96, ax2_dy/fig_dy))
    cb = fig.colorbar(img, cax=ax2, orientation='horizontal')
    ax2.tick_params(labelsize=10)
    ax2.set_xlabel(label, size=10)

    # Write
    fig.savefig(outfile)
    plt.close()


def plot_scatter(x, y, xlabel, ylabel, outfile):
    fig_dx = 4.0

    fig = plt.figure(figsize=(fig_dx, fig_dx))
    ax = fig.add_axes((0.15, 0.15, 0.80, 0.80))
    ax.set_rasterization_zorder(0)

    plt.axhline(0, color='0.6', ls='--')
    plt.plot(x, y, 'k.', mec='none', ms=3, alpha=0.2, zorder=-1)

    ax.set_xlabel(xlabel, size=10)
    ax.set_ylabel(ylabel, size=10)
    ax.tick_params(labelsize=10)

    # Write
    fig.savefig(outfile, dpi=200)
    plt.close()


def main():
    limits = (1e-16, 7e-15)
    cmap = plt.cm.gist_heat_r
    weights_file = os.path.join(REPO_DIR, 'maps', 'weights.fits')
    weights = astropy.io.fits.getdata(weights_file)


    filename = os.path.join(REPO_DIR, 'maps', 'galex_fuv.fits')
    galex_fuv = astropy.io.fits.getdata(filename) * weights
    plotname = os.path.join(REPO_DIR, 'figs', 'galex_fuv.pdf')
    label = (r'$f_\mathrm{FUV,obs} \, '
             r'\mathrm{(erg \,s^{-1} \,cm^{-2} \,\AA^{-1})}$')
    plot_map(galex_fuv, plotname, label, limits, log10, cmap)

    filename = os.path.join(REPO_DIR, 'maps', 'mod_fuv_red.fits')
    mod_fuv_red = astropy.io.fits.getdata(filename)
    plotname = os.path.join(REPO_DIR, 'figs', 'mod_fuv_red.pdf')
    label = (r'$f_\mathrm{FUV,mod} \, '
             r'\mathrm{(erg \,s^{-1} \,cm^{-2} \,\AA^{-1})}$')
    plot_map(mod_fuv_red, plotname, label, limits, log10, cmap)

    filename = os.path.join(REPO_DIR, 'maps', 'mod_fuv_int.fits')
    mod_fuv_int = astropy.io.fits.getdata(filename)
    plotname = os.path.join(REPO_DIR, 'figs', 'mod_fuv_int.pdf')
    label = (r'intrinsic $f_\mathrm{FUV,mod} \, '
             r'\mathrm{(erg \,s^{-1} \,cm^{-2} \,\AA^{-1})}$')
    plot_map(mod_fuv_int, plotname, label, limits, log10, cmap)

    fuv_ratio = mod_fuv_red / galex_fuv
    plotname = os.path.join(REPO_DIR, 'figs', 'fuv_ratio.pdf')
    label = r'$\log_{10}(f_\mathrm{FUV,mod} / f_\mathrm{FUV,obs})$'
    ratio_limits = (0.1, 10)
    ratio_cmap = plt.cm.RdBu_r
    plot_map(fuv_ratio, plotname, label, ratio_limits, np.log10, ratio_cmap)

    xlabel = (r'$\log_{10}(f_\mathrm{FUV,obs} / '
              r'\mathrm{(erg \,s^{-1} \,cm^{-2} \,\AA^{-1})})$')
    ylabel = r'$\log_{10}(f_\mathrm{FUV,mod} / f_\mathrm{FUV,obs})$'
    plotname = os.path.join(REPO_DIR, 'figs', 'fuv_ratio_plot.pdf')
    plot_scatter(np.log10(galex_fuv), np.log10(fuv_ratio),
                 xlabel, ylabel, plotname)


    filename = os.path.join(REPO_DIR, 'maps', 'galex_nuv.fits')
    galex_nuv = astropy.io.fits.getdata(filename) * weights
    plotname = os.path.join(REPO_DIR, 'figs', 'galex_nuv.pdf')
    label = (r'$f_\mathrm{NUV,obs} \, '
             r'\mathrm{(erg \,s^{-1} \,cm^{-2} \,\AA^{-1})}$')
    plot_map(galex_nuv, plotname, label, limits, log10, cmap)

    filename = os.path.join(REPO_DIR, 'maps', 'mod_nuv_red.fits')
    mod_nuv_red = astropy.io.fits.getdata(filename)
    plotname = os.path.join(REPO_DIR, 'figs', 'mod_nuv_red.pdf')
    label = (r'$f_\mathrm{NUV,mod} \, '
             r'\mathrm{(erg \,s^{-1} \,cm^{-2} \,\AA^{-1})}$')
    plot_map(mod_nuv_red, plotname, label, limits, log10, cmap)

    filename = os.path.join(REPO_DIR, 'maps', 'mod_nuv_int.fits')
    mod_nuv_int = astropy.io.fits.getdata(filename)
    plotname = os.path.join(REPO_DIR, 'figs', 'mod_nuv_int.pdf')
    label = (r'intrinsic $f_\mathrm{NUV,mod} \, '
             r'\mathrm{(erg \,s^{-1} \,cm^{-2} \,\AA^{-1})}$')
    plot_map(mod_nuv_int, plotname, label, limits, log10, cmap)

    nuv_ratio = mod_nuv_red / galex_nuv
    plotname = os.path.join(REPO_DIR, 'figs', 'nuv_ratio.pdf')
    label = r'$\log_{10}(f_\mathrm{NUV,mod} / f_\mathrm{NUV,obs})$'
    ratio_limits = (0.1, 10)
    ratio_cmap = plt.cm.RdBu_r
    plot_map(nuv_ratio, plotname, label, ratio_limits, np.log10, ratio_cmap)

    xlabel = (r'$\log_{10}(f_\mathrm{NUV,obs} / '
              r'\mathrm{(erg \,s^{-1} \,cm^{-2} \,\AA^{-1})})$')
    ylabel = r'$\log_{10}(f_\mathrm{NUV,mod} / f_\mathrm{NUV,obs})$'
    plotname = os.path.join(REPO_DIR, 'figs', 'nuv_ratio_plot.pdf')
    plot_scatter(np.log10(galex_nuv), np.log10(nuv_ratio),
                 xlabel, ylabel, plotname)


if __name__ == '__main__':
    main()
