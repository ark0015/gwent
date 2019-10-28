import numpy as np

import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.image import NonUniformImage
import matplotlib.ticker as ticker
from matplotlib import cm

import astropy.constants as const
import astropy.units as u
from astropy.cosmology import z_at_value
from astropy.cosmology import WMAP9 as cosmo

def Plot_SNR(source,instrument,var_x,sample_x,var_y,sample_y,SNRMatrix,display=True,dl_axis=False,smooth_contours=False,figloc=None):
    """Plots the SNR contours from calcSNR

    Parameters
    ----------
    source : object
        Instance of a gravitational wave source class
    instrument : object
        Instance of a gravitational wave detector class
    var_x : str
        x-axis variable
    sample_x : array
        samples at which SNRMatrix was calculated corresponding to the x-axis variable
    var_y : str
        y-axis variable
    sample_y : array
        samples at which SNRMatrix was calculated corresponding to the y-axis variable
    SNRMatrix : array-like
        the matrix at which the SNR was calculated corresponding to the particular x and y-axis variable choices

    display : bool, optional
        Option to turn off display if saving multiple plots to a file
    dl_axis : bool, optional
        Option to turn on the right hand side labels of luminosity distance
    smooth_contours : bool, optional
        Option to interpolate contours to a finer mesh size to appear smooth instead of tiered contours
    figloc : str, optional
        If not None, saves the figure to the path figloc (ie. figloc = '/path/to/save/location/figname.type')

    """

    axissize = 16
    labelsize = 18
    textsize = 12
    textcolor = 'k'
    linesize = 4
    figsize = (10,8)

    colormap = 'viridis'
    logSNR = np.log10(SNRMatrix)

    logLevels_min = np.log10(np.array([5]))
    logLevels_max = np.ceil(np.amax(logSNR))
    if logLevels_max < logLevels_min:
        raise ValueError('All SNRs are lower than 5.')
    print_logLevels = np.append(logLevels_min,np.arange(1,logLevels_max+1))
    if smooth_contours:
        logLevels = np.linspace(logLevels_min,logLevels_max,100)[:,0]
    else:
        logLevels = print_logLevels

    if var_x in source.var_dict.keys():
        if isinstance(source.var_dict[var_x]['min'],u.Quantity):
            xlabel_min = source.var_dict[var_x]['min'].value
            xlabel_max = source.var_dict[var_x]['max'].value
        else:
            xlabel_min = source.var_dict[var_x]['min']
            xlabel_max = source.var_dict[var_x]['max']
    elif var_x in instrument.var_dict.keys():
        if isinstance(instrument.var_dict[var_x]['min'],u.Quantity):
            xlabel_min = instrument.var_dict[var_x]['min'].value
            xlabel_max = instrument.var_dict[var_x]['max'].value
        else:
            xlabel_min = instrument.var_dict[var_x]['min']
            xlabel_max = instrument.var_dict[var_x]['max']
    else:
        raise ValueError(var_x + ' is not a variable in the source or the instrument.')

    if var_y in source.var_dict.keys():
        if isinstance(source.var_dict[var_y]['min'],u.Quantity):
            ylabel_min = source.var_dict[var_y]['min'].value
            ylabel_max = source.var_dict[var_y]['max'].value
        else:
            ylabel_min = source.var_dict[var_y]['min']
            ylabel_max = source.var_dict[var_y]['max']
    elif var_y in instrument.var_dict.keys():
        if isinstance(instrument.var_dict[var_y]['min'],u.Quantity):
            ylabel_min = instrument.var_dict[var_y]['min'].value
            ylabel_max = instrument.var_dict[var_y]['max'].value
        else:
            ylabel_min = instrument.var_dict[var_y]['min']
            ylabel_max = instrument.var_dict[var_y]['max']
    else:
        raise ValueError(var_y + ' is not a variable in the source or the instrument.')

    #Can't take log of astropy variables
    if isinstance(sample_x,u.Quantity):
        sample_x = sample_x.value
    if isinstance(sample_y,u.Quantity):
        sample_y = sample_y.value

    #Set whether log or linearly spaced axes
    if var_x in ['M','z','L','A_acc']:
        xaxis_type = 'log'
    else:
        xaxis_type = 'lin'

    if var_y in ['M','z','L','A_acc']:
        yaxis_type = 'log'
    else:
        yaxis_type = 'lin'

    #########################
    #Make the Contour Plots
    fig, ax1 = plt.subplots(figsize=figsize)

    #Set axis scales based on what data sampling we used 
    if yaxis_type == 'lin' and xaxis_type == 'log':
        CS1 = ax1.contourf(np.log10(sample_x),sample_y,logSNR,logLevels,cmap = colormap)
        ax1.contour(np.log10(sample_x),sample_y,logSNR,print_logLevels,colors = 'k',alpha=1.0)
        ax1.set_xlim(np.log10(xlabel_min),np.log10(xlabel_max))
        ax1.set_ylim(ylabel_min,ylabel_max)
        x_labels = np.logspace(np.log10(xlabel_min),np.log10(xlabel_max),np.log10(xlabel_max)-np.log10(xlabel_min)+1)
        y_labels = np.linspace(ylabel_min,ylabel_max,ylabel_max-ylabel_min+1)
        ax1.set_yticks(y_labels)
        ax1.set_xticks(np.log10(x_labels))
    elif yaxis_type == 'log' and xaxis_type == 'lin':
        CS1 = ax1.contourf(sample_x,np.log10(sample_y),logSNR,logLevels,cmap = colormap)
        ax1.contour(sample_x,np.log10(sample_y),logSNR,print_logLevels,colors = 'k',alpha=1.0)
        ax1.set_xlim(xlabel_min,xlabel_max)
        ax1.set_ylim(np.log10(ylabel_min),np.log10(ylabel_max))
        x_labels = np.linspace(xlabel_min,xlabel_max,xlabel_max-xlabel_min+1)
        y_labels = np.logspace(np.log10(ylabel_min),np.log10(ylabel_max),np.log10(ylabel_max)-np.log10(ylabel_min)+1)
        ax1.set_xticks(x_labels)
        ax1.set_yticks(np.log10(y_labels))
    elif yaxis_type == 'lin' and xaxis_type == 'lin':
        CS1 = ax1.contourf(sample_x,sample_y,logSNR,logLevels,cmap = colormap)
        ax1.contour(sample_x,sample_y,logSNR,print_logLevels,colors = 'k',alpha=1.0)
        ax1.set_xlim(xlabel_min,xlabel_max)
        ax1.set_ylim(ylabel_min,ylabel_max)
        x_labels = np.linspace(xlabel_min,xlabel_max,xlabel_max-xlabel_min+1)
        y_labels = np.linspace(ylabel_min,ylabel_max,ylabel_max-ylabel_min+1)
        ax1.set_xticks(x_labels)
        ax1.set_yticks(y_labels)
    else:
        CS1 = ax1.contourf(np.log10(sample_x),np.log10(sample_y),logSNR,logLevels,cmap = colormap)
        ax1.contour(np.log10(sample_x),np.log10(sample_y),logSNR,print_logLevels,colors = 'k',alpha=1.0)
        ax1.set_xlim(np.log10(xlabel_min),np.log10(xlabel_max))
        ax1.set_ylim(np.log10(ylabel_min),np.log10(ylabel_max))
        x_labels = np.logspace(np.log10(xlabel_min),np.log10(xlabel_max),np.log10(xlabel_max)-np.log10(xlabel_min)+1)
        y_labels = np.logspace(np.log10(ylabel_min),np.log10(ylabel_max),np.log10(ylabel_max)-np.log10(ylabel_min)+1)
        ax1.set_yticks(np.log10(y_labels))
        ax1.set_xticks(np.log10(x_labels))

    #Set axes labels and whether log or linearly spaced
    if var_x == 'M':
        ax1.set_xlabel(r'$M_{\rm tot}$ $[\mathrm{M_{\odot}}]$',fontsize = labelsize)
        ax1.set_xticklabels([r'$10^{%i}$' %x if int(x) > 1 else r'$%i$' %(10**x) for x in np.log10(x_labels)],fontsize = axissize)
    elif var_x == 'q':
        x_labels = x_labels[::2]
        ax1.set_xticks(x_labels)
        ax1.set_xlabel(r'$\mathrm{Mass~Ratio}$',fontsize = labelsize)
        ax1.set_xticklabels([r'$%i$' %int(x) for x in x_labels],fontsize = axissize)
    elif var_x == 'z':
        ax1.set_xlabel(r'$\mathrm{Redshift}$',fontsize = labelsize)
        ax1.set_xticklabels([x if int(x) < 1 else int(x) for x in x_labels],fontsize = axissize)
    elif var_x in ['chi1','chi2']:
        print('here')
        x_labels = np.arange(round(xlabel_min*10),round(xlabel_max*10)+1,1)/10
        x_labels = x_labels[::2]
        ax1.set_xticks(x_labels)
        ax1.set_xlabel(r'$\mathrm{Spin}$',fontsize = labelsize)
        ax1.set_xticklabels([r'$%.1f$' %x for x in x_labels],fontsize = axissize)
    elif var_x == 'L':
        ax1.axvline(x=np.log10(2.5*u.Gm.to('m')),linestyle='--',color='k',label='Proposed arm length')
        ax1.set_xlabel(r'${\rm Armlength}$ $[\mathrm{m}]$',fontsize = labelsize)
        ax1.set_xticklabels([r'$10^{%i}$' %x if int(x) > 1 else r'$%i$' %(10**x) for x in np.log10(x_labels)],fontsize = axissize)
        ax1.legend(loc='lower right')
    elif var_x == 'A_acc':
        ax1.axvline(x=np.log10(3e-15),linestyle='--',color='k',label='Proposed Acceleration Noise Amplitude')
        ax1.set_xlabel(r'$\mathrm{Acceleration~Noise~Amplitude}$ $[\mathrm{m~s^{-2}}]$',fontsize = labelsize)
        ax1.set_xticklabels([r'$10^{%.0f}$' %x for x in np.log10(x_labels)],fontsize = axissize)
        ax1.legend(loc='lower left')
    elif var_x == 'A_IFO':
        ax1.axvline(x=np.log10(10e-12),linestyle='--',color='k',label='Proposed Optical Metrology Noise Amplitude')
        ax1.set_xlabel(r'$\mathrm{Optical~Metrology~Noise~Amplitude}~[\mathrm{m}]$',fontsize = labelsize)
        ax1.set_xticklabels([r'$10^{%.0f}$' %x for x in np.log10(x_labels)],fontsize = axissize)
        ax1.legend(loc='lower left')
    elif var_x == 'f_acc_break_low':
        scale = 10**round(np.log10(xlabel_min))
        x_labels = np.arange(round(xlabel_min/scale),round(xlabel_max/scale)+1,1)*scale
        ax1.set_xticks(x_labels)
        ax1.set_xlabel(r'$f_{\mathrm{acc,low}}$ $[\mathrm{mHz}]$',fontsize = labelsize)
        ax1.set_xticklabels([r'$%.1f$' %x for x in x_labels*1e3],fontsize = axissize)
    elif var_x == 'f_acc_break_high':
        scale = 10**round(np.log10(xlabel_min))
        x_labels = np.arange(round(xlabel_min/scale),round(xlabel_max/scale)+1,1)*scale
        ax1.set_xticks(x_labels)
        ax1.set_xlabel(r'$f_{\mathrm{acc,high}}$ $[\mathrm{mHz}]$',fontsize = labelsize)
        ax1.set_xticklabels([r'$%.1f$' %x for x in x_labels*1e3],fontsize = axissize)
    elif var_x == 'f_IFO_break':
        scale = 10**round(np.log10(xlabel_min))
        x_labels = np.arange(round(xlabel_min/scale),round(xlabel_max/scale)+1,1)*scale
        ax1.set_xticks(x_labels)
        ax1.set_xlabel(r'$f_{\mathrm{IFO,break}}$ $[\mathrm{mHz}]$',fontsize = labelsize)
        ax1.set_xticklabels([r'$%.1f$' %x for x in x_labels*1e3],fontsize = axissize)
    elif var_x == 'N_p':
        sample_range = max(x_labels)-min(x_labels)
        sample_rate = max(2,int(sample_range/10))
        x_labels = x_labels[::sample_rate]
        ax1.set_xticks(x_labels)
        ax1.set_xlabel(r'$\mathrm{Number~of~Pulsars}$',fontsize = labelsize)
        ax1.set_xticklabels([r'$%i$' %int(x) for x in x_labels],fontsize = axissize)
    elif var_x == 'cadence': 
        x_labels = np.arange(round(xlabel_min),round(xlabel_max)+1,5)
        ax1.set_xticks(x_labels)
        ax1.set_xlabel(r'$\mathrm{Observation~Cadence}$ $[\mathrm{yr}^{-1}]$',fontsize = labelsize)
        ax1.set_xticklabels([r'$%i$' %int(x) for x in x_labels],fontsize = axissize)
    elif var_x == 'sigma':
        scale = 10**round(np.log10(xlabel_min))
        x_labels = np.arange(round(xlabel_min/scale),round(xlabel_max/scale)+1,1)*scale
        ax1.set_xticks(x_labels)
        ax1.set_xlabel(r'$\mathrm{Pulsar~Timing~RMS}$ $[\mathrm{ns}]$',fontsize = labelsize)
        ax1.set_xticklabels([r'$%.0f$' %x for x in x_labels*1e9],fontsize = axissize)
    elif var_x == 'T_obs':
        x_labels = x_labels[::2]
        ax1.set_xticks(x_labels)
        ax1.set_xlabel(r'${\rm T_{obs}}$ $[\mathrm{yr}]$',fontsize = labelsize)
        ax1.set_xticklabels([r'$%i$' %int(x) for x in x_labels],fontsize = axissize)


    if var_y == 'M':
        ax1.set_ylabel(r'$M_{\rm tot}$ $[\mathrm{M_{\odot}}]$',fontsize = labelsize)
        ax1.set_yticklabels([r'$10^{%i}$' %y if int(y) > 1 else r'$%i$' %(10**y) for y in np.log10(y_labels)],fontsize = axissize)
    elif var_y == 'q':
        y_labels = y_labels[::2]
        ax1.set_yticks(y_labels)
        ax1.set_ylabel(r'$\mathrm{Mass~Ratio}$',fontsize = labelsize)
        ax1.set_yticklabels([r'$%i$' %int(y) for y in y_labels],fontsize = axissize)
    elif var_y == 'z':
        ax1.set_ylabel(r'$\mathrm{Redshift}$',fontsize = labelsize)
        ax1.set_yticklabels([y if int(y) < 1 else int(y) for y in y_labels],\
            fontsize = axissize)
    elif var_y in ['chi1','chi2']:
        print('here')
        y_labels = np.arange(round(ylabel_min*10),round(ylabel_max*10)+1,1)/10
        y_labels = y_labels[::2]
        ax1.set_yticks(y_labels)
        ax1.set_ylabel(r'$\mathrm{Spin}$',fontsize = labelsize)
        ax1.set_yticklabels([r'$%.1f$' %y for y in y_labels],fontsize = axissize)
    elif var_y == 'L':
        ax1.axhline(y=np.log10(2.5*u.Gm.to('m')),linestyle='--',color='k',label='Proposed arm length')
        ax1.set_ylabel(r'${\rm Armlength}$ $[\mathrm{m}]$',fontsize = labelsize)
        ax1.set_yticklabels([r'$10^{%i}$' %y if int(y) > 1 else r'$%i$' %(10**y) for y in np.log10(y_labels)],fontsize = axissize)
        ax1.legend(loc='lower right')
    elif var_y == 'A_acc':
        ax1.axhline(y=np.log10(3e-15),linestyle='--',color='k',label='Proposed Acceleration Noise Amplitude')
        ax1.set_ylabel(r'$\mathrm{Acceleration~Noise~Amplitude}$ $[\mathrm{m~s^{-2}}]$',fontsize = labelsize)
        ax1.set_yticklabels([r'$10^{%.0f}$' %y for y in np.log10(y_labels)],fontsize = axissize)
        ax1.legend(loc='lower left')
    elif var_y == 'A_IFO':
        ax1.axhline(y=np.log10(10e-12),linestyle='--',color='k',label='Proposed Optical Metrology Noise Amplitude')
        ax1.set_ylabel(r'$\mathrm{Optical~Metrology~Noise~Amplitude}~[\mathrm{m}]$',fontsize = labelsize)
        ax1.set_yticklabels([r'$10^{%.0f}$' %y for y in np.log10(y_labels)],fontsize = axissize)
        ax1.legend(loc='lower left')
    elif var_y == 'f_acc_break_low':
        scale = 10**round(np.log10(ylabel_min))
        y_labels = np.arange(round(ylabel_min/scale),round(ylabel_max/scale)+1,1)*scale
        ax1.set_yticks(y_labels)
        ax1.set_ylabel(r'$f_{\mathrm{acc,low}} [\mathrm{mHz}]$',fontsize = labelsize)
        ax1.set_yticklabels([r'$%.1f$' %y for y in y_labels*1e3],fontsize = axissize)
    elif var_y == 'f_acc_break_high':
        scale = 10**round(np.log10(ylabel_min))
        y_labels = np.arange(round(ylabel_min/scale),round(ylabel_max/scale)+1,1)*scale
        ax1.set_yticks(y_labels)
        ax1.set_ylabel(r'$f_{\mathrm{acc,high}} [\mathrm{mHz}]$',fontsize = labelsize)
        ax1.set_yticklabels([r'$%.1f$' %y for y in y_labels*1e3],fontsize = axissize)
    elif var_y == 'f_IFO_break':
        scale = 10**round(np.log10(ylabel_min))
        y_labels = np.arange(round(ylabel_min/scale),round(ylabel_max/scale)+1,1)*scale
        ax1.set_yticks(y_labels)
        ax1.set_ylabel(r'$f_{\mathrm{IFO,break}} [\mathrm{mHz}]$',fontsize = labelsize)
        ax1.set_yticklabels([r'$%.1f$' %y for y in y_labels*1e3],fontsize = axissize)
    elif var_y == 'N_p':
        sample_range = max(y_labels)-min(y_labels)
        sample_rate = max(2,int(sample_range/10))
        y_labels = y_labels[::sample_rate]
        ax1.set_yticks(y_labels)
        ax1.set_ylabel(r'$\mathrm{Number~of~Pulsars}$',fontsize = labelsize)
        ax1.set_yticklabels([r'$%i$' %int(y) for y in y_labels],fontsize = axissize)
    elif var_y == 'cadence': 
        y_labels = np.arange(round(ylabel_min),round(ylabel_max)+1,5)
        ax1.set_yticks(y_labels)
        ax1.set_ylabel(r'$\mathrm{Observation~Cadence}$ $[\mathrm{yr}^{-1}]$',fontsize = labelsize)
        ax1.set_yticklabels([r'$%i$' %int(y) for y in y_labels],fontsize = axissize)
    elif var_y == 'sigma':
        scale = 10**round(np.log10(ylabel_min))
        y_labels = np.arange(round(ylabel_min/scale),round(ylabel_max/scale)+1,1)*scale
        ax1.set_yticks(y_labels)
        ax1.set_ylabel(r'$\mathrm{Pulsar~Timing~RMS}$ $[\mathrm{ns}]$',fontsize = labelsize)
        ax1.set_yticklabels([r'$%.0f$' %y for y in y_labels*1e9],fontsize = axissize)
    elif var_y == 'T_obs':
        y_labels = y_labels[::2]
        ax1.set_yticks(y_labels)
        ax1.set_ylabel(r'${\rm T_{obs}}$ $[\mathrm{yr}]$',fontsize = labelsize)
        ax1.set_yticklabels([r'$%i$' %int(y) for y in y_labels],fontsize = axissize)

    ax1.yaxis.set_label_coords(-.10,.5)

    #If true, display luminosity distance on right side of plot
    if dl_axis:
        #Set other side y-axis for lookback time scalings
        ax2 = ax1.twinx()
        #Set axis scales based on what data sampling we used 
        if yaxis_type == 'lin' and xaxis_type == 'log':
            ax2.contour(np.log10(sample_x),sample_y,logSNR,print_logLevels,colors = 'k',alpha=1.0)
        elif yaxis_type == 'log' and xaxis_type == 'lin':
            ax2.contour(sample_x,np.log10(sample_y),logSNR,print_logLevels,colors = 'k',alpha=1.0)
        else:
            ax2.contour(np.log10(sample_x),np.log10(sample_y),logSNR,print_logLevels,colors = 'k',alpha=1.0)

        dists_min = cosmo.luminosity_distance(source.var_dict['z']['min']).to('Gpc')
        dists_min = np.ceil(np.log10(dists_min.value))
        dists_max = cosmo.luminosity_distance(source.var_dict['z']['max']).to('Gpc')
        dists_max = np.ceil(np.log10(dists_max.value))
        dists = np.arange(dists_min,dists_max)
        dists = 10**dists*u.Gpc

        distticks = [z_at_value(cosmo.luminosity_distance,dist) for dist in dists]
        #Set other side y-axis for lookback time scalings
        ax2.set_yticks(np.log10(distticks))
        #ax2.set_yticklabels(['%f' %dist for dist in distticks],fontsize = axissize)
        ax2.set_yticklabels([r'$10^{%i}$' %np.log10(dist) if np.abs(int(np.log10(dist))) > 1 else '{:g}'.format(dist) for dist in dists.value],fontsize = axissize)
        ax2.set_ylabel(r'$D_{L}$ [Gpc]',fontsize=labelsize)
        cbar = fig.colorbar(CS1,ax=(ax1,ax2),pad=0.01,ticks=print_logLevels)
    else:
        #Make colorbar
        cbar = fig.colorbar(CS1,ticks=print_logLevels)

    cbar.set_label(r'$SNR$',fontsize = labelsize)
    cbar.ax.tick_params(labelsize = axissize)
    cbar.ax.set_yticklabels([r'$10^{%i}$' %x if int(x) > 1 else r'$%i$' %(10**x) for x in print_logLevels])

    if display:
        plt.show()

    #Save Figure to File
    if figloc != None:
        fig.savefig(figloc,bbox_inches='tight')