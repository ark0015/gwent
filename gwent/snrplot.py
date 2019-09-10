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
        if var_x == 'T_obs':
            xlabel_min = xlabel_min*u.yr.to('s')
            xlabel_max = xlabel_max*u.yr.to('s')
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
        if var_y == 'T_obs':
            ylabel_min = ylabel_min*u.yr.to('s')
            ylabel_max = ylabel_max*u.yr.to('s')
    else:
        raise ValueError(var_y + ' is not a variable in the source or the instrument.')

    #Can't take log of astropy variables
    if isinstance(sample_x,u.Quantity):
        sample_x = sample_x.value
    if isinstance(sample_y,u.Quantity):
        sample_y = sample_y.value

    #Set whether log or linearly spaced axes
    if var_x == 'q' or var_x == 'chi1' or var_x == 'chi2':
        xaxis_type = 'lin'
    else:
        xaxis_type = 'log'

    if var_y == 'q' or var_y == 'chi1' or var_y == 'chi2':
        yaxis_type = 'lin'
    else:
        yaxis_type = 'log'

    #########################
    #Make the Contour Plots
    fig, ax1 = plt.subplots(figsize=figsize)
    #Set other side y-axis for lookback time scalings
    ax2 = ax1.twinx()

    #Set axis scales based on what data sampling we used 
    if yaxis_type == 'lin' and xaxis_type == 'log':
        CS1 = ax1.contourf(np.log10(sample_x),sample_y,logSNR,logLevels,cmap = colormap)
        ax2.contour(np.log10(sample_x),sample_y,logSNR,print_logLevels,colors = 'k',alpha=1.0)
        ax1.set_xlim(np.log10(xlabel_min),np.log10(xlabel_max))
        ax1.set_ylim(ylabel_min,ylabel_max)
        x_labels = np.logspace(np.log10(xlabel_min),np.log10(xlabel_max),np.log10(xlabel_max)-np.log10(xlabel_min)+1)
        y_labels = np.linspace(ylabel_min,ylabel_max,ylabel_max-ylabel_min+1)
        ax1.set_yticks(y_labels)
        ax1.set_xticks(np.log10(x_labels))
    elif yaxis_type == 'log' and xaxis_type == 'lin':
        CS1 = ax1.contourf(sample_x,np.log10(sample_y),logSNR,logLevels,cmap = colormap)
        ax2.contour(sample_x,np.log10(sample_y),logSNR,print_logLevels,colors = 'k',alpha=1.0)
        ax1.set_xlim(xlabel_min,xlabel_max)
        ax1.set_ylim(np.log10(ylabel_min),np.log10(ylabel_max))
        x_labels = np.linspace(xlabel_min,xlabel_max,xlabel_max-xlabel_min+1)
        y_labels = np.logspace(np.log10(ylabel_min),np.log10(ylabel_max),np.log10(ylabel_max)-np.log10(ylabel_min)+1)
        ax1.set_xticks(x_labels)
        ax1.set_yticks(np.log10(y_labels))
    else:
        CS1 = ax1.contourf(np.log10(sample_x),np.log10(sample_y),logSNR,logLevels,cmap = colormap)
        ax2.contour(np.log10(sample_x),np.log10(sample_y),logSNR,print_logLevels,colors = 'k',alpha=1.0)
        ax1.set_xlim(np.log10(xlabel_min),np.log10(xlabel_max))
        ax1.set_ylim(np.log10(ylabel_min),np.log10(ylabel_max))
        x_labels = np.logspace(np.log10(xlabel_min),np.log10(xlabel_max),np.log10(xlabel_max)-np.log10(xlabel_min)+1)
        y_labels = np.logspace(np.log10(ylabel_min),np.log10(ylabel_max),np.log10(ylabel_max)-np.log10(ylabel_min)+1)
        ax1.set_yticks(np.log10(y_labels))
        ax1.set_xticks(np.log10(x_labels))

    #Set axes labels and whether log or linearly spaced
    if var_x == 'M':
        ax1.set_xlabel(r'$M_{\rm tot}$ $[M_{\odot}]$',fontsize = labelsize)
        ax1.set_xticklabels([r'$10^{%i}$' %x if int(x) > 1 else r'$%i$' %(10**x) for x in np.log10(x_labels)],fontsize = axissize)
    elif var_x == 'q':
        ax1.set_xlabel(r'$q$',fontsize = labelsize)
        ax1.set_xticklabels(x_labels,fontsize = axissize,rotation=45)
    elif var_x == 'z':
        ax1.set_xlabel(r'${\rm Redshift}$',fontsize = labelsize)
        ax1.set_xticklabels([x if int(x) < 1 else int(x) for x in x_labels],fontsize = axissize)
    elif var_x == 'chi1' or var_x == 'chi2':
        x_labels = np.arange(round(xlabel_min*10),round(xlabel_max*10)+1,1)/10
        ax1.set_xticks(x_labels)
        ax1.set_xlabel(r'${\rm Spin}$',fontsize = labelsize)
        ax1.set_xticklabels(x_labels,fontsize = axissize,rotation=45)
        ax1.legend(loc='lower right')
    elif var_x == 'L':
        ax1.axvline(x=np.log10(2.5*u.Gm.to('m')),linestyle='--',color='k',label='Proposed arm length')
        ax1.set_xlabel(r'${\rm Armlength}$ $[m]$',fontsize = labelsize)
        ax1.set_xticklabels([r'$10^{%i}$' %x if int(x) > 1 else r'$%i$' %(10**x) for x in np.log10(x_labels)],fontsize = axissize)
    elif var_x == 'T_obs':
        ax1.set_xlabel(r'${\rm T_{obs}}$ $[yr]$',fontsize = labelsize)
        ax1.set_xticklabels([r'$%i$' %int(x) for x in x_labels/u.yr.to('s')],fontsize = axissize)

    if var_y == 'M':
        ax1.set_ylabel(r'$M_{\rm tot}$ $[M_{\odot}]$',fontsize = labelsize)
        ax1.set_yticklabels([r'$10^{%i}$' %y if int(y) > 1 else r'$%i$' %(10**y) for y in np.log10(y_labels)],fontsize = axissize)
    elif var_y == 'q':
        ax1.set_ylabel(r'$q$',fontsize = labelsize)
        ax1.set_yticklabels(y_labels,fontsize = axissize,rotation=45)
    elif var_y == 'z':
        ax1.set_ylabel(r'${\rm Redshift}$',fontsize = labelsize)
        ax1.set_yticklabels([y if int(y) < 1 else int(y) for y in y_labels],\
            fontsize = axissize)
    elif var_y == 'chi1' or var_y == 'chi2':
        y_labels = np.arange(round(ylabel_min*10),round(ylabel_max*10)+1,1)/10
        ax1.set_yticks(y_labels)
        ax1.set_ylabel(r'${\rm Spin}$',fontsize = labelsize)
        ax1.set_yticklabels(y_labels,fontsize = axissize,rotation=45)
    elif var_y == 'L':
        ax1.axhline(y=np.log10(2.5*u.Gm.to('m')),linestyle='--',color='k',label='Proposed arm length')
        ax1.set_ylabel(r'${\rm Armlength}$ $[m]$',fontsize = labelsize)
        ax1.set_yticklabels([r'$10^{%i}$' %y if int(y) > 1 else r'$%i$' %(10**y) for y in np.log10(y_labels)],fontsize = axissize)
        ax1.legend(loc='lower right')
    elif var_y == 'T_obs':
        ax1.set_ylabel(r'${\rm T_{obs}}$ $[yr]$',fontsize = labelsize)
        ax1.set_yticklabels([r'$%i$' %int(y) for y in y_labels/u.yr.to('s')],fontsize = axissize)

    ax1.yaxis.set_label_coords(-.10,.5)

    #If true, display luminosity distance on right side of plot
    if dl_axis:
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
        #Remove y-axis labels
        ax2.tick_params(axis='y',right=False,labelright=False)
    cbar.set_label(r'$SNR$',fontsize = labelsize)
    cbar.ax.tick_params(labelsize = axissize)
    cbar.ax.set_yticklabels([r'$10^{%i}$' %x if int(x) > 1 else r'$%i$' %(10**x) for x in print_logLevels])

    if display:
        plt.show()

    #Save Figure to File
    if figloc != None:
        fig.savefig(figloc,bbox_inches='tight')