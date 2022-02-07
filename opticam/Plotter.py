import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import rgb2hex
from matplotlib import cm


def makeplots(Observation, type):
    """Function to visualize the output of the ``opticam.Observation`` class.

    This function has different outputs based on whether the instrument used is a 
    spectrograph or an imager. For an imager, the function will output a bar plot 
    that shows the S/N per bandpass and a plot that shows the different
    noise sources as a function of wavelength. If it is a spectrograph, then the 
    function will plot the S/N per wavelength and it also plots the noise sources 
    as a function of wavelength.

    Parameters
    -----------
    Observation : object
        The ``opticam.Observation``class.

    """

    if type == 'SN':
        DATA = Observation.SN
        exptime = str(Observation.exptime)
        if Observation.isImager == 1:

            filterNum = len(DATA)

            # Set colors
            cmap = cm.get_cmap("rainbow", filterNum)
            hexvals = []
            for i in range(cmap.N):
                rgb = cmap(i)[:3]
                hexvals.append(rgb2hex(rgb))

            width = 0.35  # width of SN/exptime plot bars
            width2 = 0.3  # width of noise plot bars

            filterSN = [row[0] for row in DATA]

            filter_names = [row.strip('_filter') for row in [row[1] for row in DATA]]
            filter_names = [row.replace('prime', "'") for row in filter_names]

            ind = np.arange(filterNum)  # the x locations for the bars

            fig = plt.figure(figsize=(10.0, 10.0))
            fig.suptitle(r"Observation Parameters for a " + exptime + " s Exposure",
                             y=0.95, weight="bold", fontsize=25)
            #plt.subplots_adjust(hspace=0.3)
            ax1 = fig.add_subplot(1, 1, 1)
            # ax2 = fig.add_subplot(2, 1, 2)

            ax1.bar(ind, filterSN, color=hexvals, edgecolor="black", width=width)
            ax1.set_xlabel(r"$\lambda$ ( $\AA$ )", fontsize=18)
            ax1.set_ylabel(r"$\frac{S}{N}$", fontsize=18, rotation=0, labelpad=20)
            ax1.set_title(r"$\frac{S}{N}$  for " + exptime + " seconds", fontsize=20)
            ax1.set_xticks(ind)
            ax1.set_xticklabels(filter_names, fontsize=18)
            ax1.tick_params(axis="y", labelsize=13)
            ax1.set_ylim(0, 1.2 * np.max(filterSN))

            for i, v in enumerate(filterSN):
                ax1.text(ind[i], v + 1, str(round(v, 2)), horizontalalignment="center",
                         verticalalignment="bottom",
                         fontsize=15)
            ax1.grid(True)

            plt.show()
            # plotname = "SNfromTime" + exptime + "s.png"
            # plotname2 = "Noise" + plotname[2:]
            # plt.savefig(plotname)

        if Observation.isImager == 0:
            dipersionNum = len(DATA)
            Wavelengths = [row[0] for row in DATA]
            SNs = [row[1] for row in DATA]
            dispersionNames = [row[2] for row in DATA]

            median = []
            fig = plt.figure(figsize=(10.0, 10.0))
            fig.suptitle(r"Observation Parameters for a " + exptime + " s Exposure", 
                        y=0.95, weight="bold", fontsize=25)
           # plt.subplots_adjust(hspace=0.3)
            ax1 = fig.add_subplot(1, 1, 1)
            # ax2 = fig.add_subplot(2, 1, 2)

            for i, row in enumerate(dispersionNames):
                x = Wavelengths[i]
                y1 = np.nan_to_num(SNs[i])
                median.append(np.median(y1))
                ax1.plot(x, y1, label=row)

            ax1.set_xlabel(r"$\lambda$ ( $\AA$ )", fontsize=18)
            ax1.set_ylabel(r"$\frac{S}{N}$", fontsize=18, rotation=0, labelpad=20)
            ax1.set_title(r"$\frac{S}{N}$  vs.  $\lambda$  for " + exptime + " seconds", fontsize=20)
            ax1.tick_params(axis="both", labelsize=13)
            ax1.set_ylim(0, 3 * np.max(median))
            ax1.legend()
            ax1.grid(True)

    if type == 'Time':
        DATA = Observation.Time
        SN = str(Observation.SigToNoise)
        if Observation.isImager == 1:
            filterNum = len(DATA)
            ind = np.arange(filterNum)  # the x locations for the bars

            # Set colors
            cmap = cm.get_cmap("rainbow", filterNum)
            hexvals = []
            for i in range(cmap.N):
                rgb = cmap(i)[:3]
                hexvals.append(rgb2hex(rgb))

            width = 0.35  # width of SN/exptime plot bars
            filterTime = [row[0] for row in DATA]

            filter_names1 = [row.replace('_time','') for row in [row[1] for row in DATA]]
            filter_names = [row.replace('prime', "'") for row in filter_names1]
            fig = plt.figure(figsize=(10.0, 10.0))
            fig.suptitle(r"Observation Parameters for S/N = " + SN,
                         y=0.95, weight="bold", fontsize=25)
            # plt.subplots_adjust(hspace=0.3)
            ax1 = fig.add_subplot(1, 1, 1)
            #ax2 = fig.add_subplot(2, 1, 2)

            ax1.bar(ind, filterTime, color=hexvals, edgecolor="black", width=width)
            ax1.set_xlabel(r"$\lambda$ ( $\AA$ )", fontsize=18)
            ax1.set_ylabel(r"$t$ ($s$)", fontsize=18, rotation=0, labelpad=25)
            ax1.set_title(r"Exposure Time for  $\frac{S}{N}=$" + SN, fontsize=20)
            ax1.set_xticks(ind)
            ax1.set_xticklabels(filter_names, fontsize=18)
            ax1.tick_params(axis="y", labelsize=13)
            ax1.set_ylim(0, 1.2 * np.max(filterTime))
            for i, v in enumerate(filterTime):
                ax1.text(ind[i], v + 1, str(round(v, 2)), horizontalalignment="center",
                         verticalalignment="bottom",
                         fontsize=15)
            ax1.grid(True)

            plt.show()
        if Observation.isImager == 0:
            dipersionNum = len(DATA)
            Wavelengths = [row[0] for row in DATA]
            Times = [row[1] for row in DATA]
            dispersionNames = [row[2] for row in DATA]

            median = []
            fig = plt.figure(figsize=(10.0, 10.0))
            fig.suptitle(r"Observation Parameters for  SN = " + SN, y=0.95,
                         weight="bold", fontsize=25)
            # plt.subplots_adjust(hspace=0.3)
            ax1 = fig.add_subplot(1, 1, 1)
            #ax2 = fig.add_subplot(2, 1, 2)

            for i, row in enumerate(dispersionNames):
                x = Wavelengths[i]
                y1 = np.nan_to_num(Times[i])
                median.append(np.median(y1))
                ax1.scatter(x, y1, label=row)

            ax1.set_xlabel(r"$\lambda$ ( $\AA$ )", fontsize=18)
            ax1.set_ylabel(r"$t$ ($s$)", fontsize=18, rotation=0, labelpad=25)
            ax1.set_title(r"Exposure Time vs.  $\lambda$  for  $\frac{S}{N}=$" + SN, fontsize=20)
            ax1.tick_params(axis="both", labelsize=13)
            ax1.set_ylim(0, 2 * np.max(median))
            ax1.legend()
            ax1.grid(True)
            plt.show()

    return fig