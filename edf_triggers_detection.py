# -*- coding: utf-8 -*-
"""
Created on Fri Feb  1 13:01:53 2019

@author: GARDY
"""
import numpy as np
import matplotlib.pyplot as plt

def edf_triggers_detection(triggers_array, sfreq, trig_way = "up", plot_triggers = False, plot_timescale = "ms"):
    """
    Find triggers in MKR2+ macro electrode channel by computing the mean of the signal and
    establishing a threshold at 2* this mean. 
    
    Careful ! Sometimes 2 triggers can be detected for 1 real trigger because one is detected
    when the signal goes above the threshold and once again when it goes back under. 
    
    Parameters
    ----------
    triggers_array : np.array or list
        from MRK2+ channel
        
    sfreq : int
        in Hertz
        
    trig_way (optional) : str
        it will usually be "up"
    
    plot_triggers : bool
        if you want to plot (1) the MKR2+ channel with (2) the mean of the signal,
        (3) the threshold used to detect triggers and (3) the vertical lines
        showing the detected triggers.
    
    plot_timescale : str
        can be "s" for seconds or "ms" for miliseconds
    
    Returns
    --------
    edf_trigs : np.array
        array containing the detected triggers 
        (seconds)
    """    
    
    # get the mean (positive) the signal value of the triggers channel
    sup_0 = np.where(triggers_array > 0 )[0]
    mean_trig_value = np.mean(triggers_array[sup_0])
    
    # consider there is a trigger if a value is 2 * superior to the mean (positive) signal
    trig_detection = mean_trig_value * 3

    # get the triggers
    if trig_way == "up":
        edf_trigs = (np.array(np.where(triggers_array > trig_detection)) / sfreq)[0] # this value can change for each file, don't be surprised.
    else:
        edf_trigs = (np.array(np.where(triggers_array < trig_detection)) / sfreq)[0] # this value can change for each file, don't be surprised.
    
    print(edf_trigs)
    
    # plot the triggers channel
    if plot_triggers == True:
        if plot_timescale == "ms":
            x = np.arange(0,len(triggers_array)) / sfreq / 1000
            plt.plot(x,triggers_array)
            plt.hlines(mean_trig_value, xmin = 0, xmax = (len(triggers_array)/sfreq/1000), linewidth = 4, color = "red")
            plt.hlines(trig_detection, xmin = 0, xmax = (len(triggers_array)/sfreq/1000), linewidth = 4, color = "green")
            plt.xlabel("time (miliseconds)")
        elif plot_timescale == "s":
            x = np.arange(0,len(triggers_array)) / sfreq / 1000
            plt.plot(x,triggers_array)
            plt.hlines(mean_trig_value, xmin = 0, xmax = (len(triggers_array)/sfreq), linewidth = 4, color = "red")
            plt.hlines(trig_detection, xmin = 0, xmax = (len(triggers_array)/sfreq), linewidth = 4, color = "green")
            plt.xlabel("time (seconds)")
            
        plot_detected_triggers = False
        if plot_detected_triggers == True:
            if plot_timescale == "s":
                for i in edf_trigs:
                    plt.vlines(i, ymin = min(triggers_array), ymax = max(triggers_array), linewidth = 2, color = "purple")
            elif plot_timescale == "ms":
                for i in edf_trigs:
                    plt.vlines((i/1000), ymin = min(triggers_array), ymax = max(triggers_array), linewidth = 2, color = "purple")        
        plt.show()
    
    return(edf_trigs)
    
    
if __name__ == "__main__":    
    triggers_array = imagery_eeg.mne_raw._data[-2]
    sfreq = 256
    
    edf_trigs = edf_triggers_detection(triggers_array, sfreq, plot_triggers=True, plot_timescale="ms")
    edf_trigs = imagery_eeg.savetest2[87]
    
    plt.plot(triggers_array)
    plt.vlines(imagery_eeg.add_events.edf_trigs[87]*256, ymin = -1, ymax = 1)
    plt.show()