# -*- coding: utf-8 -*-
"""
Created on Fri Feb  1 14:32:15 2019

@author: GARDY
"""

# Imports from computer
from load_signal import io_eeg_to_mne
from add_events import add_events_to_mne_sig
import mne

if __name__ == "__main__":    

    # Load eeg data    
    filepath = r"F:\GardyL\Data_storage\Imagery\IMAGERY_AB28\MACRO\Bipolaire\Imagery1.edf"

    mne_raw = io_eeg_to_mne(filepath = filepath, rescale_voltage = True, downsampling_freq = 256)[0]
    mne_raw._data.shape
    
    # Add events to the eeg data
    events_path = r"F:\GardyL\Data_storage\Imagery\IMAGERY_AB28\Triggers1.csv"

    add_events = add_events_to_mne_sig(mne_raw, data_type = ".edf")
    add_events.import_events(mne_raw, events_path)
    add_events.eventsTime_to_eventsTimestamps(mne_raw)
    add_events.delete_false_triggers(mne_raw)
    add_events.microMacro_realign_eventstimes(mne_raw, plot_triggers=False)
    add_events.add_events_to_data(mne_raw)    

    # Plot the eeg data
    scalings = {'mag': 2, 'grad': 2}
        
    mne_raw.plot(scalings=scalings, title='Data from arrays', show=True, 
                      block=True, events = mne.find_events(mne_raw), event_id = None,
                      duration = 10, n_channels = 40, lowpass = None,
                      highpass = None, order = None)        
    
