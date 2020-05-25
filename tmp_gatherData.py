# -*- coding: utf-8 -*-
"""
Created on Tue Feb  5 15:31:15 2019

@author: gardy
"""

# General imports
import mne
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Imports from computer
from imagery_main import *

def merge_datasets(blocks_dict):

    res = []
    
    for bl in blocks_dict:
        for blist in blocks_dict[bl][0]:
            # Create object
            imagery_eeg = imagery_sig(patient_ID = "AB28")
            
            # Load imagery data   
            imagery_eeg.load_data(filepath = blocks_dict[bl][1])
        
            # Add events to the imagery eeg data
            imagery_eeg.import_triggers(events_path = blocks_dict[bl][2], remove_events_after_timestamp = blocks_dict[bl][3])
            
            # Load behavioural data
            behaviour_path = r"F:\GardyL\Data_storage\Imagery\IMAGERY_AB28\AB28_Imagery.xlsx"
            imagery_eeg.import_behaviour(behaviour_path, bloc_number = [blist])
        
            # Find blocs beginning & end
            imagery_eeg.find_bloc_EndBeginning()
            
            # Add the triggers semantics
            imagery_eeg.triggers_semantic()    
            
            res.append(imagery_eeg.triggers_dataframe_dict)
        
        res[0]["bloc_1"].keys()

    return(res)


if __name__ == "__main__": 
    
    blocks_dict = {}
    blocks_dict["bloc_1"] = [
            [1],
            r"F:\GardyL\Data_storage\Imagery\IMAGERY_AB28\MACRO\Bipolaire\Imagery1.edf",
            r"F:\GardyL\Data_storage\Imagery\IMAGERY_AB28\Triggers1.csv",
            []
                ]
    
    blocks_dict["bloc_234"] = [
            [2,3,4],
            r"F:\GardyL\Data_storage\Imagery\IMAGERY_AB28\MACRO\Bipolaire\Imagery234.edf",
            r"F:\GardyL\Data_storage\Imagery\IMAGERY_AB28\Triggers2.csv",
            1450
            ]

myres = merge_datasets(blocks_dict)







        
    # Plot the eeg data
    scalings = {'mag': 2, 'grad': 2}
        
    imagery_eeg.mne_raw.plot(scalings=scalings, title='Data from arrays', show=True, 
                      block=True, events = mne.find_events(imagery_eeg.mne_raw, shortest_event=0), event_id = None,
                      duration = 10, n_channels = 40, lowpass = None,
                      highpass = None, order = None)        
    