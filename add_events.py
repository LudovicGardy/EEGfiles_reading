# -*- coding: utf-8 -*-
"""
Created on Wed Jan  2 11:46:51 2019

@author: GARDY
"""
# General imports
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import mne

# Imports from computer
from load_signal import io_eeg_to_mne
from edf_triggers_detection import edf_triggers_detection

class add_events_to_mne_sig:
    def __init__(self, mne_raw, filepath, rescale_voltage = True):
        self.events_dict = {}
        self.file_ext = filepath[-4:]
        self.sfreq = mne_raw.info["sfreq"]
        self.time_stamps = len(mne_raw._data[0])
        self.recording_duration_seconds = round(self.time_stamps / self.sfreq, 2)
        self.recording_duration_minutes = round(self.recording_duration_seconds / 60, 2)

        print("--")
        print("raw data shape: ", mne_raw._data.shape)
        print("")
        print("\n\nrecording duration = {} seconds ({} minutes)".format(self.recording_duration_seconds, self.recording_duration_minutes))
        print("..")

    def import_events(self, mne_raw, events_path):
        """
        Events are already loaded in .trc files (done in utils_eeg_io.py).
        Events must be loaded here for .edf files.
        
        This method can be used to import any kind of events.
        
        Parameters
        ----------
        mne_raw : mne signal
            mne structure
        
        events_path : str
            
        Returns
        --------
        self.events_array_sec : numpy 1D array
            contains all the event times (sec)
            
        self.event_index : numpy 1D array
            contains all the event indexes (range)
            
        self.nb_events : int
            amount of events           
        """        
        
        if self.file_ext == ".edf":
            if events_path[-4:] == ".csv":
                self.events_dict = pd.read_csv(events_path)
            else:
                self.events_dict = pd.read_excel(events_path)  
                
            self.events_array_sec = np.array(self.events_dict[list(self.events_dict.keys())[0]])
            self.nb_events = len(self.events_array_sec)
            self.event_index = np.arange(1,(self.nb_events + 1)) 
            
        elif self.file_ext == ".TRC":
            self.events_array_sec = mne.find_events(mne_raw)
        
        print("--")
        print("raw data shape: ", mne.find_events(mne_raw).shape)
        print("events_array_sec length = ", len(self.events_array_sec))
        print("..")

    def delete_false_triggers(self, mne_raw):
        """
        Delte 1 trigger on 2 because of the double-triggers:
        One trigger appears when the white square appears and an other when the
        white square disappears.
        
        Parameters
        ----------
        mne_raw : mne signal
            mne structure

        Returns
        --------
        self.events_array_sec : numpy 1D array (UPDATE)
            contains all the event times (sec)
            
        self.event_index : numpy 1D array (UPDATE)
            contains all the event indexes (range)
            
        self.nb_events : int
            amount of events           
        """        
        
        tmp = []

        for i in range(len(self.events_array)):
            if np.mod(i, 2) == 0:
                tmp.append(self.events_array[i])
        self.events_array = np.array(tmp)      
                
        print("--")
        print("events_array (time stamps) length is now: ", len(self.events_array))
        print("..")
       
    def microMacro_realign_eventstimes(self, mne_raw, trig_way = "up", plot_triggers = True):
        """
        Align micro triggers on the macro triggers. Correct the events_array 
        by adding the shift between these two values.

        Parameters
        ----------
        mne_raw : mne signal
            mne structure

        trig_way (optional) : str
        
        plot_triggers (optional) : bool
        
        Returns
        --------
        self.events_array : numpy 1D array (UPDATE)
            contains all the event times (sec)
        """     

        print("--")
        
        # find the triggers channel
        trigChanl = np.where(np.array(mne_raw.info["ch_names"]) == "MKR2+")[0]
        
        edf_trigs = edf_triggers_detection(triggers_array = mne_raw._data[trigChanl][0], sfreq = self.sfreq)
        print("")
        print("edf_trigs length = ", len(edf_trigs))
        print("")
        
        self.check_edf_trigs = edf_trigs
        
        # calculate the shift between the 1st [mcaro.edf] trig and the 1st [micro.ns5] trig
        self.first_edf_trig = edf_trigs[0] * self.sfreq
        shift = abs(self.events_array[0] - self.first_edf_trig)
        print("first edf trig found at: ", self.first_edf_trig)
        print("first event recorded at: ", self.events_array[0])
        print("the shift between these two is: ", shift)

        # realign the micro triggers on the macro triggers
        # realign the micro triggers on the macro triggers
        if self.events_array[0] > self.first_edf_trig:
            self.events_array = self.events_array - shift
        elif self.events_array[0] < self.first_edf_trig:
            self.events_array = self.events_array + shift
        print("")
        print("first event recorded corrected = ", self.events_array[0])   
        print("--> 1st event is now aligned with 1st edf trigger")
        print("..")        

        
    def eventsTime_to_eventsTimestamps(self, mne_raw):    
        """
        Transfrom the <event time> in <event timestamp>.

        Parameters
        ----------
        mne_raw : mne signal
            mne structure
        
        Returns
        --------            
        self.events_array_timestamp : numpy 1D array
            contains all the event times stamps
            
        self.events_array : numpy 1D array 
            contains all the event times stamps
        """     
        
        if self.file_ext != ".TRC":
            self.events_array_timestamp = self.events_array_sec * self.sfreq
            self.events_array = self.events_array_sec * self.sfreq

        elif self.file_ext == ".TRC":
            for i in range(len(self.events_array_sec[0])):
                # from miliseconds to seconds
                self.events_array_sec[i][0] = self.events_array_sec[i][0] * 1e-3
                # from time to time_stamp
                self.events_array_sec[i][0] = self.events_array_sec[i][0] * self.sfreq

        print("--")
        print("events_array (time stamps) length = ", len(self.events_array))       
        print("..")

    def eventstime_periods_toKeep(self, mne_raw):
        
        def onclick(event):
            print('%s click: button=%d, x=%d, y=%d, xdata=%f, ydata=%f' %
                  ('double' if event.dblclick else 'single', event.button,
                   event.x, event.y, event.xdata, event.ydata))
        
            
            if event.dblclick:
                ax.vlines(x = event.xdata, ymin = min(y_data), ymax = max(y_data), linestyle = "--", color = "magenta")

                count.append(event.xdata)
                if len(count) % 2 != 0:
                    ax.text(event.xdata,event.ydata,"begin \n({})".format(int(round(event.xdata,2))))
                    self.periods_dict["beginnings"].append(int(round(event.xdata,4)))
                else:
                    ax.text(event.xdata,event.ydata,"end \n({})".format(int(round(event.xdata,2))))
                    self.periods_dict["endings"].append(int(round(event.xdata,4)))
        
                fig.canvas.draw()
                fig.canvas.flush_events()

        self.periods_dict = {}
        self.periods_dict["beginnings"] = []
        self.periods_dict["endings"] = []
        
        event_channel = len(mne_raw._data) - 2
        x_data = np.arange(0,self.time_stamps)
        x_sec = x_data / self.sfreq
        y_data = mne_raw._data[event_channel]
        
        plt.ion()
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(x_sec, y_data)
        ax.set_xlabel("time (sec), sfreq = {}".format(self.sfreq))
        ax.set_title("[double click] for the beginning, [double cick] for the end")
        fig.suptitle("select the tie time periods you want to keep (events)")
        
        count = []
            
        cid = fig.canvas.mpl_connect('button_press_event', onclick)
        plt.show()               

    def remove_bad_events(self):
        
        loopsize = len(self.periods_dict["beginnings"])
        
        self.clean_array = []
        
        for i in range(0,loopsize):
            b = self.periods_dict["beginnings"][i] * self.sfreq
            e = self.periods_dict["endings"][i] * self.sfreq
        
            _b = np.where(self.events_array > b)[0][0]
            _e = np.where(self.events_array < e)[0][-1]
            
            self.clean_array.append(self.events_array[_b:_e])

        print("--")
        print("events array is now clean.")
        print("nb blocs = {}".format(len(self.clean_array)))
        
        print("length of events array was = {}".format(len(self.events_array)))        
        print("..")
        
    def cleanArray_to_eventsArray(self):
        
        try:
            print("--")
            self.events_array = np.concatenate(self.clean_array)
            
            print("length of events array is now = {}".format(len(self.events_array)))        
            print("..")
     
            #self.nb_events = len(self.events_array)
            #self.event_index = np.arange(1,(self.nb_events + 1))     
        except NameError:
            print("")
            print("Events array is probably already clean")
            print("")
        
    def add_events_to_data(self, mne_raw):
        """
        Now that the events are clean and transformed into timestamps,
        add it to the mne structure.
        
        Parameters
        ----------
        mne_raw : mne signal
            mne structure
        
        Returns
        --------            
        self.events_forplot: numpy 1D array
            mne.find_events(mne_raw)
        """    
        self.nb_events = len(self.events_array)
        self.event_index = np.arange(1,(self.nb_events + 1))     
       
        if self.file_ext != ".TRC":
            # Create events matrix
            event_matrix = np.ndarray((self.nb_events,3),dtype = int)
            event_matrix[:,0] = (self.events_array)
            event_matrix[:,1] = 0
            event_matrix[:,2] = self.event_index
            # Add events to data
            mne_raw.add_events(event_matrix)    
        
        print("--")
        print("raw data shape: ", mne.find_events(mne_raw).shape)
        print("..")

    def plot_data(self,mne_raw, ev_dict, duration, n_channels, lowpass = None, highpass = None, order = None):
        """
        Plot raw data using mne GUI.
        """
        
        self.events_forplot = mne.find_events(mne_raw)
                
        #events = mne.read_events(op.join(data_path, 'sample_audvis_raw-eve.fif'))
        scalings = {'mag': 2, 'grad': 2}
        mne_raw.plot(scalings=scalings, title='Data from arrays', show=True, 
                          block=True, events = self.events_forplot, event_id = ev_dict,
                          duration = duration, n_channels = n_channels, lowpass = lowpass,
                          highpass = highpass, order = order)
        
if __name__ == "__main__":    
        
    # Change it
    #-----------     

    filepath = r"F:\GardyL\Data_storage\Imagery\IMAGERY_AB28\MACRO\Bipolaire\Imagery1.edf"
    filepath = r"F:\GardyL\Data_storage\Imagery\IMAGERY_AB28\MACRO\Bipolaire\Imagery234.edf"
    filepath = r"F:\GardyL\Scenario RM02\CHU Service Epilepsie\SEEG\22. 9h30 CRISE 2 (bleu)\EEG_114.TRC"

    filepath = r"F:\GardyL\Data_storage\OddballAndSabData\029_FF_SAB_OddBall\FF29-sab1-2.edf"

    sfreq = 256
    choose_event = []
    file_day =  []
    file_part =  []
    
    task = "SAB" 

    # Events path (not needed if .trc file <<events already are loaded>>)
    #----------------------------------------------------------------------
    events_path = r"F:\GardyL\Data_storage\Imagery\IMAGERY_AB28\Triggers1.csv"
    events_path = r"F:\GardyL\Data_storage\Imagery\IMAGERY_AB28\Triggers2.csv"

    events_path = r"F:\GardyL\Data_storage\OddballAndSabData\029_FF_SAB_OddBall\comportement\events\Triggers1.xlsx"

    # "event_name": event_position
    ev_dict = {"start": 1, "end" : 126}
    
    # Run
    #-----
    mne_raw = io_eeg_to_mne(filepath = filepath, rescale_voltage = True, downsampling_freq = 256)[0]
    mne_raw._data.shape
    
    add_events = add_events_to_mne_sig(mne_raw, filepath)
    add_events.import_events(mne_raw, events_path)
    add_events.eventsTime_to_eventsTimestamps(mne_raw)
    if task != "SAB":
        add_events.delete_false_triggers(mne_raw)
    add_events.microMacro_realign_eventstimes(mne_raw, plot_triggers=True)
    add_events.eventstime_periods_toKeep(mne_raw)
    add_events.remove_bad_events()
    add_events.cleanArray_to_eventsArray()
    add_events.add_events_to_data(mne_raw)    

    scalings = {'mag': 2, 'grad': 2}
        
    mne_raw.plot(scalings=scalings, title='Data from arrays', show=True, 
                      block=True, events = mne.find_events(mne_raw), event_id = None,
                      duration = 10, n_channels = 40, lowpass = None,
                      highpass = None, order = None)     
    
    # Print events times (and not time stamps)
    events_time = mne.find_events(mne_raw)[:,0] / sfreq
    print(events_time)
