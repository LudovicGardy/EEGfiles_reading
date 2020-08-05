"""
Original file by Deudon M.,
modified by GardyL. (19/12/2018)
"""

import numpy as np
import re
import mne
import neo.io
import matplotlib.pyplot as plt
import pandas as pd

def io_eeg_to_mne(filepath, rescale_voltage = False, read_data = True, export_data = False, downsampling_freq = [], events_array = [], load_events = False):
    """ Import an EEG file to MNE as a Raw instance.

    Supported formats : 'edf', 'eeg', 'trc'

    Parameters
    ----------
    filepath : str
        EEG filepath

    read_data : bool (default: True)
        If True, read the data

    Returns
    -------
    mne_raw : MNE RAW instance
        Output MNE structure

    ch_names : list
        Channel names

    Note
    -----
    Events in mne_raw event channel must be in timestamps (not in seconds)
    """

    def init():
        """
        Initialization
        """
        possible_ext = ['.edf', '.fif', '.trc', '.ns5', '.nsx', '.bdf']
        file_ext = re.search('\.\w+$', filepath)
        if file_ext:
            file_ext = file_ext[0].lower()
        else:
            raise ValueError('Could not detect file extension of file {}'.format(filepath))
        if file_ext not in possible_ext:
            raise ValueError('The file {} has not a supported extension. Extensions must be in {}'.format(filepath,
                                                                                                          possible_ext))

        return(file_ext)

    def add_events_to_mneraw(mne_raw, events_array = []):
        """
        If file_ext != TRC, add events to the object mne_raw.

        Parameters
        ----------
        mne_raw : MNE RAW instance
            Output MNE structure

        events_array : list
            events list in seconds will be transformed into timestamps

        Returns
        -------
        mne.find_events(mne_raw) : 3D matrix
            events timestamps (col 0) + index (col 2)
        """
        events_array_timestamps = events_array * mne_raw.info["sfreq"]
        events_array = events_array_timestamps

        nb_events = len(events_array)
        event_index = np.arange(1,(nb_events + 1))

        # Create events matrix
        event_matrix = np.ndarray((nb_events,3),dtype = int)
        event_matrix[:,0] = (events_array)
        event_matrix[:,1] = 0
        event_matrix[:,2] = event_index
        # Add events to data
        mne_raw.add_events(event_matrix)

        print("--")
        print("raw data shape: ", mne.find_events(mne_raw).shape)
        print("..")

    def load_data_from_file(file_ext, filepath, events_array = [], stim_chan = None):
        """
        Load data from file.

        Parameters
        ----------
        file_ext : str
            file extension is necessary for chosing the loading function

        filepath : str
            data file (possible formats are: edf, trc, ns5, fif)

        events_array : list (optional)
            for plotting vertial lines on the EEG graph when events occur

        Returns
        -------
        mne_raw : MNE RAW instance
            Output MNE structure

        ch_names : list
            Channel names

        events_dict : dict (optional)
            Only if .TRC file
        """

        mne_raw = []
        events_dict = {}

        # From edf
        if file_ext in ['.edf']:
            try:
                mne_raw = mne.io.read_raw_edf(filepath, preload=read_data, stim_channel = stim_chan)
            except:
                mne_raw = mne.io.read_raw_edf(filepath, preload=True, stim_channel = stim_chan)
            ch_names = mne_raw.ch_names

            # rescale voltage *1e6 (for plotting)
            if rescale_voltage == True:
                not_eeg_list = ["STI 014", "MKR2+"]

                eeg_bool = []
                for i in ch_names:
                    if i in not_eeg_list:
                        eeg_bool.append(True)
                    else:
                        eeg_bool.append(False)

                mne_raw._data[eeg_bool] = mne_raw._data[eeg_bool] * 1e6

            if load_events == True:
                add_events_to_mneraw(mne_raw, events_array)

        # From bdf
        elif file_ext in ['.bdf']:
            try:
                mne_raw = mne.io.read_raw_bdf(filepath, preload=read_data, stim_channel = stim_chan)
            except:
                mne_raw = mne.io.read_raw_bdf(filepath, preload=True, stim_channel = stim_chan)
            ch_names = mne_raw.ch_names

            # rescale voltage *1e6 (for plotting)
            if rescale_voltage == True:
                not_eeg_list = ["STI 014", "MKR2+"]

                eeg_bool = []
                for i in ch_names:
                    if i in not_eeg_list:
                        eeg_bool.append(True)
                    else:
                        eeg_bool.append(False)

                mne_raw._data[eeg_bool] = mne_raw._data[eeg_bool] * 1e6

            if load_events == True:
                add_events_to_mneraw(mne_raw, events_array)

        # From fif
        elif file_ext == '.fif':
            try:
                mne_raw = mne.io.read_raw_fif(filepath, preload=read_data)
            except:
                mne_raw = mne.io.read_raw_fif(filepath, preload=True)
            ch_names = mne_raw.ch_names

            if load_events == True:
                add_events_to_mneraw(mne_raw, events_array)

        # From trc
        elif file_ext == '.trc':
            trc_reader = neo.io.MicromedIO(filename=filepath)
            header = trc_reader.header
            ch_names = [header['signal_channels'][i][0] for i in range(trc_reader.signal_channels_count())]
            ch_names.append('STI 014')
            ch_types =['eeg' for _ in ch_names[0:-1]] + ['stim']

            if read_data:
                bl = trc_reader.read(lazy=False)[0]
                seg = bl.segments[0]
                len_seg = len(seg.analogsignals)
                sfreq = int(seg.analogsignals[0].sampling_rate.magnitude)

                # add stim channel (events)
                if len(trc_reader._raw_events) > 0:
                    raw_events = trc_reader._raw_events
                    events_time_list = []
                    events_label_list = []

                    for time, label in raw_events[1]:
                        events_time_list.append(time)
                        events_label_list.append(label)

                    n_events = len(events_time_list)

                    events = np.empty([n_events, 3])
                    events[:, 0] = events_time_list
                    events[:, 2] = range(1, (n_events + 1))

                    events_dict = {"time": events_time_list, "label": events_label_list, "n_event" : range(0,n_events)}

                n_chan = []
                for i in range(0,len_seg):
                    n_chan.append(len(seg.analogsignals[i][0]))
                n_chan = sum(n_chan)

                # add stim channel (events)
                if "stim" in ch_types:
                    n_chan += 1

                data_list = []
                for i in range(0, len_seg):
                    seg_new = seg.analogsignals[i].T
                    n_pnts, n_chan = len(seg_new[0]), len(seg_new)
                    data = np.zeros((n_chan, n_pnts), dtype=float)
                    for i, asig in enumerate(seg_new):
                        # We need the ravel() here because Neo < 0.5 gave 1D, Neo 0.5 gives
                        # 2D (but still a single channel).
                        data[i] = asig.magnitude.ravel()
                    data_list.append(data)

                for i in range(0,len_seg):
                    if i == 0:
                        data_full = data_list[0]
                    else:
                        data_full = np.concatenate((data_full, data_list[i]), axis = 0)

                # add stim channel (events)
                data_full = np.vstack((data_full, np.zeros((1, data_full.shape[1]))))

                info = mne.create_info(ch_names=ch_names, ch_types=ch_types, sfreq=sfreq)
                mne_raw = mne.io.RawArray(data_full, info)
                mne_raw.add_events(events)

        # From ns5, nsx (new version : 15/01/2020)
        elif file_ext == '.ns5':
            nsx_reader = neo.io.BlackrockIO(filename=filepath)
            header = nsx_reader.header
            #print(header)
            ch_names = [header['signal_channels'][i][0] for i in range(nsx_reader.signal_channels_count())]
            ch_types =['seeg' for _ in ch_names]

            if read_data:
                bl = nsx_reader.read(lazy=False)[0]
                seg = bl.segments[0]
                len_seg = len(seg.analogsignals)
                #sfreq = int(seg.analogsignals[0].sampling_rate.magnitude)
                sfreq = 30000

                for i in range(0,len_seg):
                    if i == 0:
                        print("First channel (len = {}) removed.".format(len(seg.analogsignals[i]))) # /*****/ should I always do it ???????? or only when header['nb_segments'] > 1
                    else:
                        print("Channel {} (len = {}) added to data.".format(i, len(seg.analogsignals[i])))

                n_chan = []
                #for i in range(0,len_seg): # /*****/ ???????????????
                for i in range(1,len_seg):
                    n_chan.append(len(seg.analogsignals[i][0]))
                n_chan = sum(n_chan)

                data = []
                #for i in range(0,len_seg): # /*****/ ???????????????
                for i in range(1, len_seg):
                    seg_new = seg.analogsignals[i].T
                    n_pnts = len(seg_new[0])
                    #data = np.zeros((n_chan, n_pnts), dtype=float)
                    data.append(seg_new[0].magnitude.ravel())
                data = np.array(data)

                info = mne.create_info(ch_names=ch_names, ch_types=ch_types, sfreq=sfreq)
                mne_raw = mne.io.RawArray(data, info)

        ## From ns5 or nsx (old version)
        #elif file_ext in ['.ns5', 'nsx']:
        #    nsx_reader = neo.io.BlackrockIO(filename=filepath)
        #    header = nsx_reader.header
        #    ch_names = [header['signal_channels'][i][0] for i in range(nsx_reader.signal_channels_count())]
        #    if read_data:
        #        bl = nsx_reader.read(lazy=False)[0]
        #        seg = bl.segments[0]
        #        n_pnts, n_chan = len(seg.analogsignals[0]), len(seg.analogsignals)
        #        data = np.zeros((n_chan, n_pnts), dtype=float)
        #        for i, asig in enumerate(seg.analogsignals):
        #            # We need the ravel() here because Neo < 0.5 gave 1D, Neo 0.5 gives
        #            # 2D (but still a single channel).
        #            data[i, :] = asig.magnitude.ravel()
        #        sfreq = int(seg.analogsignals[0].sampling_rate.magnitude)
        #        info = mne.create_info(ch_names=ch_names, sfreq=sfreq)
        #        mne_raw = mne.io.RawArray(data, info)

            if load_events == True:
                add_events_to_mneraw(mne_raw, events_array)

        return(mne_raw, ch_names, events_dict)

    def rescale_voltage(mne_raw, rescale_voltage):
        """
        Rescale voltage *1e-6 (from volt to microvolt. for plotting.)

        Parameters
        ----------
        mne_raw : MNE RAW instance
            Output MNE structure

        rescale_voltage : bool
            Multiply voltage by 1e-6 if TRUE
        """

        if rescale_voltage == True:
            not_eeg_list = ["STI 014", "MKR2+"]

            eeg_bool = []
            for i in ch_names:
                if i in not_eeg_list:
                    eeg_bool.append(True)
                else:
                    eeg_bool.append(False)

            mne_raw._data[eeg_bool] = mne_raw._data[eeg_bool] * 1e-6

    def export_to_fif(mne_raw, export_data):
        """
        Export mne data to fif format
        """

        # If want to export data (to .fif)
        if export_data == True:
            mne_raw.save(export_path, overwrite= True)

        if downsampling_freq:
            mne_raw.resample(sfreq = downsampling_freq)


    file_ext = init()
    mne_raw, ch_names, events_dict = load_data_from_file(file_ext, filepath, events_array= events_array)
    rescale_voltage(mne_raw, rescale_voltage)
    export_to_fif(mne_raw, export_data)

    return(mne_raw, ch_names, events_dict)

if __name__ == "__main__":

    #-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_#
    #      You should change these       #
    #-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_#

    # Data path parameters
    sub_num = 5
    sess_num = 3
    run_num = 3
    chantype = "micro"

    #-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_#
    #     You might not change these     #
    #-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_#

    # Filter parameters
    apply_filter = False
    filter_boundaries = [200, 600]

    # Open plot signal GUI
    open_signalVis_GUI = True

    # Signal parameters
    export_data = False
    rescale_voltage = False
    load_events = False
    plot_mne = True
    events_array = []

    #-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_#
    #         Do not change these        #
    #-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_#

    # Set file path according to channel type (micro/macro)
    if chantype == "macro":
        filepath = r"F:\GardyL\Data_storage\EPIFAR_storage\BIDS_data\macro_bipolar_BIDS\sub-{:03}\ses-{:02}\eeg\sub-{:03}_ses-{:02}_task-EPIFAR_run-{:02}_eeg.edf".format(sub_num, sess_num, sub_num, sess_num, run_num)
    elif chantype == "micro":
        filepath = r"F:\GardyL\Data_storage\EPIFAR_storage\BIDS_data\micro_5KHz_BIDS\sub-{:03}\ses-{:02}\eeg\sub-{:03}_ses-{:02}_task-EPIFAR_run-{:02}_eeg.edf".format(sub_num, sess_num, sub_num, sess_num, run_num)

    # Load data
    mne_struct = io_eeg_to_mne(filepath = filepath, rescale_voltage = rescale_voltage, load_events = load_events) # downsampling_freq = 2000
    mne_raw = mne_struct[0]
    mne_raw._data.shape
    mne_raw.info["ch_names"]

    # Get some info from the data
    events_dict = mne_struct[2]
    sfreq = mne_raw.info["sfreq"]
    time_stamps = len(mne_raw.get_data()[0])
    recording_duration_seconds = round(time_stamps / sfreq, 2)
    recording_duration_minutes = round(recording_duration_seconds / 60, 2)

    print("\n\nrecording duration = {} seconds ({} minutes)".format(recording_duration_seconds, recording_duration_minutes))

    # Bandpass filter on signal
    if apply_filter:
        mne_raw.filter(l_freq = filter_boundaries[0], h_freq = filter_boundaries[1])

    # Plot data
    if open_signalVis_GUI:
        #events = mne.find_events(mne_raw, stim_channel=None)
        scalings = 1e-3 # {'mag': 2, 'grad': 2}
        ev_dict = {"test_1": 1, "test_2":2, "test_3":3, "test_4":4, "test_5":5, "test_6":6 }
        mne_raw.plot(scalings=scalings, title='Data from arrays', show=True, block=True,
                        event_id= ev_dict, duration = 0.6, n_channels = 12)
