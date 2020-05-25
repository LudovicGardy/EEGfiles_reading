# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 12:03:08 2019

@author: GARDY
"""
# general imports
import numpy as np
import matplotlib.pyplot as plt
from dtw import dtw
from fastdtw import fastdtw
from scipy.spatial.distance import euclidean

# from computer
from l2_norm import l2_norm

class pattern_recognition:
    def __init__(self):
        self.color = "black"
    
    def exact_pattern_finder(self, signal, pattern_input):
        """
        Calculate the distance between a query (the pattern) and
        a part of the signal (same length)
        
        Parameters
        ----------
        signal : list or numpy array
            full signal
        
        pattern_input : list or numpy array
            query
            
            
        Returns
        --------
        self.distances_results
            to comment
            
        self.mindist
            to comment


        self.where_mindist
            to comment
          
        """       
        
        begin = 0
        pattern_length = len(pattern_input)
        end = len(signal) - pattern_length

        self.color = "orange"       
        self.distances_result = []
        
        d = []
        
        for i in range(begin,end):
            a = i
            b = i + pattern_length
            
            sig = signal[a:b]
            print(a, " to ",b, ". (len = {})".format(end))
                    
            d_tmp = abs(np.array(sig) - np.array(pattern_input))
                        
            d.append(d_tmp)
            
        for i in d:    
            self.distances_result.append(sum(i))

        self.distances_result = np.array(self.distances_result)
        self.mindist = np.min(self.distances_result)
        #self.where_mindist = np.array(np.where(self.distances_result == np.min(self.distances_result))[0])
        self.where_mindist = np.array(np.where(self.distances_result == 0)[0])
            
    def real_triggers_begin_end(self,pattern_input, begin_or_end):
        
        self.begin_or_end = begin_or_end
        
        if self.begin_or_end == "begin":
            nb_999_in_pattern = len(np.where(pattern_input == 999)[0])
            self.where_mindist =  self.where_mindist + nb_999_in_pattern + 1

        if self.begin_or_end == "end":
            nb_999_in_pattern = len(np.where(pattern_input == 999)[0])
            self.where_mindist =  self.where_mindist + len(pattern_input) - nb_999_in_pattern 
        
    def plot_distance_result(self, pattern_input):
        plt.figure()
        plt.plot(self.distances_result, color = self.color)
        plt.show()

        print("")
        print("--------------------")
        print("min dist = {}".format(self.mindist))
        print("")
        print("where min dist = {}".format(self.where_mindist))
        print("--------------------")
        print("")        

if __name__=="__main__":
    
    pattern = {
    "begin" : np.array([999,999,999,0,0,0,0,0,0,0,0,0,0]),
    "end" : np.array([0,0,0,0,0,0,0,0,0,0, 999,999,999])  
    }
    
    # signal & pattern
    #------------------
    sig_full = imagery_eeg.check_trigs_BeginningEnding
    begin_or_end = "end"
    
    pattern = pattern[begin_or_end]

    f, ax = plt.subplots()
    ax.plot(sig_full)
    ax.set_title("noise + hidden pattern")
    plt.show()
    
    # Create object
    p_recog = pattern_recognition()
    
    # dtw method
    p_recog.exact_pattern_finder(sig_full, pattern)
    p_recog.real_triggers_begin_end(pattern,begin_or_end)
    p_recog.plot_distance_result(pattern)
