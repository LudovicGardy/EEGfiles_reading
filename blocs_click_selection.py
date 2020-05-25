# -*- coding: utf-8 -*-
"""
Created on Mon Feb 18 16:22:29 2019

@author: gardy
"""

import numpy as np
import matplotlib.pyplot as plt

def clickable_graph(x_coords, y_coords):
    
    def onclick(event):
        print('%s click: button=%d, x=%d, y=%d, xdata=%f, ydata=%f' %
              ('double' if event.dblclick else 'single', event.button,
               event.x, event.y, event.xdata, event.ydata))
        coords.append((event.x, event.xdata))
      
        ax.vlines(x = event.xdata, ymin = min(y_coords), ymax = max(y_coords), linestyle = "--", color = "magenta")
        
        count.append(event.xdata)
        if len(count) % 2 != 0:
            ax.text(event.xdata,event.ydata,"begin \n({})".format(round(event.xdata,2)))
        else:
            ax.text(event.xdata,event.ydata,"end \n({})".format(round(event.xdata,2)))

        fig.canvas.draw()
        fig.canvas.flush_events()
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(x_coords, y_coords)
    
    count = []
    
    global coords    
    coords = []

    cid = fig.canvas.mpl_connect('button_press_event', onclick)
    plt.show()
    
if __name__ == "__main__":
    x = np.arange(-10,10)
    y = x ** 2
    clickable_graph(x, y)
    
    
    
