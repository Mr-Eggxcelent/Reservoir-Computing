# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

#Refer to Naveen Venkatesan article on publication quality plotting
##https://towardsdatascience.com/an-introduction-to-making-scientific-publication-plots-with-python-ea19dfa7f51e


import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
import numpy as np
from pylab import cm
from sklearn.metrics import mean_squared_error

mpl.rc('legend',**{'fontsize':12})
plt.rcParams['font.size'] = 18
plt.rcParams['axes.linewidth'] = 2


# Load data
signal_1 = pd.read_csv('Results/Record_3NF/outputsignal.csv',header=None)
signal_2 = pd.read_csv('Results/Record_3NF/outputsignal_two.csv',header=None)
signal_3 = pd.read_csv('Results/Record_3NF/outputsignal_three.csv',header=None)

target_1=pd.read_csv('../Data/volterra.csv',header=None)
target_2=pd.read_csv('../Data/2ndOrder.csv',header=None)
target_3=pd.read_csv('../Data/NARMA.csv',header=None)
input_signal=pd.read_csv('../Data/input.csv',header=None)

data_x = np.linspace(0, 5,signal_1.shape[0])
target_x = np.linspace(0, 5,5000)          
                                        

def plot_signal(signal,target,x_limit_min,x_limit_max,y_limit_min,y_limit_max):
    
    fig = plt.figure(figsize=(5, 3))
    axs = fig.add_axes([0, 0, 1, 1])
    
    axs.xaxis.set_tick_params(which='major', size=10, width=2, direction='in', top='on')
    axs.xaxis.set_tick_params(which='minor', size=7, width=2, direction='in', top='on')
    axs.yaxis.set_tick_params(which='major', size=10, width=2, direction='in', right='on')
    axs.yaxis.set_tick_params(which='minor', size=7, width=2, direction='in', right='on')
      
    axs.plot(target_x, target[130000:135000],'-b',label='target')
    axs.plot(data_x, signal,'--r',label='MC output')
     
    axs.set_xlim(x_limit_min, x_limit_max)
    axs.set_ylim(y_limit_min, y_limit_max)
    axs.set_xlabel('time(s)')
    axs.legend(loc="lower right");
    plt.show()
      
    signal_MSE=mean_squared_error(target[0][130000:135000], signal[0])
    return signal_MSE


def main():
    MSE=plot_signal(signal_2,target_2,0,5,-5,5)
    print("Signal MSE:%.8f"%(MSE))
    
    
if __name__ == "__main__":
    main()


