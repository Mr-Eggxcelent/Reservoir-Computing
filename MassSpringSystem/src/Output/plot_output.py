# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


# Load data
signal_1 = pd.read_csv('Results/outputsignal_washout.csv',header=None)
target_signal_1 = pd.read_csv('Results/targetsignal_washout.csv',header=None)

signal_2 = pd.read_csv('Results/outputsignal_two_washout.csv',header=None)
target_signal_2 = pd.read_csv('Results/targetsignal_two_washout.csv',header=None)

signal_3 = pd.read_csv('Results/outputsignal_three_washout.csv',header=None)
target_signal_3 = pd.read_csv('Results/targetsignal_three_washout.csv',header=None)

merged = pd.read_csv('merged_test_output.csv',header=None)
feedback= pd.read_csv('merged_feedback.csv',header=None)

data_x = np.linspace(0, 1,20000)
merged_x = np.linspace(0, 1,120000)
target_x = np.linspace(0, 1,20000)          
                                        

fig, axs = plt.subplots(5)
fig.suptitle('Vertically stacked subplots')

axs[0].plot(merged_x, merged[0],'-r',label='Test Output')
axs[0].plot(merged_x, feedback[0],'--r',label='Closed Loop Feedback')
axs[0].plot(merged_x, merged[1],'-g',label='Test Output Second Eq')
axs[0].plot(merged_x, feedback[1],'--y',label='Closed Loop Feedback Second Eq')
leg=axs[0].legend();

axs[1].plot(data_x, signal_1[0],'-y',label='x_1')
axs[1].plot(data_x, signal_1[1],'--y',label='x_2')
##axs[1].plot(target_x, target_signal_1[0],'--r',label='target eq 1')
##axs[1].plot(target_x, target_signal_1[1],'--y',label='target eq 2')
leg=axs[1].legend();

axs[2].plot(data_x, signal_2[0],'-g',label='x_1')
axs[2].plot(data_x, signal_2[1],'--g',label='x_2')
##axs[2].plot(target_x, target_signal_2[0],'--r',label='target eq 1')
##axs[2].plot(target_x, target_signal_2[1],'--y',label='target eq 2')
leg=axs[2].legend();

axs[3].plot(data_x, signal_3[0],'-b',label='x_1')
axs[3].plot(data_x, signal_3[1],'--b',label='x_2')
##axs[3].plot(target_x, target_signal_3[0],'--r',label='target eq 1')
##axs[3].plot(target_x, target_signal_3[1],'--y',label='target eq 2')
leg=axs[3].legend();


axs[4].plot(signal_1[0], signal_1[1],'-y',label='Van der Pol')
axs[4].plot(signal_2[0], signal_2[1],'-g',label='Quad')
axs[4].plot(signal_3[0], signal_3[1],'-b',label='Optional Signal')
leg=axs[4].legend();

plt.show()

