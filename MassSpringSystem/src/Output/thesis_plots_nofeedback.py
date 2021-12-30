# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from sklearn.metrics import mean_squared_error


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
                                        

fig, axs = plt.subplots(1)
fig.suptitle('Volterra')
axs.plot(data_x, signal_1,'--r',label='MC output')
axs.set_xlabel('time(s)')
axs.plot(target_x, target_1[130000:135000],'-b',label='target')
leg=axs.legend(loc="lower right");

signal1_MSE=mean_squared_error(target_1[0][130000:135000], signal_1[0])


fig, axs = plt.subplots(1)
fig.suptitle('2ndOrder')
axs.plot(data_x, signal_2,'--r',label='MC output')
axs.set_xlabel('time(s)')
axs.plot(target_x, target_2[130000:135000],'-b',label='target')
leg=axs.legend(loc="lower right");

signal2_MSE=mean_squared_error( target_2[0][130000:135000], signal_2[0])


fig, axs = plt.subplots(1)
fig.suptitle('NARMA')
axs.plot(data_x, signal_3,'--r',label='MC output')
axs.plot(target_x, target_3[130000:135000],'-b',label='target')
axs.set_xlabel('time(s)')
leg=axs.legend(loc="lower right");

signal3_MSE=mean_squared_error(target_3[0][130000:135000], signal_3[0])


fig, axs = plt.subplots(1)
axs.plot(data_x, input_signal[130000:135000],'-b')
axs.set_xlabel('time(s)')
##leg=axs.legend(loc="lower right");


print("Volterra:%.8f \n 2ndOrder:%.8f \n NARMA:%.8f " %(signal1_MSE,signal2_MSE,signal3_MSE))


plt.show()
