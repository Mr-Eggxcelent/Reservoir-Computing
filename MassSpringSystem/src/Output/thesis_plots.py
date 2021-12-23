# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from sklearn.metrics import mean_squared_error
from matplotlib.collections import LineCollection

# Load data
signal_1 = pd.read_csv('Results/Record_1/outputsignal_washout.csv',header=None)
target_f_1 = pd.read_csv('Results/Record_1/targetsignal_washout.csv',header=None)

signal_2 = pd.read_csv('Results/Record_1/outputsignal_two_washout.csv',header=None)
target_f_2 = pd.read_csv('Results/Record_1/targetsignal_two_washout.csv',header=None)

signal_3 = pd.read_csv('Results/Record_1/outputsignal_three_washout.csv',header=None)
target_f_3 = pd.read_csv('Results/Record_1/targetsignal_three_washout.csv',header=None)

merged = pd.read_csv('Results/Record_1/merged_test_output.csv',header=None)
merged_openloop = pd.read_csv('Results/Record_1/merged_feedback.csv',header=None)

target_1=pd.read_csv('../Data/vanderpol.csv',header=None)
target_2=pd.read_csv('../Data/quad.csv',header=None)
target_3=pd.read_csv('../Data/lissajous.csv',header=None)

data_x = np.linspace(0, 15,15000)
merged_x = np.linspace(0, 150,150000)
merged_x_openloop = np.linspace(0, 150,150000)
target_x = np.linspace(0, 15,15000)          

data_test = np.linspace(0, 20,20000)
                                        
##fig, axs = plt.subplots(1)
##fig.suptitle('Van der Pol')

##axs.plot(data_x, signal_1[0][5000:20000],'--b',label='$x_1$')
##axs.plot(data_x, signal_1[1][5000:20000],'--g',label='$x_2$')
##axs.plot(data_x, target_1[0][165000:180000],'b',label='$x_1$ target')
##axs.plot(data_x, target_1[1][165000:180000],'g',label='$x_2$ target')
##axs.set_xlabel('time(s)')
##axs.plot(signal_1[0][5000:20000], signal_1[1][5000:20000],'--r',label='MC output')
##axs.plot(target_1[0][165000:180000],target_1[1][165000:180000],'b',label='target')
##axs.set_xlabel('$x_1$')
##axs.set_ylabel('$x_2$')
##leg=axs.legend(loc="lower right");

signal1_MSE=mean_squared_error(target_1[0][165000:180000], signal_1[0][5000:20000])
signal1_col2_MSE=mean_squared_error(target_1[1][165000:180000], signal_1[1][5000:20000])


##fig, axs = plt.subplots(1)
##fig.suptitle('Quadratic')
##axs.plot(data_x, signal_2[0][5000:20000],'--b',label='$x_1$')
##axs.plot(data_x, signal_2[1][5000:20000],'--g',label='$x_2$')
##axs.plot(data_x, target_2[0][165000:180000],'b',label='$x_1$ target')
##axs.plot(data_x, target_2[1][165000:180000],'g',label='$x_2$ target')
##axs.set_xlabel('time(s)')
##axs.plot(signal_2[0][5000:20000], signal_2[1][5000:20000],'--r',label='MC output')
##axs.plot(target_2[0][165000:180000],target_2[1][165000:180000],'b',label='target')
##axs.set_xlabel('$x_1$')
##axs.set_ylabel('$x_2$')
##leg=axs.legend(loc="lower right");


signal2_MSE=mean_squared_error( target_2[0][165000:180000], signal_2[0][5000:20000])
signal2_col2_MSE=mean_squared_error(target_2[1][165000:180000], signal_2[1][5000:20000])


##fig, axs = plt.subplots(1)
##fig.suptitle('Lissajous Figure')
##axs.plot(data_x, signal_3[0][5000:20000],'--b',label='$x_1$')
##axs.plot(data_x, signal_3[1][5000:20000],'--g',label='$x_2$')
##axs.plot(data_x, target_3[0][165000:180000],'b',label='$x_1$ target')
##axs.plot(data_x, target_3[1][165000:180000],'g',label='$x_2$ target')
##axs.set_xlabel('time(s)')
##axs.plot(signal_3[0][5000:20000], signal_3[1][5000:20000],'--r',label='MC output')
##axs.plot(target_3[0][165000:180000],target_3[1][165000:180000],'b',label='target')
##axs.set_xlabel('$x_1$')
##axs.set_ylabel('$x_2$')
##leg=axs.legend(loc="lower right");

signal3_MSE=mean_squared_error(target_3[0][165000:180000], signal_3[0][5000:20000])
signal3_col2_MSE=mean_squared_error(target_3[1][165000:180000], signal_3[1][5000:20000])


fig, axs = plt.subplots(1)
axs.plot(merged_x, merged[0],'b',label='$x_1$')
axs.plot(merged_x, merged[1],'g',label='$x_2$')
axs.set_xlabel('time(s)')
leg=axs.legend(loc="lower right");


##fig, axs = plt.subplots(1)
##axs.plot(data_test, target_2[0][160000:180000],'b',label='$x_1$ target')
##axs.plot(data_test, target_2[1][160000:180000],'g',label='$x_2$ target')
##points_s2_c1 = np.array([data_test, signal_2[0]]).T.reshape(-1,1,2)
##points_s2_c2 = np.array([data_test, signal_2[1]]).T.reshape(-1,1,2)
##
##
##segments_s2_c1 = np.concatenate([points_s2_c1[:-1],points_s2_c1[1:]], axis=1)
##segments_s2_c2 = np.concatenate([points_s2_c2[:-1],points_s2_c2[1:]], axis=1)
##
##
##lc_3 = LineCollection(segments_s2_c1, cmap="viridis", linewidth=3)
##lc_4 = LineCollection(segments_s2_c2, cmap="viridis", linewidth=3)
##
##
##lc_3.set_array(data_test)
##lc_4.set_array(data_test)
##
##axs.add_collection(lc_3)
##axs.add_collection(lc_4)
##
##line = axs.add_collection(lc_3)
##leg=axs.legend(loc="lower right");
##
##fig.colorbar(line, ax=axs,orientation='horizontal')
##axs.autoscale_view()
##axs.set_xlabel('time(s)')


##fig, axs = plt.subplots(1)
##axs.plot(target_2[0][160000:180000], target_2[1][160000:180000],'b',label='$x_1$ target')
##points_cycle = np.array([signal_2[0], signal_2[1]]).T.reshape(-1,1,2)
##segments_cycle = np.concatenate([points_cycle[:-1],points_cycle[1:]], axis=1)
##lc_5 = LineCollection(segments_cycle, cmap="viridis", linewidth=3)
##lc_5.set_array(data_test)
##axs.set_xlabel('$x_1$')
##axs.set_ylabel('$x_2$')
##axs.add_collection(lc_5)
##line = axs.add_collection(lc_5)
##fig.colorbar(line, ax=axs, orientation='horizontal')
##leg=axs.legend(loc="lower right");
##axs.autoscale_view()


print("Van der Pol:%.8f \t %.8f \n Quadratic:%.8f \t %.8f \n Lissajous:%.8f \t %.8f " %(signal1_MSE,signal1_col2_MSE,signal2_MSE,signal2_col2_MSE,signal3_MSE,signal3_col2_MSE))

plt.show()
