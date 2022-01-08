import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from sklearn.metrics import mean_squared_error
from matplotlib.collections import LineCollection
from matplotlib.ticker import MaxNLocator
from matplotlib import colors

plt.rcParams['font.size'] = 18
plt.rcParams['axes.linewidth'] = 2

# Load data
signal_1 = pd.read_csv('Results/quad_limit_cycles/outputsignal_washout.csv',header=None)
target_f_1 = pd.read_csv('Results/quad_limit_cycles/targetsignal_washout.csv',header=None)

signal_2 = pd.read_csv('Results/quad_limit_cycles/outputsignal_two_washout.csv',header=None)
target_f_2 = pd.read_csv('Results/quad_limit_cycles/targetsignal_two_washout.csv',header=None)

signal_3 = pd.read_csv('Results/quad_limit_cycles/outputsignal_three_washout.csv',header=None)
target_f_3 = pd.read_csv('Results/quad_limit_cycles/targetsignal_three_washout.csv',header=None)

feedback_signal = pd.read_csv('Results/quad_limit_cycles/targetsignal_washout.csv',header=None)

merged = pd.read_csv('Results/quad_limit_cycles/merged_test_output.csv',header=None)
merged_openloop = pd.read_csv('Results/quad_limit_cycles/merged_feedback.csv',header=None)

target_1=pd.read_csv('../Data/quad_1.csv',header=None)
target_2=pd.read_csv('../Data/quad_2.csv',header=None)
target_3=pd.read_csv('../Data/quad_3.csv',header=None)

#Change values of according to size of arrays
#Or use len(df.index) or df.shape[0] where df is the dataframe of choice
merged_x = np.linspace(0, 180,180000)   
data_test = np.linspace(0, 30,30000)
                                        


def timed_limit_cycle(signal,target,x_limit_min,x_limit_max,y_limit_min,y_limit_max,data_points):
    
    fig = plt.figure(figsize=(6.4, 4.8))
    axs = fig.add_axes([0, 0, 1, 1])
    axs.xaxis.set_tick_params(which='major', size=10, width=2, direction='in', top='on')
    axs.xaxis.set_tick_params(which='minor', size=7, width=2, direction='in', top='on')
    axs.yaxis.set_tick_params(which='major', size=10, width=2, direction='in', right='on')
    axs.yaxis.set_tick_params(which='minor', size=7, width=2, direction='in', right='on')
    axs.set_xlim(x_limit_min,x_limit_max)
    axs.set_ylim(y_limit_min,y_limit_max)
    
    # cmap = colors.ListedColormap(['green', 'orange','navy'])
    axs.plot(target[0][130000:180000], target[1][130000:180000],'--r',label='target',linewidth=1)
    points_cycle = np.array([signal[0], signal[1]]).T.reshape(-1,1,2)
    segments_cycle = np.concatenate([points_cycle[:-1],points_cycle[1:]], axis=1)
    lc_5 = LineCollection(segments_cycle, cmap="viridis", linewidth=3)
    lc_5.set_array(data_points)
    
    axs.set_xlabel('$x_1$')
    axs.set_ylabel('$x_2$')
    axs.add_collection(lc_5)
    line = axs.add_collection(lc_5)
    
    fig.colorbar(line, ax=axs, orientation='horizontal')
    axs.legend(loc="lower right");
    #axs.autoscale_view()
    
def timed_signal(signal,target,x_limit_min,x_limit_max,y_limit_min,y_limit_max,data_points):
    
    fig = plt.figure(figsize=(6.4, 4.8))
    axs = fig.add_axes([0, 0, 1, 1])
    axs.xaxis.set_tick_params(which='major', size=10, width=2, direction='in', top='on')
    axs.xaxis.set_tick_params(which='minor', size=7, width=2, direction='in', top='on')
    axs.yaxis.set_tick_params(which='major', size=10, width=2, direction='in', right='on')
    axs.yaxis.set_tick_params(which='minor', size=7, width=2, direction='in', right='on')
    axs.xaxis.set_major_locator(MaxNLocator(integer=True))
    axs.yaxis.set_major_locator(MaxNLocator(integer=True))
    axs.set_xlim(x_limit_min,x_limit_max)
    axs.set_ylim(y_limit_min,y_limit_max)
    plt.axvline(x=12,color='r')
    
    axs.plot(data_points, target[0][160000:190000],'--r',label='$x_1$ target',linewidth=1.5)
    axs.plot(data_points, target[1][160000:190000],'violet',label='$x_2$ target',linewidth=1.5,linestyle="dashed")
    
    points_s2_c1 = np.array([data_points, signal[0]]).T.reshape(-1,1,2)
    points_s2_c2 = np.array([data_points, signal[1]]).T.reshape(-1,1,2)
    segments_s2_c1 = np.concatenate([points_s2_c1[:-1],points_s2_c1[1:]], axis=1)
    segments_s2_c2 = np.concatenate([points_s2_c2[:-1],points_s2_c2[1:]], axis=1)
    
    lc_3 = LineCollection(segments_s2_c1, cmap="viridis", linewidth=4)
    lc_4 = LineCollection(segments_s2_c2, cmap="viridis", linewidth=4)
    lc_3.set_array(data_points)
    lc_4.set_array(data_points)
    
    axs.add_collection(lc_3)
    axs.add_collection(lc_4)
    
    line = axs.add_collection(lc_3)
    axs.legend(loc="lower right");
    
    fig.colorbar(line, ax=axs,orientation='horizontal')
    #axs.autoscale_view()
    axs.set_xlabel('time(s)')

def main():
    
    feedback_signal=signal_1
    target_signal=target_1
    
    timed_signal(feedback_signal,target_signal,0,30,-2,2,data_test)
    timed_limit_cycle(feedback_signal,target_signal,-3,25,-70,8,data_test)
    

    
    
if __name__ == "__main__":
    main()




##print("Van der Pol:%.8f \t %.8f \n Quadratic:%.8f \t %.8f \n Lissajous:%.8f \t %.8f " %(signal1_MSE,signal1_col2_MSE,signal2_MSE,signal2_col2_MSE,signal3_MSE,signal3_col2_MSE))
plt.show()
