# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


##File just checks if the signals being produced are correct
# Load data
data   = pd.read_csv('input.csv',header=None)
data_2  = pd.read_csv('volterra.csv',header=None)
data_3   = pd.read_csv('2ndOrder.csv',header=None)
data_4  = pd.read_csv('NARMA.csv',header=None)

numbers = np.linspace(0, len(data.index),len(data.index))
numbers_2 = np.linspace(0, len(data_2.index),len(data_2.index))

# Plot
plt.figure(figsize=(6.8, 4.2))
x = range(len(data[0]))
plt.plot(numbers[0:5000], data[1][0:5000],color="b")
plt.plot(numbers_2[0:5000], data_2[0][0:5000],color="r")
plt.plot(numbers_2[0:5000], data_3[0][0:5000],color="g")
plt.plot(numbers_2[0:5000], data_4[1][0:5000],color="y")
#plt.xticks(x, data)
plt.show()

