# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


# Load data
##data = pd.read_csv('outputsignal_two.csv',header=None)
target = pd.read_csv('NARMA.csv',header=None)
numbers = np.linspace(0, 499001,499001)

print(numbers)

# Plot
plt.figure(figsize=(6.8, 4.2))
x = range(len(target[0]))
plt.plot(numbers, target)
#plt.xticks(x, data)
plt.show()

