import pandas
from sklearn import preprocessing
import numpy as np
import matplotlib.pyplot as plt


#file you wish to preprocess
file="input_2.csv"

df= pandas.read_csv(file)
scaler = preprocessing.StandardScaler().fit(df)
X_scaled= scaler.transform(df)

x_axis = np.linspace(0, len(df.index),len(df.index))
x_axis_2 = np.linspace(0, len(X_scaled),len(X_scaled))

# Plot
plt.figure(figsize=(6.8, 4.2))
plt.plot(x_axis[0:5000], df.iloc[0:5000, 0],'b')
plt.plot(x_axis_2[0:5000], X_scaled[0:5000],'g')

##np.savetxt(file, X_scaled, fmt ="%s", delimiter=",")

plt.show()

print("DONE")

