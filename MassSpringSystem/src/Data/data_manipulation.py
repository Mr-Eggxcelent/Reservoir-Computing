import pandas
from sklearn import preprocessing
import numpy as np
import matplotlib.pyplot as plt

df = pandas.read_csv('matrix.csv')
scaler = preprocessing.StandardScaler().fit(df)
X_scaled = scaler.transform(df)

df_2 = pandas.read_csv('matrix_2.csv')
scaler_2 = preprocessing.StandardScaler().fit(df_2)
X_scaled_2 = scaler_2.transform(df_2)

df_3 = pandas.read_csv('matrix_3.csv')
scaler_3 = preprocessing.StandardScaler().fit(df_3)
X_scaled_3 = scaler_3.transform(df_3)

df_4 = pandas.read_csv('targetSignal.csv')
scaler_4 = preprocessing.StandardScaler().fit(df_4)
X_scaled_4 = scaler_4.transform(df_4)

numbers = np.linspace(0, 469999,469999)
numbers_2 = np.linspace(0, 399999,399999)
numbers_3 = np.linspace(0, 29999,29999)
numbers_4 = np.linspace(0, 20000,20000)

# Plot
plt.figure(figsize=(6.8, 4.2))
##plt.plot(numbers, X_scaled[:,1],'b')
##plt.plot(numbers_2, X_scaled_2[:,6],'r')
##plt.plot(numbers_2, X_scaled_2[:,1],'r')
##plt.plot(numbers_2, X_scaled_2[:,2],'r')
##plt.plot(numbers_2, X_scaled_2[:,3],'r')
##plt.plot(numbers_3, X_scaled_3[:,6],'y')
#plt.xticks(x, data)
##plt.plot(numbers_4, X_scaled_4[190000:210000,0],'r')
plt.plot(numbers_2, X_scaled_2[:,0],'r')
plt.plot(numbers_2, X_scaled_2[:,1],'y')
plt.plot(numbers_2, X_scaled_2[:,2],'r')
plt.plot(numbers_2, X_scaled_2[:,3],'b')
plt.plot(numbers_2, X_scaled_2[:,4],'r')
plt.plot(numbers_2, X_scaled_2[:,5],'r')
plt.plot(numbers_2, X_scaled_2[:,6],'r')
plt.plot(numbers_2, X_scaled_2[:,7],'r')
plt.plot(numbers_2, X_scaled_2[:,8],'r')
plt.plot(numbers_2, X_scaled_2[:,9],'r')
plt.plot(numbers_2, X_scaled_2[:,10],'r')
plt.plot(numbers_2, X_scaled_2[:,11],'r')
plt.plot(numbers_2, X_scaled_2[:,12],'r')
plt.plot(numbers_2, X_scaled_2[:,13],'r')
plt.plot(numbers_2, X_scaled_2[:,14],'r')
plt.plot(numbers_2, X_scaled_2[:,15],'r')
plt.plot(numbers_2, X_scaled_2[:,16],'r')
plt.plot(numbers_2, X_scaled_2[:,17],'r')
plt.plot(numbers_2, X_scaled_2[:,18],'r')
plt.plot(numbers_2, X_scaled_2[:,19],'r')




plt.show()

#df= pandas.DataFrame(X_scaled)
#df.to_csv('normalized.csv',index=False)
