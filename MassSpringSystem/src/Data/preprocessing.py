import pandas
from sklearn import preprocessing
import numpy as np
import matplotlib.pyplot as plt

##df = pandas.read_csv('matrix.csv')
##scaler = preprocessing.StandardScaler().fit(df)
##X_scaled = scaler.transform(df)
##
##df_2 = pandas.read_csv('matrix_2.csv')
##scaler_2 = preprocessing.StandardScaler().fit(df_2)
##X_scaled_2 = scaler_2.transform(df_2)
##
##df_3 = pandas.read_csv('matrix_3.csv')
##scaler_3 = preprocessing.StandardScaler().fit(df_3)
##X_scaled_3 = scaler_3.transform(df_3)

df_4 = pandas.read_csv('vanderpol.csv')
scaler_4 = preprocessing.StandardScaler().fit(df_4)
X_scaled_4 = scaler_4.transform(df_4)

numbers = np.linspace(0, 469999,469999)
numbers_2 = np.linspace(0, 399999,399999)
numbers_3 = np.linspace(0, 29999,29999)
numbers_4 = np.linspace(0, 499998,499998)

# Plot
plt.figure(figsize=(6.8, 4.2))
##plt.plot(numbers, X_scaled[:,1],'b')
##plt.plot(numbers_2, X_scaled_2[:,6],'r')
##plt.plot(numbers_2, X_scaled_2[:,1],'r')
##plt.plot(numbers_2, X_scaled_2[:,2],'r')
##plt.plot(numbers_2, X_scaled_2[:,3],'r')
##plt.plot(numbers_3, X_scaled_3[:,6],'y')

plt.plot(numbers_4, X_scaled_4[:,0],'r')
plt.plot(numbers_4, df_4.iloc[:, 0],'b')
np.savetxt("vanderpol.csv", X_scaled_4, fmt ="%s", delimiter=",")

plt.show()

print("DONE")

#df= pandas.DataFrame(X_scaled)
#df.to_csv('normalized.csv',index=False)
