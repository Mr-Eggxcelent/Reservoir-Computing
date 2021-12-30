import numpy as np
from matplotlib import pyplot as plt
from scipy import integrate
import csv
from scipy.integrate import solve_ivp
import math


def lissajous (t=0):
    a=1
    b=2
    phase=np.pi/2
    x=np.sin(a*t+phase)
    y=np.sin(b*t)
    return np.array((x, y))


def morris_lecar(u, I_ext,t=0):
    # Taken from Table 1 of doi:10.1016/j.neucom.2005.03.006
    C_M = 20
    phi = 0.066667
    g_K = 8
    g_L = 2
    g_Ca = 4.0
    V_L = -60
    V_K = -84
    V_Ca = 120
    V_1 = -1.2
    V_2 = 18
    V_3 = 12
    V_4 = 17
    # Taken from Table 2 (class I) of doi:10.1016/j.neucom.2005.03.006
    # References from doi:10.1016/j.neucom.2005.03.006
    (V, N) = u
    M_inf = 0.5*(1 + np.tanh((V - V_1)/V_2)) # (2)
    N_inf = 0.5*(1 + np.tanh((V - V_3)/V_4)) # (3)
    tau_N = 1/(np.cosh((V - V_3)/(2*V_4))) # (4)
    # (1)
    dVdt = (-g_L*(V - V_L) - g_Ca*M_inf*(V - V_Ca) - g_K*N*(V - V_K) + I_ext)/C_M
    dNdt = phi*((N_inf - N)/tau_N)
    return np.array((dVdt, dNdt))

ml30 = lambda u,t: morris_lecar(u, 80,t)


def lokta_volterra(X, t=0):
    a = 1.1
    b = 0.4
    c = 0.4
    d = 0.1
    return np.array([ a*X[0] -   b*X[0]*X[1] ,
                  -c*X[1] + d*b*X[0]*X[1] ])


def second_order_system(y0,y1,u):
    return np.array([(0.4)*y0 +(0.4)*y0*y1+ 0.6*(u**3)+0.1])



def ode_van_der_Pol_sd(x,epsilon,dt,u):
    x_new_1= x[0] + dt*x[1]
    x_new_2= x[1] - dt*x[0] + dt*epsilon*(1-x[0]**2)*x[1] + dt*u

    return np.array((x_new_1,x_new_2))
    

def quad(x,t=0):
    x_new_1=    x[0]+x[1]- 5*x[0]*(x[0]**2+x[1]**2)
    x_new_2= -2*x[0]+x[1]- x[1]*(x[0]**2+x[1]**2)

    return np.array((x_new_1,x_new_2))
    

step_size = 0.001
init_phase = 20
data_time = 500
total_time = data_time+init_phase
tot_timestep=total_time/step_size
t = np.linspace(0,total_time,int(tot_timestep))

##Simply uncomment the signal you want. Will change them to functions later.

#### Quad
##X0= np.array([-0.28,-0.539])
##X_1, infodict = integrate.odeint(quad, X0, t, full_output=True)
##
####get rid of initial phase
##discard=init_phase/step_size+1
##X_trunc = X_1[int(init_phase/step_size+1):,:]
##
##plt.xlabel('timesteps []');
##plt.ylabel ('[ ]');
##plt.title('state variables x_1 and x_2')
##plt.plot(t[int(init_phase/step_size+1):int(init_phase/step_size+1)+6000],X_trunc[:6000,0])
##plt.plot(t[int(init_phase/step_size+1):int(init_phase/step_size+1)+6000],X_trunc[:6000,1])
####plt.plot(X_trunc_1[:,0],X_trunc_1[:,1])
##plt.show()

##np.savetxt("quad.csv", X_trunc, fmt ="%s", delimiter=",")

##### Vanderpol
##X= np.zeros((2,t.size))
##X[:,0] = np.array([3,3]); 
##idx = 0;
### simulate van der Pol equations
##for i in range(t.size-1):
##    idx = idx+1;
##    X[:,idx] = ode_van_der_Pol_sd((X[:,idx-1]),0.8,step_size,0.0);
##X_trunc = X[:,int(init_phase/step_size+1):]
##
##plt.xlabel('timesteps []');
##plt.ylabel ('[ ]');
##plt.title('state variables x_1 and x_2')
##plt.plot(t[int(init_phase/step_size+1):int(init_phase/step_size+1)+9000],X_trunc[0,:9000])
##plt.plot(t[int(init_phase/step_size+1):int(init_phase/step_size+1)+9000],X_trunc[1,:9000])
####plt.plot(X_trunc[0,:],X_trunc[1,:])
##plt.show()

##np.savetxt("vanderpol.csv", np.transpose(X_trunc), fmt ="%s", delimiter=",")

##Lokta Volterra
##X0= np.array([5,5])
##X_1, infodict = integrate.odeint(lokta_volterra, X0, t, full_output=True)
##
###get rid of initial phase
##discard=init_phase/step_size+1
##X_trunc_1 = X_1[int(init_phase/step_size+1):,:]
##
##plt.xlabel('timesteps []');
##plt.ylabel ('[ ]');
##plt.title('state variables x_1 and x_2')
##plt.plot(t[int(init_phase/step_size+1):],X_trunc_1[:,0])
##plt.plot(t[int(init_phase/step_size+1):],X_trunc_1[:,1])
####plt.plot(X_trunc_1[:,0],X_trunc_1[:,1])
##plt.show()

##np.savetxt("lokta_volterra.csv", X_trunc_1, fmt ="%s", delimiter=",")

##morris_lecar
##X0= np.array([-25.050458314,0.3])
##X_1, infodict = integrate.odeint(ml30, X0, t, full_output=True)
##
####get rid of initial phase
##discard=init_phase/step_size+1
##X_trunc_1 = X_1[int(init_phase/step_size+1):,:]
##
##plt.xlabel('timesteps []');
##plt.ylabel ('[ ]');
##plt.title('state variables x_1 and x_2')
####plt.plot(t[int(init_phase/step_size+1):],X_trunc_1[:,0])
####plt.plot(t[int(init_phase/step_size+1):],X_trunc_1[:,1])
##plt.plot(X_trunc_1[:,0],X_trunc_1[:,1])
##plt.show()     

##np.savetxt("morris_lecar.csv", X_trunc_1, fmt ="%s", delimiter=",")

####Lissajous 
##X_1=lissajous(t)
### get rid of initial phase
##discard=init_phase/step_size+1
##X_trunc = X_1[:,int(init_phase/step_size+1):]
##
##
##plt.xlabel('timesteps []');
##plt.ylabel ('[ ]');
##plt.title('state variables x_1 and x_2')
####plt.plot(t[int(init_phase/step_size+1):],X_trunc[0,:])
####plt.plot(t[int(init_phase/step_size+1):],X_trunc[1,:])
##plt.plot(np.transpose(X_trunc)[:,0],np.transpose(X_trunc)[:,1])
##plt.show()     

##np.savetxt("lissajous.csv",np.transpose(X_trunc), fmt ="%s", delimiter=",")

##### 2nd order
##X= np.zeros(t.size) 
##f1 = 2.11; 
##f2 = 3.73;
##f3 = 4.33;
##u=np.zeros(t.size)
##u = 0.1*(np.sin(2*math.pi*f1*t)*np.sin(2*math.pi*f2*t)*np.sin(2*math.pi*f3*t));
##
##for i in range(2,t.size):
##    X[i] = second_order_system((X[i-1]),X[i-2],u[i-1]);
##
### get rid of initial phase
##discard=init_phase/step_size+1
##X_trunc_2 = X[int(init_phase/step_size+1):]
##
##plt.xlabel('timesteps []');
##plt.ylabel ('[ ]');
##plt.title('state variables x_1 and x_2')
##plt.plot(t[int(init_phase/step_size+1):int(init_phase/step_size+1)+5000],X_trunc_2[:5000])
##plt.show()

##np.savetxt("2ndOrder.csv", X_trunc_2, fmt ="%s", delimiter=",")


print("DONE")
