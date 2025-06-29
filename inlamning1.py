import numpy as np
from matplotlib import pyplot as plt
#%% Uppgift 1
def approx_ln(x,n):
    #don't know if this is necessary
    #if type(x)==type(1e1)or type(x)==type(0):
    #    if x==0:
    #        return np.inf
    #    elif x<0:
    #        return np.nan
    a=(1+x)/2#              a_0=\frac{1+x}[2}
    g=np.sqrt(x)#           g_0=\sqrt{x}
    for i in range(n):
        a=(a+g)/2#          a_{i+1}=\frac{a_i \cdot g_i}{2}
        g=np.sqrt(a*g)#     g_{i+1}=\sqrt{a_{i+1} \cdot g_i}
    return (x-1)/a

#%% Uppgift 2
def plota2Ln(n):
    x=np.exp(np.linspace(-10,10,100))
    fig, ax = plt.subplots(1,2)
    fig.suptitle(f"Amount of iterations: {n}")
    ax[0].set_title("Values")
    ax[1].set_title("Difference")
    ax[0].set_xlabel("x-value")
    ax[1].set_xlabel("x-value")
    ax[0].set_ylabel("Value")
    ax[0].set_ylabel("Error")
    approx_ln_value=approx_ln(x,n)
    log_of_x=np.log(x)
    ax[0].plot(x,approx_ln_value,label="approx_ln(x)")
    ax[0].plot(x,log_of_x,label="ln(x)")
    ax[1].plot(x,np.abs(log_of_x-approx_ln_value),label="|approx_ln(x)-ln(x)|")

    ax[0].legend()
    ax[1].legend()
    ax[0].set_xscale('log')
    ax[1].set_xscale('log')
    plt.show()

plota2Ln(2)
plota2Ln(3)
plota2Ln(5)


#%% Uppgift 3
def plotaUppgift3(x):
    #formating
    plt.title(f"Error of the approximation of ln(x) when x={x}")
    plt.xlabel("Iteration n")
    plt.ylabel("Error")
    plt.legend()
    plt.yscale('log')

    #calculations
    n=np.array([i+1 for i in range(100)])
    approx_ln_value_of_n=np.array([approx_ln(x,i) for i in n])
    logValue=np.log(x)

    #plot
    plt.plot(n,np.abs(logValue-approx_ln_value_of_n),label="|approx_ln_n(1.41)-ln(1.41)|")
    plt.show()
    #we see when n \approx 20 the error \approx 5.55e-17 
    #and that is the lowest precision we can get with 0 excluded,
    #beacuse the code and the datatype is float (real double). 

plotaUppgift3(1.41)
#%% Uppgift 4


#%% Uppgift 5


#%% Uppgift 6
