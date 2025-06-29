import numpy as np
from matplotlib import pyplot as plt
#%% Uppgift 1
def approx_ln(x,n):
    #don't know if this is necessary
    if x==0:
        return np.inf
    elif x<0:
        return np.nan
    a=(1+x)/2#              a_0=\frac{1+x}[2}
    g=np.sqrt(x)#           g_0=\sqrt{x}
    for i in range(n):
        a=(a+g)/2#          a_{i+1}=\frac{a_i \cdot g_i}{2}
        g=np.sqrt(a+g)#     g_{i+1}=\sqrt{a_{i+1} \cdot g_i}
    return (x-1)/a

#%% Uppgift 2
def plota2Ln(n):
    x=np.exp(np.linspace(-10,10,100))
    fig, ax = plt.subplots(2,1)
    ax[0].set_title("Values")
    ax[1].set_title("Difference")
    ax[0].set_xlabel("Iteration")
    ax[1].set_xlabel("Iteration")
    ax[0].set_xlabel("Value")
    ax[0].set_xlabel("Error")
    approx_ln_value=approx_ln(x,n)
    log_of_x=np.log(x)
    ax[0].plot(approx_ln_value,label="approx_ln(x)")
    ax[0].plot(log_of_x,label="ln(x)")
    ax[1].plot(np.abs(log_of_x-approx_ln_value),"|approx_ln(x)-ln(x)|")

    ax[0].legend()
    ax[1].legend()
    plt.show()

plota2Ln(2)
#%% Uppgift 3


#%% Uppgift 4


#%% Uppgift 5


#%% Uppgift 6
