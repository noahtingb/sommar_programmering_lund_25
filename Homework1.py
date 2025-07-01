#%%imports
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
    #transform x to an np.array
    if type(x)==type(float(.1)):
        x=np.array([x])
    
    #initiate
    a=np.zeros((n+1,x.shape[0]))

    #a_0
    a[0,:]=(1+x)/2#              a_0=\frac{1+x}[2}
    g=np.sqrt(x)#           g_0=\sqrt{x}

    #a for n>=i+1>0
    for i in range(n):
        a[i+1,:]=(a[i,:]+g)/2#          a_{i+1}=\frac{a_i \cdot g_i}{2}
        g=np.sqrt(a[i+1,:]*g)#     g_{i+1}=\sqrt{a_{i+1} \cdot g_i}
    return np.array([list((x-1)/a[i,:]) for i in range(n+1)])

#%% Uppgift 2
def plota2Ln(n):
    x=np.exp(np.linspace(-10,10,100))#a linspace for x
    fig, ax = plt.subplots(1,2)#two subplots
    
    #labls and titles and scales
    fig.suptitle(f"Amount of iterations: {n}")
    ax[0].set_title("Values")
    ax[1].set_title("Difference")
    ax[0].set_xlabel("x-value")
    ax[1].set_xlabel("x-value")
    ax[0].set_ylabel("Value")
    ax[0].set_ylabel("Error")
    ax[0].set_xscale('log')
    ax[1].set_xscale('log')
    
    #approximate ln and calculate ln
    approx_ln_value=approx_ln(x,n)[n,:]
    log_of_x=np.log(x)

    #plot the correct ln and the approximated ln dependent on x
    ax[0].plot(x,approx_ln_value,label="approx_ln(x)")
    ax[0].plot(x,log_of_x,label="ln(x)")

    #plot the absoulut difference
    ax[1].plot(x,np.abs(log_of_x-approx_ln_value),label="|approx_ln(x)-ln(x)|")

    #legends
    ax[0].legend()
    ax[1].legend()

    #show
    plt.show()

#dont understand exactly how you want the plots
plota2Ln(2) 
plota2Ln(3) 
plota2Ln(5) 


#%% Uppgift 3
def plotaUppgift3(x):
    #formating
    plt.title(f"Error of the approximation of ln(x) when x={x}")
    plt.xlabel("Iteration n")
    plt.ylabel("Error")
    plt.yscale('log')

    #calculations
    n=np.array([i+1 for i in range(100)])
    approx_ln_value_of_n=approx_ln(np.array([x]),100)[1:,0]
    logValue=np.log(x)

    #plot
    plt.plot(n,np.abs(logValue-approx_ln_value_of_n),label="|approx_ln_n(1.41)-ln(1.41)|")
    plt.legend()
    plt.show()
    #we see when n \approx 20 the error \approx 5.55e-17 
    #and that is the lowest precision we can get with 0 excluded,
    #beacuse the code and the datatype is float (real double). 

#calc x
plotaUppgift3(1.41)

#%% Uppgift 4
def fast_approx_ln(x,n_max,**kwarg):#kwarg for the return format (the comments below are in Swedish)
    #transform x to an np.array
    if type(x)==type(.1):#lättare med anpassning ifall vi använder np array för x (C snabbare språk)
        x=np.array([x])#omvandla till float till np array 

    #initiate an 3d_array for x with (k,n,x)
    d=np.zeros((n_max+1,n_max+1,x.shape[0]))#inisierar en (n+1)x(n+1)x(x.shape) "matris"
    
    #create d[0,:,x] or a in the last assigment
    d[0,0,:]=(1+x)/2#initial värde för d[0,0] (a[0])
    g=np.sqrt(x)#    initialvärde för g för skapandet av vektorn a
    for n in range(1,n_max+1):#loppa 
        d[0,n,:]=(d[0,n-1,:]+g)/2#  skapa a[n] från a[n-1] och g 
        g=np.sqrt(d[0,n,:]*g)#     uppdatera g (g[n] från g[n-1] och a[n])
    
    #calculate d for k>0 according to the formula up to d[n_max,n_max]
    for k in range(1,n_max+1):#loppa k mellan 1 och n_max
        for n_i in range(k,n_max+1):#loppa n  mellan k och n_max (n>=k) men alla n_max>=n>k behövs
            d[k,n_i,:]=(d[k-1,n_i,:]-4**(-k)*d[k-1,n_i-1,:])/(1-4**(-k))#d[k,n] från d[k-1,n] och d[k-1,n-1]
    
    #return the approximated dependent on which chossen format
    if kwarg.get("output_np")==True:#om outputen ska vara i en np
        if kwarg.get("output_array")==True:#return an array with all x and the different iterations
            return np.array([(x-1)/d[i,i,:] for i in range(n_max+1)])
        return (x-1)/d[n_max,n_max,:]#else return the array with all x
    return float(((x-1)/d[n_max,n_max,:])[0])#return


#%% Uppgift 5
def plotaUppgift5(**kwarg):
    #formating
    plt.title(f"Error behivour of the accelerated Carlsson method for the log")#set a title
    plt.xlabel("x")#xlabel
    plt.ylabel("error")#error
    plt.yscale('log')
    plt.ylim((1e-19,1e-5))

    #amount of points
    points=kwarg.get("points",1000)
    if type(points)!=type(int(1)):points=1000

    #calculations
    x_linspace=np.linspace(20/points,20,points)
    approxed_logValues=fast_approx_ln(x_linspace,6,output_np=True,output_array=True)
    correct_logValues=np.log(x_linspace)

    #plot
    for i in range(2,7):
        plt.scatter(x_linspace,np.abs(correct_logValues-approxed_logValues[i,:]),s=5,label=f"Iteration {i}")

    #legend+show
    plt.legend()
    plt.show()

plotaUppgift5()