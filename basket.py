import numpy as np
import matplotlib.pyplot as plt

class Basketball:
    def __init__(self,init_angle,z=1,size=5000,init_speed=9,AIR_resistans=True):
        self.x=np.zeros(size)
        self.y=np.zeros(size)
        self.v_x=np.zeros(size)
        self.v_y=np.zeros(size)
        
        self.mass=0.6
        self.diameter=0.24

        self.g=9.81
        if AIR_resistans:
            self.rho=1.23
        else:
            self.rho=0
        self.cw=0.45
        self.pi=np.pi

        self.hfac=(z)/(size-1)

        self.y_b=.0
        self.x_b=.0
        self.set_inits(init_angle=init_angle/180*self.pi,init_speed=init_speed)
        self.init_basket()
        self.fillV()
        self.fill_xy()
        
        self.label=None

    def set_inits(self,init_x=0,init_y=1.75,init_angle=0,init_speed=0):
        self.x[0]=init_x
        self.y[0]=init_y
        self.v_x[0]=np.cos(init_angle)*init_speed
        self.v_y[0]=np.sin(init_angle)*init_speed
        
    def init_basket(self,init_y_b=3.05,init_x_b=2):    
        self.y_b=3.05
        self.x_b=init_x_b
    
    def fillV(self,startindex=1):
        self.v_x[startindex]=self.v_x[startindex-1]+self.hfac*self.get_x_acc(startindex-1)
        self.v_y[startindex]=self.v_y[startindex-1]+self.hfac*self.get_y_acc(startindex-1)
        for i in range(startindex+1,self.x.shape[0]):
            self.v_x[i]=self.v_x[i-1]+3/2*self.hfac*self.get_x_acc(i-1)-1/2*self.hfac*self.get_x_acc(i-2)#Euler
            self.v_y[i]=self.v_y[i-1]+3/2*self.hfac*self.get_y_acc(i-1)-1/2*self.hfac*self.get_y_acc(i-2)#Euler
 
    def fill_xy(self,startindex=1):
        self.x=self.integrate(self.x,self.v_x*self.hfac,startindex=startindex)
        self.y=self.integrate(self.y,self.v_y*self.hfac,startindex=startindex)

    def integrate(self,array,der_array,startindex=1):
        for i in range(startindex,array.shape[0]):
            array[i]=array[i-1]+(der_array[i]+der_array[i-1])/2
        return array

    def get_cos_alpha(self,index):
        return self.v_x[index]/(self.v_x[index]**2+self.v_y[index]**2)**0.5
    
    def get_sin_alpha(self,index):
        return self.v_y[index]/(self.v_x[index]**2+self.v_y[index]**2)**0.5
    
    def get_x_acc(self,index):
        return (-self.getF(index)*self.get_cos_alpha(index))/self.mass
    
    def get_y_acc(self,index):
        return (-self.getF(index)*self.get_sin_alpha(index))/self.mass-self.g
    
    def getF(self,index):
        return 1/2*self.pi/4 * self.rho*self.cw * self.diameter**2*(self.v_x[index]**2+self.v_y[index]**2)

    def increaseSize(self,newSize):
        dif=newSize-self.x.shape[0]
        self.x=np.array(list(self.x)+[0]*dif)
        self.y=np.array(list(self.y)+[0]*dif)
        self.v_x=np.array(list(self.v_x)+[0]*dif)
        self.v_y=np.array(list(self.v_y)+[0]*dif)

    def get_x_last(self):
        return self.x[self.x.shape[0]-1]
    
    def get_y_last(self):
        return self.y[self.y.shape[0]-1]
    
    def get_yder_last(self):
        return (self.y[self.y.shape[0]-1]-self.y[self.y.shape[0]-2])*self.hfac
    
class ploterC:
    def __init__(self):
        pass
    
    def plota(self,basketballs):
        for basketball in basketballs:
            if basketball!=None:
                if basketball.label==None:
                    plt.plot(basketball.x,basketball.y)
                else:
                    plt.plot(basketball.x,basketball.y,label=basketball.label)

        plt.scatter([2],[3.05])
        plt.ylim(bottom=0)
        plt.xlabel("x")
        plt.ylabel("y")
        plt.legend()
        plt.show()

class Newton:
    def __init__(self,init_a,init_z,h_a=1e-3,h_z=1e-3,amountofbasketballs=1):
        self.h_a=h_a
        self.h_z=h_z
        self.xb=2
        self.yb=3.05
        self.Class=Basketball
        self.basketballs=[None]*max(amountofbasketballs,1)
        self.a=[None]*max(amountofbasketballs,1)
        self.z=[None]*max(amountofbasketballs,1)
        
        self.a[0]=init_a
        self.z[0]=init_z        
        self.basketballs[0]=self.Class(init_a,init_z)

    def increase(self,index,inlarningspeed=1):
        error_last=self.error(self.basketballs[index-1])

        a_der=(self.error(self.Class(self.a[index-1]+self.h_a,self.z[index-1]))-error_last)/self.h_a
        z_der=(self.error(self.Class(self.a[index-1],self.z[index-1]+self.h_z))-error_last)/self.h_z

        #decrease lambda
        #a_der*=min(1,abs(7/a_der))
        #z_der*=min(1,abs(7/z_der))
        print(a_der,z_der)
        if abs(a_der)<1e-2:
            a_der=1e-2*a_der/abs(a_der)
        if abs(z_der)<1e-2:
            z_der=1e-2*z_der/abs(z_der)

        #print(f"\t{a_der} {z_der}")
        self.a[index]=self.a[index-1]-inlarningspeed*error_last/a_der#partielderivata
        self.z[index]=self.z[index-1]-inlarningspeed*error_last/z_der#partielderivata
        
        self.a[index]=max(0,min(90,self.a[index]))
        self.z[index]=max(0.2,self.z[index])

        self.basketballs[index]=Basketball(self.a[index],self.z[index])

    def error(self,basketball):
        #print(basketball.get_x_last(),self.xb,basketball.get_y_last(),self.yb)
        return np.abs(basketball.get_x_last()-self.xb)+np.abs(basketball.get_y_last()-self.yb)
        return (np.abs(basketball.get_x_last()-self.xb)**2+np.abs(basketball.get_y_last()-self.yb)**2)
    
    def maxlearn(self,dif,maxdif):
        return min(max(dif,-maxdif),maxdif)
    
    def increaseSize(self,newSize):
        dif=newSize-len(self.basketballs)
        if dif>0:
            self.basketballs=list(self.basketballs)+[None]*dif
            self.a=list(self.a)+[None]*dif
            self.z=list(self.z)+[None]*dif

    

#%%
amountofbasketballs=10
basketballs=[None]*amountofbasketballs
basketballs[0]=Basketball(45,1)
basketballs[1]=Basketball(30,1)
basketballs[2]=Basketball(60,1)
basketballs[3]=Basketball(80,2)
basketballs[4]=Basketball(82,2)
basketballs[5]=Basketball(81,1.6)
basketballs[6]=Basketball(81.206,1.567)
basketballs[7]=Basketball(82,2,AIR_resistans=False)
basketballs[7].label="Air Resistans=0"
basketballs[6].label="Analytic"

ploter=ploterC()
ploter.plota(basketballs)

n=Newton(65,1.5,amountofbasketballs=160)
for i in range(len(n.basketballs)-1):
    inlspeed=.2*(1-i/len(n.basketballs))+0.1
    n.increase(i+1,inlspeed)
    print(i,n.a[i],n.z[i],n.error(n.basketballs[i]))

basketballs=[None]*10
basketballs[0]=Basketball(n.a[len(n.a)-1],n.z[len(n.z)-1])
basketballs[0].label="Analytic"
basketballs[1]=Basketball(70,1.5)
basketballs[2]=Basketball(85,1.5)
basketballs[1].label="70 deg"
basketballs[2].label="85 deg"



ploter.plota(basketballs)

print("done")
#%%

def f1(t,y):
    return 0
def Euler(f,T=5,N=10000):
    h=T/N
    u=np.zeros(N)
    t=np.linspace(0,T,N,endpoint=False)
    u[0]=f(0,0)

    for i in range(N-1):
        u[i+1]=u[i]+h*f(t[i],u[i])
    return t,u
