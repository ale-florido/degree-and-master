import numpy as np
from gr_pyhole import metric
class MyMetric(metric.SphericalMetric ) :
    ID ='MY'  # short ID of this metric
    #GAMMA-METRIC
    def __init__ (self , M=1) :
        super ( MyMetric , self ) . __init__ ()
        self .M=M
        self.gamma =1.5  # =5 parece mas achatado
        self.m=self.M/self.gamma
    def update (self) :
        r, th, p = self.x[1:4]
        F= (1-2*self.m/r)**(self.gamma)
        u=r**2-2*self.m*r
        d=u+self.m**2*np.sin(th)**2
        G=(u/d)**(self.gamma**2-1)
        H=(u**(self.gamma**2))/(d**(self.gamma**2-1))

        self .tt =-1/F
        self .rr=F/G
        self.thth=F/H
        self.pp=F/(u*np.sin(th)**2)

        self.dr_tt=-self.tt*self.gamma*(1-2*self.m/r)**(-1)*2*self.m/r**2
        q=(self.gamma**2-1)*(r-self.m)*self.m**2*np.sin(th)**2
        self.dr_rr=2*F*G**(-1)*(self.m*self.gamma+q/d)/(r**2-2*self.m*r)
        self.dth_rr=F*G**(-1)*(self.gamma**2-1)*2*self.m**2*np.sin(th)*np.cos(th)/d
        a=(self.gamma**2-1)/d-self.gamma**2/(r**2-2*self.m*r)
        self.dr_thth =2*F*H**(-1)*(self.m*self.gamma*(1-2*self.m/r)**(-1)/r**2+(r-self.m)*a)
        self.dth_thth=2*self.m**2*F*H**(-1)*((self.gamma**2-1)*np.sin(th)*np.cos(th))/d
        self.dr_pp =2*F*(self.m*self.gamma-r+self.m)/((r**2-2*self.m*r)**2*np.sin(th)**2)
        self.dth_pp =-2*F*np.cos(th)/((r**2-2*self.m*r)*np.sin(th)**3)

        return r<1.001*2*self.m