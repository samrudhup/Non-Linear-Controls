import numpy as np
from math import sin
from math import cos

#Defining a class dynamics
class dynamics():
    def __init__(self, alpha=0.25*np.ones(2,dtype=np.float32), beta=0.1*np.ones(2, dtype=np.float32), gamma=0.01*np.ones(5,dtype=np.float32)):
        
        #defining the gaines
        self.alpha=np.diag(alpha)
        self.beta=np.diag(beta)
        self.Gamma=np.diag(gamma)
        
        #defining the rigid body parameters
        self.m=np.array([2.0, 2.0],dtype=np.float32)
        self.l=np.array([0.5, 0.5],dtype=np.float32)
        self.mBounds=np.array([1.0, 3.0], dtype=np.float32)
        self.lBounds=np.array([0.25, 0.75],dtype=np.float32)
        self.g=9.8
        
        #defining the desired trajectory parameters
        self.phidMag=np.array([np.pi/8, np.pi/4],dtype=np.float32)
        self.freq=0.2
        self.a=np.pi/2
        self.b=np.array([np.pi/2, np.pi/4],dtype=(np.float32))
        
        #Initialize state
        self.phi,_,_=self.getDesiredstate(0.0) #seting the inititial angles to the desired angles
        self.phiD=np.zeros(2,dtype=np.float32) #initial angular velocity 
        self.phiDD=np.zeros(2,dtype=(np.float32))   #initial angular acceleration
        
        # unknown parameters
        self.theta = self.getTheta(self.m, self.l)
        self.thetaH= self.getTheta(self.mBounds[0]*np.ones(2,dtype=np.float32),self.lBounds[0]*np.ones(2,dtype=np.float32))
        self.tau=np.zeros(2,np.float32)
        
        
    def getDesiredstate(self,t):
        #To get the desired anglular displacements
        # self.phiMag corrected in all these equations 
        phid=np.array([self.phidMag[0]*sin(2*np.pi*self.freq*t - self.a) - self.b[0],
                     self.phidMag[1]*sin(2*np.pi*self.freq*t-self.a)+self.b[1]], dtype=np.float32)
        
        #To get the desired angular velocity
        phiDd=np.array([2*np.pi*self.freq*self.phidMag[0]*cos(2*np.pi*self.freq*t- self.a),
                       2*np.pi*self.freq*self.phidMag[1]*cos(2*np.pi*self.freq*t- self.a)], dtype=np.float32)
        
        #To get the desired angular acceleration
        phiDDd=np.array([-((2*np.pi*self.freq)**2)*self.phidMag[0]*sin(2*np.pi*self.freq*t - self.a),
                        -((2*np.pi*self.freq)**2)*self.phidMag[1]*sin(2*np.pi*self.freq*t - self.a)], dtype=np.float32)
        
        return phid, phiDd, phiDDd
   
    
    def getTheta(self,m,l):
    
        theta = np.array([(m[0]+m[1])*l[0]**2+m[1]*l[1]**2,
                          m[1]*l[0]*l[1],
                          m[1]*l[1]**2,
                          (m[0]+m[1])*l[0],
                          m[1]*l[1]],dtype=np.float32)
        return theta
    
    #Returns the inertial matrix
    def getM(self,m,l,phi):
        
        m1=m[0]
        m2=m[1]
        c2=cos(phi[1])
        l1=l[0]
        l2=l[1]
        #there was an error in the M equation-corrected
        M= np.array([[(m1*l1**2) + m2*((l1**2) + 2*l1*l2*c2 +l2**2),(m2*(l1*l2*c2 + l2**2))],
                    [m2*(l1*l2*c2 + l2**2), m2*l2**2]], dtype=np.float32)
        return M
                    
    # Returns the Centripital Coriolis matirx                
    def getC(self,m,l,phi,phiD):
        m2=m[1]
        s2=sin(phi[1])
        l1=l[0]
        l2=l[1]
        #corrected C equation.
        C=np.array([-2*m2*l1*l2*s2*phiD[0]*phiD[1] - m2*l1*l2*s2*(phiD[1]**2),
                    m2*l1*l2*s2*phiD[0]**2], dtype=np.float32)
        return C
    #Returns the Grativational Matrix    
    
    def getG(self, m, l, phi):
        m1=m[0]
        m2=m[1]
        c1=cos(phi[0])
        c12=cos(phi[0]+phi[1])
        l1=l[0]
        l2=l[1]
                        
        G=np.array([((m1+m2)*self.g*l1*c1 + m2*self.g*l2*c12),
                    (m2*self.g*l2*c12)], dtype=np.float32)    
        return G
    
    #Returns the inertial matrix regressor
    def getYM(self, vphi, phi):
        
        c2=cos(phi[1])
        
        #YM equation had errors - corrected
        
        YM=np.array([[vphi[0], 2*c2*vphi[0] + c2*vphi[1], vphi[1], 0.0, 0.0],
                      [0.0, c2*vphi[0], vphi[0]+vphi[1], 0.0, 0.0]], dtype=np.float32)
        
        return YM
    
    #Returns the centripital coriolis matrix regressor
    def getYC(self, phi, phiD):
        
        s2=sin(self.phi[1])
        
        YC= np.array([[ 0.0, -2*s2*phiD[0]*phiD[1] - s2*phiD[1]**2, 0.0, 0.0, 0.0],
                       [0.0, s2*phiD[0]**2, 0.0, 0.0, 0.0]], dtype=np.float32)
        return YC
    
    #Returns the Gravitational matrix regressor
    def getYG(self, phi):
        c1=cos(phi[0])
        c12=cos(phi[0]+phi[1])
        
        YG=np.array([[ 0.0, 0.0, 0.0, self.g*c1, self.g*c12],
                      [0.0, 0.0, 0.0, 0.0, self.g*c12]], dtype=(np.float32))
        return YG
    
    #Returns the M_dot matrix regressor
    def getYM_dot(self, phi, phiD, r):
        s2=sin(phi[1])
        
        YM_dot=np.array([[ 0.0, -2*s2*phiD[1]*r[0] - s2*phiD[1]*r[1], 0.0, 0.0, 0.0],
                          [0.0, -s2*phiD[1]*r[0], 0.0, 0.0, 0.0]], dtype=(np.float32))
        return YM_dot
    
    #returns the state
    def getState(self, t):
        
        return self.phi, self.phiD, self.phiDD, self.thetaH, self.theta
    
    #returns the error states
    def getErrorStates(self, t):
        
        #gets the desired states
        phid, phiDd, phiDDd= self.getDesiredstate(t)
        
        #gets the errors
        e = phid - self.phi
        eD = phiDd - self.phiD
        r = eD + self.alpha@e
        
        #calculate thetatilda
        
        thetatilda=self.theta-self.thetaH
        
        return e, eD, r, thetatilda
    
    #returns the inputs and update law
    def getTauThetaHD(self, t):
        #gets the desired states
        _,_, phiDDd= self.getDesiredstate(t)
        
        #get the error
        e,eD,r,_=self.getErrorStates(t)
        
        #get the regressors
        vphi = phiDDd + self.alpha@eD
        YM = self.getYM(vphi, self.phi)
        YC =self.getYC(self.phi, self.phiD)
        YG =self.getYG(self.phi)
        YM_dot=self.getYM_dot(self.phi, self.phiD, r)
        Y = YM +YC +YG + 0.5*YM_dot
        
        #calculating the contoroller's update law
        taufb=e + self.beta@r
        tauff=Y@self.thetaH
        tau = taufb + tauff
        
        #update the update law
        thetaHD = self.Gamma@Y.T@r
        return tau, thetaHD, tauff, taufb
    
    #steping the system
    def step(self, dt, t):
        #get the dynamics
        M = self.getM(self.m, self.l, self.phi)
        C =self.getC(self.m, self.l, self.phi, self.phiD)
        G = self.getG(self.m, self.l, self.phi)
        
        #get the input and the update law
        tau, thetaHD,_,_ =self.getTauThetaHD(t)
        
        #calculate the dynamics using the input
        self.phiDD =np.linalg.inv(M)@(-C-G+tau)
        
        #update the internal states
        self.phi += dt*self.phiD
        self.phiD += dt*self.phiDD
        self.thetaH += dt*thetaHD
            
            
    
    
                    
    
    
    
                          
        