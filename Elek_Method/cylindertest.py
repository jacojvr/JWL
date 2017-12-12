'''
cylindertest.py:
------------------------------------------------------------------------------

  ENERGY BALANCE METHOD AS MENTIONED IN
    >> TECHNICAL_NOTE 2016-07-20-1 <<

------------------------------------------------------------------------------

  P.M. Elek et. al. (2015)
  "Determination of Detonation Products Equation of State from Cylinder Test:
   Analytical Model and Numerical Analysis". Thermal Science 19(1) pp 38-45
   
------------------------------------------------------------------------------

'''
import numpy as np
from scipy.optimize import curve_fit,differential_evolution


def is_number(s):
    '''supplementary function to check if string is a number
    '''
    try:
        complex(s) # for int, long, float and complex
    except ValueError:
        return False
    return True


class jwl:
    
    def __init__(self,txtname):
        '''
        
This function reads the cylinder experiment text file
    
    the text file should be structured as follows :
    
    
    **Cylinder properties:
        *radii (r10,r20)
        [mm],[mm]
        #value,#value
        *Cylinder Density
        [kg/m3]
        #value
        *Cylinder Johnson-Cook (A,B,C,n)
        [MPa],[MPa],-,-
        #value,#value,#value,#value
        
    **Explosive properties:
        *Density of the explosive
        [kg/m3]
        #value
        *Detonation velocity
        [mm/1e-6s]
        #value
        *Pressure at CJ point
        [GPa]
        #value
        *Detonation Heat
        [GPa]
        #value
    
    **Expansion data
        *Time, Radius
        [1e-6s], [mm]
        #value,#value
        #value,#value
        ...
        ...
        #value,#value
'''
        #
        # The text file float data:
        filedata = open(txtname).read().split('**Expansion') #split in two (single vals and lists)
        # The data dictionary
        self.props = {}
        # The cylinder values / keys in order they appear in the text file
        Ckeys = ['Rin','Rout','RhoCu','JC_A','JC_B','JC_C','JC_n']
        # The explosive values / keys as they appear in the text file
        Ckeys += ['Rho0','D','PCJ','E0']
        # extract all single value floats into data dictionary using the keynames above:
        kn=0
        for ln in filedata[0].split('\n'):
            #could be that the data is split using (, ; or spaces)
            if ln.split(',').__len__()>1: 
                lnsplit = ln.split(',')
            elif ln.split(';').__len__()>1:
                lnsplit = ln.split(';')
            else:
                lnsplit = ln.split(' ')
            for val in lnsplit:
                if is_number(val):
                    if kn==Ckeys.__len__():
                        print "*"*80+"\n*WARNING - MISMATCH between expected keys and input text file"
                        print "\ttry help(readfile) for expected input file format\n"+"*"*80
                        yn = raw_input("\t>>> would you like to see help on the input file structure? (Y/n)").lower()
                        if 'n' in yn:
                            return
                        else:
                            help(self.__init__)
                            return
                        
                        
                    self.props[Ckeys[kn]]=float(val)
                    kn+=1
        if not 'E0' in self.props.keys():
            print "*"*80+"\n*WARNING - MISMATCH between expected keys and input text file"
            print "\ttry help(readfile) for expected input file format\n"+"*"*80
            yn = raw_input("\t>>> would you like to see help on the input file structure? (Y/n)").lower()
            if 'n' in yn:
                return
            else:
                print self.__init__.__doc__
                return
        # the expansion list
        xydata = []
        for ln in filedata[1].split('\n'):
            #could be that the data is split using (, ; or spaces)
            if ln.split(',').__len__()>1:
                lnsplit = ln.split(',')
            elif ln.split(';').__len__()>1:
                lnsplit = ln.split(';')
            else:
                lnsplit = ln.split(' ')
            for val in lnsplit:
                if  is_number(val):
                    xydata+=[float(val)]
                
        self.xdata = np.array(xydata[::2])
        self.ydata = np.array(xydata[1::2])
        #
        # constants:
        r10,r20,rho0,rhoCu = [self.props[kn] for kn in ['Rin','Rout','Rho0','RhoCu']]
        #
        # mass per length of cylinder and explosive [kg/m]
        self.props['M'] = 1e-6*np.pi*rhoCu*(r20**2-r10**2)
        self.props['C'] = 1e-6*np.pi*rho0*r10**2






    def fit_r2(self,timeshift=False):
        '''
        
    This function fits the outer radius expansion function to the xdata vs ydata.
        
    >> By default the function does not take a time shift into account,
       it is possible however to call it using:
       
            fit_r2(True),
       
       thereby allowing determination of an initial time shift 
       applied once all possible function formulations are checked.
       
        '''
        # constants:
        # cylinder properties:
        r10,r20,rhoCu,M = [self.props[kn] for kn in ['Rin','Rout','RhoCu','M']]
        # cylinder yield stress parameters (Johnson-Cook)
        JC_A,JC_B,JC_C,JC_n = [self.props[kn] for kn in ['JC_A','JC_B','JC_C','JC_n']]
        # detonation properties
        rho0,D,E0,C = [self.props[kn] for kn in  ['Rho0','D','E0','C']]

## here get the expansion function that best fits the data and assign to the jwl object
        ## Eq.(5) in Elek et al.
        def expand_F1(t,a0,v_inf,sig,t0=0):
            if t is None: t=self.xdata
            if not timeshift: t0=0.
            g_t = (1+t-t0)**sig-1
            f1t = v_inf*(t-t0)*g_t
            f1b = 2*v_inf*sig/a0 + g_t
            return f1t/f1b
        #d(F1)/dt
        def veloc_F1(t,a0,v_inf,sig):
            if t is None: t=self.xdata
            g_t = (1+t)**sig-1
            dg_t = sig*(1+t)**(sig-1)
            f1t = v_inf*t*g_t
            df1t = v_inf*g_t+v_inf*t*dg_t
            f1b = 2*v_inf*sig/a0 + g_t
            df1b = 2*v_inf*sig*(sig+1)/a0 + dg_t
            return df1t/f1b-f1t*f1b/(f1b**2)
        #d2(F1)/dt2
        def accel_F1(t,a0,v_inf,sig):
            if t is None: t=self.xdata
            g_t = (1+t)**sig-1
            dg_t = sig*(1+t)**(sig-1)
            f1t = v_inf*t*g_t
            df1t = v_inf*g_t+v_inf*t*dg_t
            f1b = 2*v_inf*sig/a0 + g_t
            df1b = 2*v_inf*sig*(sig+1)/a0 + dg_t
            return df1t/f1b-f1t*f1b/(f1b**2)
## Eq.(7) in Elek et al.
        def expand_F2(t,a1,b1,a2,b2,t0=0):
            if t is None: t=self.xdata
            if not timeshift: t0=0.
            fa1 = a1*(b1*(t-t0)-(1-np.exp(-b1*(t-t0))))
            fa2 = a2*(b2*(t-t0)-(1-np.exp(-b2*(t-t0))))
            return fa1+fa2
        #d(F1)/dt
        def veloc_F2(t,a1,b1,a2,b2):
            if t is None: t=self.xdata
            dfa1 = a1*b1*(1-np.exp(-b1*t))
            dfa2 = a2*b2*(1-np.exp(-b2*t))
            return dfa1+dfa2
        #d2(F1)/dt2
        def accel_F2(t,a1,b1,a2,b2):
            if t is None: t=self.xdata
            #t+=t0
            ddfa1 = a1*b1*b1*np.exp(-b1*t)
            ddfa2 = a2*b2*b2*np.exp(-b2*t)
            return ddfa1+ddfa2
# Alternative from Sanchidrian paper - F3:
        def expand_F3(t, a, b, u, t0=0):  
            if t is None: t=self.xdata 
            if not timeshift: t0=0.
            k1 = b*(t-t0)
            k2 = 1.0-np.exp(-k1)
            r = a*(t-t0)+(1.0/b)*(u-a)*k2
            return r
        #d(F3)/dt
        def veloc_F3(t, a, b, u):  
            if t is None: t=self.xdata    
            v = a*(1-np.exp(-b*t))+u*np.exp(-b*t)
            return v
        #d2(F3)/dt2
        def accel_F3(t, a, b, u):  
            if t is None: t=self.xdata 
            g1 = a*b*np.exp(-b*t)
            g2 =  u*b*np.exp(-b*t)
            return g1 - g2
        
        Fv1,Fv2,Fv3 = 1e3,1e3,1e3
        try:
            popt1, pcov1 = curve_fit(expand_F1, self.xdata, self.ydata,p0=np.array([1,1,1,0]))
            ydata1 = expand_F1(self.xdata,popt1[0],popt1[1],popt1[2],popt1[3])
            Fv1 = np.sum(np.power(ydata1-self.ydata,2))
            print " Fitting function 1 squared error = %.4f"%Fv1
        except:
            print "No function 1 evaluation"
            
        try:
            popt2, pcov2 = curve_fit(expand_F2, self.xdata, self.ydata,p0=np.array([1,1,1,1,0]))
            ydata2 = expand_F2(self.xdata,popt2[0],popt2[1],popt2[2],popt2[3],popt2[4])
            Fv2 = np.sum(np.power(ydata2-self.ydata,2))
            print " Fitting function 2 squared error = %.4f"%Fv2
        except:
            print "No function 2 evaluation"
            
        try:
            popt3, pcov3 = curve_fit(expand_F3, self.xdata, self.ydata,p0=np.array([1,1,1,0]))
            ydata3 = expand_F3(self.xdata,popt3[0],popt3[1],popt3[2],popt3[3])
            Fv3 = np.sum(np.power(ydata3-self.ydata,2))
            print " Fitting function 3 squared error = %.4f"%Fv3
        except:
            print "No function 3 evaluation"
            
            
        if (Fv1<Fv2)&(Fv1<Fv3):
            print " * Using fitting function 1 with:"
            print "\t a0 = %.4f, v_infinite= %.4f, sigma= %.4f and t0 = %.4f"%tuple(popt1)
            self.r2 = lambda t=None: expand_F1(t,popt1[0],popt1[1],popt1[2])#,popt1[3])
            self.dr2 = lambda t=None: veloc_F1(t,popt1[0],popt1[1],popt1[2])
            self.ddr2 = lambda t=None: accel_F1(t,popt1[0],popt1[1],popt1[2])
            self.xdata -= popt1[-1] # apply time shift after fitting function
        elif (Fv2<Fv1)&(Fv2<Fv3):
            print " * Using fitting function 2 with:"
            print "\t a1 = %.4f, b1= %.4f, a2= %.4f, b2= %.4f  and t0 = %.4f"%tuple(popt2)
            self.r2 = lambda t=None: expand_F2(t,popt2[0],popt2[1],popt2[2],popt2[3])#,popt2[4])
            self.dr2 = lambda t=None: veloc_F2(t,popt2[0],popt2[1],popt2[2],popt2[3])
            self.ddr2 = lambda t=None: accel_F2(t,popt2[0],popt2[1],popt2[2],popt2[3])
            self.xdata -= popt2[-1]
        else:
            print " * Using fitting function 3 with:"
            print "\t a = %.4f, b= %.4f, u= %.4f and t0 = %.4f"%tuple(popt3)
            self.r2 = lambda t=None: expand_F3(t,popt3[0],popt3[1],popt3[2])#,popt3[3])
            self.dr2 = lambda t=None: veloc_F3(t,popt3[0],popt3[1],popt3[2])
            self.ddr2 = lambda t=None: accel_F3(t,popt3[0],popt3[1],popt3[2])
            self.xdata -= popt3[-1]
        # inclination centreline angle as a function of outer expansion
        # Eq.(11) Elek et al.
        def theta(t):
            if t is None: t=self.xdata
            r2 = (r20+self.r2(t))
            dr2 = self.dr2(t)
            r2m1 = 0.5*(r20**2-r10**2)
            rc0 = 0.5*(r20+r10)
            rcterm = r2**2-r2m1
            drc = (r2*dr2)/np.sqrt(rcterm)
            return np.arctan(1e3*drc/self.props['D']) # D in [m/s] and rc in [mm/microseconds]
        # apparent acceleration:
        # Eq.(10b) Elek et al
        def a_a(t):
            if t is None: t=self.xdata
            r2 = r20+self.r2(t)
            dr2 = self.dr2(t)
            ddr2 = self.ddr2(t)
            r2m1 = 0.5*(r20**2-r10**2)
            rc0 = 0.5*(r20+r10)
            rcterm = r2**2-r2m1
            a01 = -0.5*r2*dr2*rcterm**(-3./2)
            a02 = (dr2**2)/np.sqrt(rcterm)
            a03 = r2*ddr2/np.sqrt(rcterm)
            return a01+a02+a03
        #
        self.r1 = lambda t=None: np.sqrt((self.r2(t)+r20)**2-r20**2+r10**2)-r10
        self.dr1 = lambda t=None: (self.r2(t)+r20)*self.dr2(t)/np.sqrt((self.r2(t)+r20)**2-r20**2+r10**2)
        # Eq.(12) in Elek et al.
        self.velocity = lambda t=None: 0.002*self.props['D']*np.sin(theta(t)/2.) # in [mm/microseconds]
        self.acceleration = lambda t=None: a_a(t)*np.cos(theta(t))**3 # in [mm/microseconds^2]
#
# Function to determine initial pressure and expansion ratio
        r1 = self.r1()+r10
        r2 = self.r2()+r20
        #
        #strain and strain rate functions
        strainfun = lambda t=None: 1.+(self.r1(t)-self.r2(t))/(r20-r10)
        self.strain = lambda t=None: np.log(strainfun(t))
        self.strainrate = lambda t=None : 1e6*(self.dr1(t)-self.dr2(t))/((r20-r10)*strainfun(t))
        #
        # Johnson-Cook Yield stress function: (assume no effictive elastic strain)
        self.yieldstress = lambda t=None: (JC_A+JC_B*self.strain(t)**JC_n)*(1.+JC_C*np.log(self.strainrate(t)))
        #
        # Pressure Eq.(13) in Elek et al. 2011 
        #self.pressure = 1e3*(M+C/2.)*self.acceleration()/(2.*np.pi*r1)+self.yieldstress()*(r2/r1-1)/1000.
        # Pressure Eq.(14) in Elek et al. 2015
        self.pressure = 1e3*M*self.acceleration()/(2.*np.pi*r1)+self.yieldstress()*(r2/r1-1)/1000.






    def get_expansion(self):
        ''' 
        
    This function uses the current estimated pressures to 
    determine the associated expansion ratios
    
        '''
        try:
            pressure = self.pressure
        except AttributeError:
            self.fit_r2()
            pressure = self.pressure
        # Get the properties used in the current expansion ratio calculation
        M,C,r10,D,rho0 = [self.props[kn] for kn in ['M','C','Rin','D','Rho0']]
        r1 = self.r1()+r10
        # Geometric expansion:
        ARA = (r1/r10)**2
        # expansion ratio according to Elek(2011)
        #self.expansion = ARA*(1-ARA*1e9*self.pressure/(rho0*D**2)+0.5*(M/C+0.5)*(1e3*self.velocity()/D)**2)
        # expansion ratio according to Elek(2015)
        self.expansion = ARA*(1-ARA*1e9*pressure/(rho0*D**2)+(0.5*M/C)*(1e3*self.velocity()/D)**2)






    def get_energy(self):
        '''
    
    This function calculates the various energy components
    according to Elek et al. (2015)
    
        Ekin    Equation (27)
        Wdef    Equation (29)
        E1      Equation (30)
        
    and total specific internal energy:
    
        Equation (31):
        
        energy = E0+E1-Ekin-Wdef
        
        '''
        # Here if the expansion ratio is not defined, the get_expansion function is first called
        try:
            expansion = self.expansion
        except AttributeError:
            self.get_expansion()
            expansion = self.expansion
        pressure = self.pressure
        #
        # Get the relevant property values fron the "props" dictionary:
        # cylinder properties:
        r10,r20,rhoCu,M = [self.props[kn] for kn in ['Rin','Rout','RhoCu','M']]
        # detonation properties
        rho0,D,E0,C = [self.props[kn] for kn in  ['Rho0','D','E0','C']]
        # For the kinetic / Gurney energy
        w2 = r20**2-r10**2
        # get r1 and r2 values:
        r1 = self.r1()+r10
        r2 = self.r2()+r20
        # velocity of r1:
        v1 = self.dr1()
        # Gurney energy - Eq.(27) Elek et al 2015:
        r1sq = r1**2
        self.Ekin = (0.0005*rho0*v1**2)*((M*r1sq/(C*w2))*np.log(1+w2/r1sq)+0.5) # in [GPa]
        # Deformation work:
        self.Wdef = 0.005*self.yieldstress()*((r20/r10)**2-1)*np.log((r1+r2)/(r10+r20)) # in [GPa]
        # Geometric expansion:
        ARA = (r1/r10)**2
        # particle velocity:
        u = D*(1-expansion/ARA)
        # E1 energy term:
        self.E1 = pressure*u/D - (0.5*rho0*u**2)/1e9 # in [GPa]
        # total energy:
        self.energy = E0+self.E1-self.Ekin-self.Wdef







    def fit_JWL(self,**kwargs):
        '''
        
    The method used to calibrate the Jones-Wilkins-Lee equation of state 
    for explosive products as described in the paper:
    
        P.M. Elek et. al.
        "Determination of Detonation Products Equation of State from Cylinder Test:
        Analytical Model and Numerical Analysis". Thermal Science 19(1) pp 38-45, 2015.
        
    Note that the allowable JWL parameter bounds may be given 
    as keyword arguments or a **dictionary:
    
    example:
    
        fit_JWL( A=(0,1000), B=(0,100), C=(0.001,10), R1=(1,10), ... etc )
        
        OR
        
        fit_JWL( **{'A':(0,1000),'B':(0,100),'C':(0.001,10),'R1':(1,10),...etc} )
        
    Default allowable ranges will be assigned for parameters not explicitly
    mentioned during the function call:
    
        **kwargs = {
                    'A':(0,1000),
                    'B':(0,100),
                    'C':(0.001,10),
                    'R1':(1,10),
                    'R2':(1,5),
                    'w':(1.0,1.6)
                    }
        
        '''
        # determine eveyting required if not yet determined
        try:
            expansion = self.expansion
        except AttributeError:
            self.get_expansion()
            expansion = self.expansion
        pressure = self.pressure
        #
        # default parameter bounds:
        params = {'A':(0,1000),
          'B':(0,100),
	  'C':(0.001,10),
	  'R1':(1,10),
	  'R2':(1,5),
	  'w':(1.0,1.6),
	  # Values not selected automatically - defined in the input file / assumed known
	  #'E0':(5,8), 
	  #'PCJ':(20.,30.)
	  }
        #
        for kn in kwargs.keys():
            if (kn in params.keys()):
                lbcorrect = False
                if (np.array(kwargs[kn]).size==2):
                    lb,ub = np.array(kwargs[kn])
                    if is_number(lb)&is_number(ub):
                        if lb<ub:
                            lbcorrect = True
                            params[kn]=(lb,ub)
                if not lbcorrect:
                    print "\n**Parameter "+kn+" should be descrived using a tuple of lower(lb) and upper(ub) bound values as: "+kn+"=(lb,up)"
            else:
                print "\n**Parameter <<"+kn+">> associated with the bounds prescribed in fit_elek not understood\n\tRecognised keywords are:\n\n>> A, B, C, R1, R2, w(=omega+1>1)"
        #
        # Parameter bounds used on optkeys/x_names
        optkeys = ['R1','R2','w']
        bounds = [params[kn] for kn in optkeys]
        #
        # Detonation product properties defined in the input file
        PCJ,E0,rho0,D = [self.props[kn] for kn in ['PCJ','E0','Rho0','D']]
        # Tec-note Eq (14) Capman-Jouguet Pressure:
        PCJ0 = rho0*D*D*1e-9
        VCJ = 1.-PCJ/PCJ0 # NOTE VCJ<1 
        # Given the values at the Chapman-Jouget point above, the JWL values of A, B and C can be determined given a guess on the R1, R2 and w = omega+1 values
        def getABC(R1,R2,w,printerr=True):
            # Initialise matrix and vector
            MAT,VEC = np.zeros((3,3)),np.zeros(3)
            # base terms of JWL:
            eR1,eR2,Vw = np.exp(-R1*VCJ),np.exp(-R2*VCJ),VCJ**(-w) # NOTE w = omega+1 following original DPS implementation
            # CJ-Eq1 - Elek et al Eq.(33)
            MAT[0] = [eR1,eR2,Vw]
            VEC[0] = PCJ
            # CJ-Eq2 - Elek et al Eq.(34)
            wx = w-1. # NOTE remember that w=omega+1
            MAT[1] = [eR1/R1,eR2/R2,VCJ*Vw/wx]
            VEC[1] = E0+0.5*PCJ*(1-VCJ)
            # CJ-Eq3 - Elek et al Eq.(35)
            MAT[2] = [R1*eR1,R2*eR2,w*VCJ**(-w-1)]
            VEC[2] = PCJ0
            #
            return np.linalg.solve(MAT,VEC)
        #
        ## JWL function for pressure
        def JWLfunc(V, A, B, C, R1, R2, w):
            term1 = A*np.exp(-R1*V)
            term2 = B*np.exp(-R2*V)
            TT = V**(-1.0*(1+w))
            term3 = C*TT
            p = term1 + term2 + term3
            return p
        #
        # The JWL energy equation:
        def JWLderE(V, A, B, C, R1, R2,w):
            wx=w-1.
            term1 = (A/R1)*np.exp(-R1*V)
            term2 = (B/R2)*np.exp(-R2*V)
            term3 = (C/wx)/(V**wx)
            E = term1 + term2 + term3
            return E
        #
        # objective function used during the parameter identification
        # NOTE: optkeys = [R1,R2,w]
        def objective(x,returnJWL = False):
            R1,R2,w = np.array(x)
            ## make sure that R1>3*R2
            R1t,R2t = np.max([R1,R2]),np.min([R1,R2])
            R1,R2 = R1t,R2t
            ## check R2*3<=R1 and w>1:
            if (R2*3>R1)|(w<1.):
                return 1e10
            # calculate the associated A, B and C JWL parameters
            jwlA,jwlB,jwlC = getABC(R1,R2,w,False)
            # construct the current JWL dictionary:
            jwl_vals={'A':jwlA,'B':jwlB,'C':jwlC,'R1':R1,'R2':R2,'w':w}
            if returnJWL:
                return jwl_vals
            # check if any of the JWL parameters are out of the allowed bounds:
            if any([jwl_vals[kn]<params[kn][0] or jwl_vals[kn]>params[kn][1] for kn in params.keys()]):
                return 1e10 # invalid parameter selection so return the penalty value
            #
            # Using the values A,B,C,R1,R2,w and E0:
            # evaluate the fit between the energy values:
            #     Etotal = E0+E1-Ekin-Wdef (Eq. 31)
            #            and
            #     E_JWL  = A/R1*e^(-R1*V)+B/R2*e^(-R2*V)+C/w*V^-w:
            #
            E_JWL = JWLderE(self.expansion,*tuple([jwl_vals[kn] for kn in ['A','B','C','R1','R2','w']]))
            # percentage difference:
            FV = 100*np.sum(np.abs(E_JWL/self.energy-1))
            return FV
        #
        #
        Pconverged = False
        Pcounter = 0 # maximum allowed iterations to converged pressure values = 20
        # differential evolution parameters
        diff_evol_param = {
            'maxiter':100000, ## maximum allowable iterations
            'popsize':100, ## initial population size
            'disp':False, ## should information be displayed
            'tol':0.001, 
            'strategy':'best1bin'}
        #
        convHist = []
        tol=1e-3
        Pmax = 50
        while (Pcounter<Pmax)&(not Pconverged):
            pressure = self.pressure.copy()
            self.get_expansion()
            self.get_energy()
            Pcounter +=1
            reslts = differential_evolution(objective,bounds,**diff_evol_param)
            jwl_vals = objective(reslts.x,True)
            # if the differential evolution did not find any results, display an error and return
            if jwl_vals == 1e10:
                print "*"*80+"\n*\t WARNING: the current parameter bounds could not find any viable match to the experimental data set\n\t!!! PLEASE ALTER THE ALLOWABLE PARAMETER RANGE !!!\n"+"*"*80
                yn = raw_input("\n\t>>> would you like to see help on the energy characterisation function? (Y/n)").lower()
                if 'n' in yn:
                    return
                else:
                    print self.fit_JWL.__doc__
                    return
            
            
            self.pressure = JWLfunc(self.expansion, *tuple([jwl_vals[kn] for kn in ['A','B','C','R1','R2','w']]))
            pressdiff = np.average(np.abs(self.pressure/pressure-1))
            convHist+=[pressdiff]
            print "Iteration %i Jones-Wilkins-Lee parameters"%Pcounter 
            print "\tAverage pressure difference = %.4f"%(100*pressdiff)+' % '+'vs %.4f'%(100*tol) + ' % '+ 'desired tolerance '
            print "\tAverage internal energy error = %.4f"%(reslts.fun/pressure.size)+' %'
            print "\tA = %.4f, B = %.4f, C = %.4f, R1 = %.4f, R2 = %.4f, w+1 = %.4f"%tuple([jwl_vals[kn] for kn in ['A','B','C','R1','R2','w']])
            if (pressdiff<tol):
                Pconverged = True
        
        self.jwl_vals = jwl_vals
        self.jwl_vals['Convergence'] = np.array(convHist)





        
    def JWL_pressure(self,expansion=None,**kwargs):
        '''
    The Jones-Wilkins-Lee Pressure
    
    The function may be called without expansion and property values
    
        >> JWL_pressure()
        
        in this case the self defined property values and expansion is used.
        If the self defined properties or expansion is undefined, it is determined 
        before evaluationg the JWL pressure curves
    
    It is also possible to determine the pressure at (a) specfic point(s):
        if Vdata is a scalar expansion value or array of values
        
        >> JWL_pressure(Vdata)
        
        returns the pressure values at Vdata points
    
    The function may also be called using different parameter values defined using keywords
    associated with the different JWL model keys
        JWLkeys = ['A','B','C','R1','R2','w']
    
        for example:
        
        JWL_pressure( A=370., B=2., C=1., ... etc)
        
        OR
        
        JWL_pressure( **{ 'A':370, 'B':2., 'C':1., ... etc} )
        '''
        if expansion is None:
            try:
                expansion = self.expansion
            except AttributeError:
                self.get_expansion()
                expansion = self.expansion
        # All of the JWL parameters necessary to evaluate the pressure expansion curve
        JWLkeys = ['A','B','C','R1','R2','w']
        # check if a scpecific set of parameter values are defined - alternatively use the class values themselves
        if all([kn in kwargs.keys() for kn in JWLkeys]):
            # Fully defined JWL parameters for use in current JWL evaluation
            jwl_vals = kwargs
        else:
            try:
                jwl_vals = self.jwl_vals
            except AttributeError:
                self.fit_JWL()
                jwl_vals = self.jwl_vals
        # with the JWL parameters defined, extract the pressure curve at given expansion values:
        A,B,C,R1,R2,w = [jwl_vals[kn] for kn in ['A','B','C','R1','R2','w']]
        term1 = A*np.exp(-R1*expansion)
        term2 = B*np.exp(-R2*expansion)
        TT = expansion**(-1.0*(1+w))
        term3 = C*TT
        p = term1 + term2 + term3
        return p






    def JWL_energy(self,expansion=None,**kwargs):
        '''
    The Jones-Wilkins-Lee Specific Energy
    
    The function may be called without expansion and property values
    
        >> JWL_energy()
        
        in this case the self defined property values and expansion is used.
        If the self defined properties or expansion is undefined, it is determined 
        before evaluationg the JWL energy curves
    
    It is also possible to determine the energy at (a) specfic point(s):
        if Vdata is a scalar expansion value or array of values
        
        >> JWL_energy(Vdata)
        
        returns the specific energy values at Vdata points
    
    The function may also be called using different parameter values defined using keywords
    associated with the different JWL model keys
        JWLkeys = ['A','B','C','R1','R2','w']
    
        for example:
        
        JWL_energy( A=370., B=2., C=1., ... etc)
        
        OR
        
        JWL_energy( **{ 'A':370, 'B':2., 'C':1., ... etc} )
        '''
        if expansion is None:
            try:
                expansion = self.expansion
            except AttributeError:
                self.get_expansion()
                expansion = self.expansion
        # All of the JWL parameters necessary to evaluate the pressure expansion curve
        JWLkeys = ['A','B','C','R1','R2','w']
        # check if a scpecific set of parameter values are defined - alternatively use the class values themselves
        if all([kn in kwargs.keys() for kn in JWLkeys]):
            # Fully defined JWL parameters for use in current JWL evaluation
            jwl_vals = kwargs
        else:
            try:
                jwl_vals = self.jwl_vals
            except AttributeError:
                self.fit_JWL()
                jwl_vals = self.jwl_vals
        # with the JWL parameters defined, extract the pressure curve at given expansion values:
        A,B,C,R1,R2,w = [jwl_vals[kn] for kn in ['A','B','C','R1','R2','w']]
        wx=w-1.
        term1 = (A/R1)*np.exp(-R1*expansion)
        term2 = (B/R2)*np.exp(-R2*expansion)
        term3 = (C/wx)/(expansion**wx)
        E = term1 + term2 + term3
        return E