#
#
#   THIS IS AN EXAMPLE ILLUSTRATING HOW TO USE THE <<jwl>> CLASS CODED IN "cylindertest.py"
#
#
from cylindertest import jwl

# All cylinder / detonation / expansion information now contained in the text file
nameof_file = './LLNL_TNT.txt'

# Read input file intto a jwl class named "Test"
Test = jwl(nameof_file)

# May also assign different Johnson-Cook properties to the cylinder
# Elek et al only use the yield stress value so A= 89.63 B=0. C=0. and n=0 in the text file
JCookV = {}
# NOTE: Uncomment next line to use A = 89.63, B= 291.6, C= 0.025 and n= 0.31 instead
#JCookV = {'JC_A':89.63,'JC_B':291.6,'JC_C':0.025,'JC_n':0.31}
for kn in JCookV.keys():
    Test.props[kn] = JCookV[kn]
    
    
# Check the different r2 function formulations and apply the best without allowing an initial time shift >> t0=0 <<
Test.fit_r2()
# NOTE: run the same function using Test.fit_r2(True) to allow an initial time shift, i.e. t0<>0
#Test.fit_r2(True)


# a different function calculates the energy terms and total specific energy
Test.get_energy()


# NOTE: it is possible to copy pressure, expansion and energy to P1,V1 and E1 arrays before fitting / modification for example if necessary:
P1 = Test.pressure.copy()
V1 = Test.expansion.copy()
E1 = Test.energy.copy()


# Fit the JWL parameters according to the method described by Elek et al. (2015)
Test.fit_JWL()


# Copy the final values of the pressures, expansions and different energy contributions:
P2 = Test.pressure.copy()
V2 = Test.expansion.copy()
E2 = Test.energy.copy()
# use the parameters found to evaluate the JWL_pressure and JWL_energy functions
P_jwl = Test.JWL_pressure()
E_jwl = Test.JWL_energy()


# NOTE: it is also possible to evaluate the JWL pressure and energy curves given a different set of parameters, defined in a dictionary or set of keyword arguments:

# The TNT parameters according to Elek et al.
ElekJWL = {
    'A':366.42,
    'B':2.6983,
    'C':1.1480,
    'R1':4.1245,
    'R2':0.9436,
    'w':1.3135,# w=omega+1
    }
# evaluated at the same expansion points as the test itself
P_elek = Test.JWL_pressure(**ElekJWL)
E_elek = Test.JWL_energy(**ElekJWL)

# The TNT parameter according to Dobratz and Crawford
LLNL_JWL = {
    'A':371.21,
    'B':3.2306,
    'C':1.10453,
    'R1':4.15,
    'R2':0.95,
    'w':1.3,# w=omega+1
    }
# evaluated at the same expansion points as the test itself
P_llnl = Test.JWL_pressure(**LLNL_JWL)
E_llnl = Test.JWL_energy(**LLNL_JWL)







# plotting
from pylab import ion,close,figure,show
#ion()
#close('all')
#
# print the data compare to the fitting function:
fig = figure(1)
ax = fig.add_subplot(111)
ax.plot(Test.xdata,Test.ydata,'k.',label=r'$\mathrm{Data}$')
ax.plot(Test.xdata,Test.r2(),'r',label=r'$\Delta r_2$')
ax.legend(loc='upper left')
ax.set_ylabel(r'$\Delta r_2\;\mathrm{[mm]}$',fontsize=16)
ax.set_xlabel(r'$\mathrm{t\;[\mu s]}$',fontsize=16)



fig = figure(2)
ax = fig.add_subplot(111)
ax.plot(Test.xdata,Test.r2()+Test.props['Rout'],'r',label=r'$r_2$')
ax.plot(Test.xdata,Test.r1()+Test.props['Rin'],'b',label=r'$r_1$')
ax.legend(loc='upper left')
ax.set_ylabel(r'$\mathrm{Radii\;[mm]}$',fontsize=16)
ax.set_xlabel(r'$\mathrm{t\;[\mu s]}$',fontsize=16)



# Expansion ratios
fig = figure(3)
ax = fig.add_subplot(111)
ax.plot(Test.xdata,((Test.r1()+Test.props['Rin'])/Test.props['Rin'])**2,'r',label=r'$\mathrm{Geometric\;expansion}$')
ax.plot(Test.xdata,Test.expansion,'b',label=r'$\mathrm{Calculated\;expansion}$')
ax.legend(loc='upper left')
ax.set_ylabel(r'$\mathrm{Expansion\;ratio}$',fontsize=16)
ax.set_xlabel(r'$\mathrm{t\;[\mu s]}$',fontsize=16)



# Strain and strain rate
fig = figure(4)
ax1 = fig.add_subplot(111)
ax2 = ax1.twinx()
ax1.plot(Test.xdata,Test.strain(),'k')
ax2.plot(Test.xdata,Test.strainrate()/1e4,'b')
ax1.set_xlabel(r'$\mathrm{t\;[\mu s]}$',fontsize=16)
ax1.set_ylabel(r'$\mathrm{Strain}$',fontsize=16)
ax2.set_ylabel(r'$\mathrm{Strain\;rate\;[10^4\;s^{-1}]}$',fontsize=16)
ax2.yaxis.label.set_color('blue')
ax2.tick_params(axis='y',colors='b')



# Stress-strain curve
fig = figure(5)
ax= fig.add_subplot(111)
ax.plot(Test.strain(),Test.yieldstress(),'k')
ax.set_xlabel(r'$\mathrm{Strain}$',fontsize=16)
ax.set_ylabel(r'$\mathrm{Yield\;stress\;[MPa]}$',fontsize=16)



# Pressure Expansion data
fig = figure(6)
ax = fig.add_subplot(111)
#ax.plot(V1,P1,'y.',label=r'$\mathrm{Initial\;Pressure}$')
ax.plot(V2,P2,'k.',label=r'$\mathrm{Final\;Pressure}$')
ax.plot(V2,P_jwl,'r',label=r'$\mathrm{JWL\;Pressure}$')
ax.plot(V2,P_elek,'b',label=r'$\mathrm{Elek\;et\;al.(2015)}$')
ax.plot(V2,P_llnl,'g',label=r'$\mathrm{Dobratz\;and\;Crawford\;(1985)}$')
ax.set_yscale('log')
ax.legend(loc='upper right')
ax.set_xlabel(r'$\mathrm{Expansion}$',fontsize=16)
ax.set_ylabel(r'$\mathrm{Pressure\;[GPa]}$',fontsize=16)



# Energy data
fig = figure(7)
ax = fig.add_subplot(111)
#ax.plot(V1,E1,'y.',label=r'$\mathrm{Initial\;Energy}$')
ax.plot(V2,E2,'k.',label=r'$\mathrm{Final\;Energy}$')
ax.plot(V2,E_jwl,'r',label=r'$\mathrm{JWL\;Energy}$')
ax.plot(V2,E_elek,'b',label=r'$\mathrm{Elek\;et\;al.(2015)}$')
ax.plot(V2,E_llnl,'g',label=r'$\mathrm{Dobratz\;and\;Crawford\;(1985)}$')
ax.set_yscale('log')
ax.legend(loc='upper right')
ax.set_ylabel(r'$\mathrm{Specific\;Energy\;[GPa]}$',fontsize=16)
ax.set_xlabel(r'$\mathrm{Expansion\;ratio,\; V}$',fontsize=16)


# Energy contributions
fig = figure(8)
ax = fig.add_subplot(111)
Vmm = [Test.expansion[0],Test.expansion[-1]]
E0 = Test.props['E0']
ax.plot(Vmm,[E0,E0],'c',label=r'$E_0$')
ax.plot(Test.expansion,Test.Ekin,'g',label=r'$E_\mathrm{kin}$')
ax.plot(Test.expansion,Test.Wdef,'r',label=r'$W_\mathrm{def}$')
ax.plot(Test.expansion,Test.E1,'m',label=r'$E_1$')
ax.plot(Test.expansion,Test.energy,'b',label=r'$E_\mathrm{Total}$')
ax.legend(loc='upper right')
ax.set_ylabel(r'$\mathrm{Specific\;Energy\;[GPa]}$',fontsize=16)
ax.set_xlabel(r'$\mathrm{Expansion\;ratio,\; V}$',fontsize=16)
    
show()