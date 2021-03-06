#'samuel'

# simple integrator for sine-G
# Solves u_tt - u_xx = sin(u)
# Plots the field profile, energy density and lump trajectories

# last update: 15/1/15

# This program solves the sine-Gordon equation for initial conditions
# specified by the functions u0(x) and udot0(x) defined below. At the
# moment it is set up so that the initial field is static - it will be
# part of the third assignment to modify the program to make them
# move.

from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

plt.ion()  # turn off interaction mode for plotting

#######################################################

# Range of time, and number of time steps
tmin=0.0
tmax=200.0
TSTEPS=4000
dt=(tmax-tmin)/TSTEPS

# X range of the plots
xmin=-20
xmax=20

# Y ranges of the plots

yminE=-1.0
ymaxE=5.0
yminu=-15.0
ymaxu=5
ymin=-8.5
ymax=8.0

# Number of points:
M=512

# x-axis points for plots: (note p.b.c. equates xmin to xmax)
x = np.linspace(xmin,xmax,M)

# x0 and x1 are the initial positions of the kink/anti-kink pair:
x0 = -6.5
x1 = -x0

v0 = 0
v1 = -0.5

# Display only every Kth plot (to speed up the program!)
K=TSTEPS/200

# periodicity
L=(xmax-xmin)
h=L/(M-1) # x step size

# arrays for positions of solitons versus time
sol_pos_times=[]
sol_pos_dat1=[]
sol_pos_dat2=[]

def sign(x):
	if x>0: return +1
	else: return -1

def g(v):
    return 1/(1-v**2)**(1/2)


###############################################################
# Initial form of u and udot.
def u0(x):
    return np.tanh(g(v0)*(x-x0))+np.tanh(-g(v1)*(x-x1))-1

def udot0(x):
    t=0
    return ((1/(np.cosh(g(v0)*(x-v0*t-x0)))**2)*(-g(v0)*v0))+\
           ((1/(np.cosh(-g(v1)*(x-v1*t-x1)))**2)*(g(v1)*v1))

###############################################################

# Set up initial configuration at t=0
u=[]
ux=[0.0]*M
utt=[0.0]*M
udot=[]
for i in range(M):
    xx=xmin+i*h
    u.append(u0(xx))
    udot.append(udot0(xx))
u=np.array(u)
ux=np.array(ux)
utt=np.array(utt)

# Initialise plot:
plt.clf()
plt.subplot(411)
lineu, = plt.plot(x,u,color='b')
lineu.axes.set_ylim(yminu,ymaxu)
#lineu.axes.set_yticks((-8,-4,0,4,8))
plt.title('Field profile',loc='left',fontsize=13)
plt.subplot(412)
linee, = plt.plot(x,u,color='r')
linee.axes.set_ylim(yminE,ymaxE)
# linee.axes.set_yticks((0.0,1.0,2.0,3.0))
plt.title('Energy density',loc='left',fontsize=13)
plt.subplots_adjust(hspace=0.7,bottom=0.05)

#############################
#    INTEGRATION ROUTINE    #
#############################

t=0.0
txt=plt.suptitle('t='+ '%.2f' % t,fontsize=14, fontweight='bold')

timestamp = 20
pot = []
kin = []
tot = []
time  = []

c=0
OK=1
stopearly=' '
print('starting...')
try:
 for t in np.arange(tmin,tmax+dt,dt):

# Simple integration routine:
    utt=(np.roll(u,-1)+np.roll(u,1)-2*u)/h**2-(2*u*((u**2)-1))

    udot+=dt*utt
    u+=udot*dt

# Plot every Kth plot
    if c%K==0:

       txt.set_text('Time='+ '%.2f' % t)

       # update the wave profile:
       lineu.set_ydata(u)

       # update the energy density:
       ux=(np.roll(u,-1)-u)/h
       linee.set_ydata((udot**2+ux**2)/2+((u**2-1)**2)/2)


       # redraw the plot:
       plt.draw()

       # store soliton positions so can plot at end:

       i=0
       signu = sign(u[i])
       while (i<M) and (sign(u[i]) == signu): i = i + 1
       sol1=i
       if i!=M:
         signu = sign(u[i])
       while (i<M) and (sign(u[i]) == signu): i = i + 1
       sol2=i
       if sol2!=M: #only append if have found both solitons
         sol_pos_times.append(t)
         sol_pos_dat1.append(xmin+sol1*h)
         sol_pos_dat2.append(xmin+sol2*h)

        #sample to calculate energies
    if t == 30 or t == 40  or t == 50 or t == 70:
        pe = 1/2*(ux**2)+(1/2)*((u**2)-1)**2
        ke = (udot**2)/2
        print('At time=',t,'Potential Energy:',sum(pe),'Kinetic Energy:'
            ,sum(ke),'Total Energy', sum(pe+ke))



    #appends the potential,kinetic and total energy lists
    if c%K == 0:
        ke = (udot**2)/2
        pe = 1/2*(ux**2)+(1/2)*((u**2)-1)**2
        time.append(t)
        kin.append(sum(ke))
        pot.append(sum(pe))
        tot.append(sum(pe+ke))

    c+=1

# Has it blown up?
    um=list(abs(u[:]))
    um.sort()
    if not um[-1]<1000:
        print 'Unstable... decrease the step size!'
        OK=0
    if not OK: break

except KeyboardInterrupt:
 stopearly=' early'
print('stopped'+stopearly)

plt.ioff()  # turn off interaction mode for plotting

plt.subplot(4,1,(3,4))
plt.title('Positions of the solitons versus time',loc='left',fontsize=13)
plt.xlim((xmin,xmax))
plt.ylim((0,tmax))
plt.plot(sol_pos_dat1,sol_pos_times,marker='.',linestyle='none',color='c')
plt.plot(sol_pos_dat2,sol_pos_times,marker='.',linestyle='none',color='c')
plt.draw()
plt.show()


#plots the energies with respect to time
plt.subplot(4,1,(1,3))
plt.title('Energies at any point in time',loc='left',fontsize=13)
plt.xlim((tmin,tmax))
maxsize = []
maxsize.append(max(pot))
maxsize.append(max(kin))
maxsize.append(max(tot))
plt.ylim(-10,max(maxsize)+30)
plt.plot(time,pot,color='g', label='Potential Energy')
plt.plot(time,kin,color='b', label='Kinetic Energy')
plt.plot(time,tot,color='r', label='Total Energy')
plt.legend(loc='upper right')
plt.ylabel('Energy')
plt.xlabel('Time')
plt.draw()
plt.show()

