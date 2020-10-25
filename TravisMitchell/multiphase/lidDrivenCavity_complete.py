"""
Manipulated from VoF tutorial of Prof. G. Tryggvason and 
originally coded into Python by the PyCFD team. Edits made
by TM to assist understanding for MECH6480 at UQ, 2020.

This was a tutorial that I worked through at the start
of my PhD to understand multiphase flow methodlogy.

"""
__author__  = 'Travis Mitchell'
__version__ = '1.0.0'
__date__    = '10/02/2016'

from matplotlib import cm
import matplotlib.pyplot as plt
import numpy as np

######## SIMULATION PARAMETERS ########
#Domain size and body forces (i.e. gravity)
Lx=1.0
Ly=1.0
gx=0.0
gy=0.0

#Fluid properties
rho1=1.0
m0=0.05
#Velocity boundary conditions
unorth=2.
usouth=0.
veast=0.
vwest=0.

#Discretisation
nx=41
ny=41
dt=0.0004
nstep=150
maxiter=1000    #Maximum iterations in pressure equation
maxError=0.001  #Acceptable error tolerance in pressure equation
beta=1.2        #SOR relaxation variable in pressure equation

#Flow variable arrays
u=np.zeros((nx+1,ny+2)) #horizontal velocity
v=np.zeros((nx+2,ny+1)) #vertical velocity
p=np.zeros((nx+2,ny+2)) #pressure field

ut=np.zeros((nx+1,ny+2)) #non-divergent-free horizontal velocity
vt=np.zeros((nx+2,ny+1)) #non-divergent-from vertical velocity

#temp variables for pressure calculations
tmp1=np.zeros((nx+2,ny+2))
tmp2=np.zeros((nx+2,ny+2))

#velocities at center of grid for plotting
uu=np.zeros((nx+1,ny+1))
vv=np.zeros((nx+1,ny+1))
x = np.linspace(0.,Lx,nx)
y = np.linspace(0.,Ly,ny)

#Define the grid
dx=Lx/(nx) 
dy=Ly/(ny)

r=np.ones((nx+2,ny+2))*rho1 #density field (defind here for when we extend to multiphase)
            
time=0.
for steps in range(nstep):    
    #tangential velocity at boundaries
    u[:,0]=2.*usouth-u[:,1]
    u[:,-1]=2.*unorth-u[:,-2]
    v[0,:]=2.*vwest-v[1,:]
    v[-1,:]=2.*veast-v[-2,:]
    
    # Predict u-velocity                               
    for i in range(1,nx):
        for j in range(1,ny+1):
            #Slides Eqn 15
            Ax = -0.25*(((u[i+1,j]+u[i,j])**2-(u[i,j]+\
                          u[i-1,j])**2)/dx+((u[i,j+1]+u[i,j])*(v[i+1,j]+\
                          v[i,j])-(u[i,j]+u[i,j-1])*(v[i+1,j-1]+v[i,j-1]))/dy)
            #Slides Eqn 16
            Dx = m0*((u[i+1,j]-2*u[i,j]+u[i-1,j])/dx**2+\
                (u[i,j+1]-2*u[i,j]+u[i,j-1])/dy**2 )
            ut[i,j]=u[i,j]+dt*(Ax + gx + (2.0/(r[i+1,j] + r[i,j])) * Dx)


    # Predict v-velocity                               
    for i in range(1,nx+1):
        for j in range(1,ny):
            Ay = -0.25*(((u[i,j+1]+u[i,j])*(v[i+1,j]+\
                v[i,j])-(u[i-1,j+1]+u[i-1,j])*(v[i,j]+v[i-1,j]))/dx+\
                ((v[i,j+1]+v[i,j])**2-(v[i,j]+v[i,j-1])**2)/dy)
            Dy = m0*((v[i+1,j]-2*v[i,j]+v[i-1,j])/dx**2+\
                (v[i,j+1]-2*v[i,j]+v[i,j-1])/dy**2 )    
            vt[i,j]=v[i,j]+dt*(Ay + gy + (2.0/(r[i,j+1] + r[i,j])) * Dy)
    
    #Compute source term and the coefficient for p(i,j)
    #Here I use the Successive Over Relaxation (SOR) method
    #You can also try to solve this with e.g. Jacobi method 
    #if you are more familiary.
    #SOR takes a weighted average of the new pressure and the old pressure.
    rt=r.copy()
    lrg=1000.       # This is a sneaky method of setting the pressure in the ghost cell ~0 so we don't need a separate eqn for points near the boundary
    rt[:,0]=lrg
    rt[:,-1]=lrg
    rt[0,:]=lrg
    rt[-1,:]=lrg
    
    for i in range(1,nx+1):
        for j in range(1,ny+1):
            #Slides: LHS of Eqn 17
            tmp1[i,j]= (0.5/dt)*( (ut[i,j]-ut[i-1,j])/dx+(vt[i,j]-vt[i,j-1])/dy)
            #Use to calculate new value of p_i,j. In eqn 17 of the slides
            # imagine solving the equation so you have p_i,j = ...
            tmp2[i,j]=1.0/( (1./dx)*( 1./(dx*(rt[i+1,j]+rt[i,j]))+\
                1./(dx*(rt[i-1,j]+rt[i,j])))+\
                (1./dy)*(1./(dy*(rt[i,j+1]+rt[i,j]))+\
                1./(dy*(rt[i,j-1]+rt[i,j])))) 

    iter=0
    while True:
        pn=p.copy()
        iter=iter+1
        for i in range(1,nx+1):
            for j in range(1,ny+1):
                p[i,j]=(1.0-beta)*p[i,j] + beta*tmp2[i,j]*(\
                    (1./dx**2)*( p[i+1,j]/(rt[i+1,j]+rt[i,j])+\
                    p[i-1,j]/(rt[i-1,j]+rt[i,j]))+\
                    (1./dy**2)*( p[i,j+1]/((rt[i,j+1]+rt[i,j]))+\
                    p[i,j-1]/((rt[i,j-1]+rt[i,j])))-tmp1[i,j])
        if np.abs(pn-p).max()<maxError:
            break
        if iter>maxiter:
            break
    
    #CORRECT THE u-velocity
    for i in range(1,nx):
        for j in range(1,ny+1):
            #Slides: Eqn 18
            u[i,j]=ut[i,j]-dt*(2.0/dx)*(p[i+1,j]-p[i,j])/(r[i+1,j]+r[i,j])
    
    #CORRECT THE v-velocity
    for i in range(1,nx+1):
        for j in range(1,ny):
            #Slides: Eqn 19
            v[i,j]=vt[i,j]-dt*(2.0/dy)*(p[i,j+1]-p[i,j])/(r[i,j+1]+r[i,j])

    time=time+dt
    
    #Plot the velocity field - so interpolate velocity to the nodes.
    uu=0.5*(u[0:nx,1:ny+1]+u[0:nx,0:ny])
    vv=0.5*(v[1:nx+1,0:ny]+v[0:nx,0:ny])
    yy,xx=np.mgrid[0:(nx-1)*dx:nx*1j,0:(ny-1)*dx:ny*1j]
    plt.clf()
    plt.quiver(xx,yy,uu.T,vv.T)
    plt.pause(0.01)
    print steps

plt.figure()
plt.streamplot(x,y,uu.T,vv.T)
plt.contourf(x,y,np.sqrt(uu.T**2 + vv.T**2))
plt.show()
