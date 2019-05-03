import numpy as np

from matplotlib import pyplot as plt



def initial(box_length, minimum_distance):

    #constants
    lengthX = box_length
    lengthY = box_length
    lengthZ = box_length
    #minimum_distance = 0.8
    min_dist = box_length
    density=.5
    
    particle_count = int(density * lengthX * lengthY * lengthZ)

    xx = np.zeros(particle_count, float)
    yy = np.zeros(particle_count, float)
    zz = np.zeros(particle_count, float)

    for i in range(0, particle_count):
      
        #randomly positions particles
        x = np.random.uniform(0, lengthX)
        y = np.random.uniform(0, lengthY)
        z = np.random.uniform(0, lengthZ)

        if i > 0:
            j = 0
            while j < i:
            
                #print(j)
                
                rxdx=x-xx[j]
                rxdx =rxdx -round(rxdx/lengthX)*lengthX
                
                rxdy=y-yy[j]
                rxdy =rxdy -round(rxdy/lengthY)*lengthY
                
                rxdz=z-zz[j]
                rxdz =rxdz -round(rxdz/lengthZ)*lengthZ

                dist = np.sqrt(rxdx**2+rxdy**2+rxdz**2)
                
                j = j + 1

                if dist < minimum_distance:

                    j = 0
                    x = np.random.uniform(0, lengthX)
                    y = np.random.uniform(0, lengthY)
                    z = np.random.uniform(0, lengthZ)
                    
                else:

                    if dist < min_dist:
                        min_dist = dist

        xx[i] = x
        yy[i] = y
        zz[i] = z
        #print(i, x, y, z)

    #print('The minimum distance between the particles is', min_dist)
    return(xx, yy, zz, min_dist, particle_count)

def forces(sig, eps, T, rho, m, dt, box_length, minimum_distance):
    lengthX = box_length
    lengthY = box_length
    lengthZ = box_length
    xx = initial(lengthX, minimum_distance)[0]
    yy = initial(lengthY, minimum_distance)[1]
    zz = initial(lengthZ, minimum_distance)[2]
    x=xx
    y=yy
    z=zz
    min_dist = initial(box_length, minimum_distance)[3]
    particle_count = initial(box_length, minimum_distance)[4]
    u = 0
    FX = np.zeros(particle_count)
    FY = np.zeros(particle_count)
    FZ = np.zeros(particle_count)
    fx=np.zeros(particle_count)
    fy=np.zeros(particle_count)
    fz=np.zeros(particle_count)
    Fx=np.zeros(particle_count)
    Fy=np.zeros(particle_count)
    Fz=np.zeros(particle_count)
    ax = np.zeros(particle_count)
    ay = np.zeros(particle_count)
    az = np.zeros(particle_count)
    axold = np.zeros(particle_count)
    ayold = np.zeros(particle_count)
    azold = np.zeros(particle_count)
    vx = np.zeros(particle_count)
    vy = np.zeros(particle_count)
    vz = np.zeros(particle_count)
    potential = np.zeros(T)
    kinetic = np.zeros(T)
    lamda=1
    temp=np.zeros(particle_count)
    pressure=np.zeros(particle_count)
    thetotal=np.zeros(particle_count)
    theradii=np.zeros(particle_count)
    theg=np.zeros(particle_count)
    VIR=0

    T_inst=0

    Tdes=3

    KE_avg=0

   
    
    while(u < T):
        PE = 0
        KE = 0
        i = 0
        while(i < particle_count):
            j = 0
            FX[i] = 0
            FY[i] = 0
            FZ[i] = 0
            pesum=0
            p=1
            while(j < particle_count):
                if (i != j):
                    rxdx = xx[i]-xx[j]
                    rxdx = rxdx -round(rxdx/lengthX)*lengthX
                    rxdy = yy[i]-yy[j]
                    rxdy = rxdy -round(rxdy/lengthY)*lengthY
                    rxdz = zz[i]-zz[j]
                    rxdz = rxdz -round(rxdz/lengthZ)*lengthZ
                    dist = np.sqrt(rxdx**2+rxdy**2+rxdz**2)
                    if(dist <= 3*sig):
                        
                        
                        PE = 4 * eps * ((sig/dist)**12) -((sig/dist)**6)
                        pesum+=PE
                        VIR=4 *eps * ((sig**12)/(dist**12)-(sig**6)/ (dist**6))
                        
                       
                    
                   
                   
                j += 1      
           
           
            
           
            i += 1
     
       
        
        mylist=[pesum, VIR, xx, yy, zz]
           
        u += 1
    
    return(mylist)      


def pairCorrelationFunction_3D(x, y, z, S, rMax, dr):

    #find how particles are distributed
    #if there is a peak this means there are many particles at that location 
    #counts the number of particles within a certain radius 
    #radii is the distance
    #g is the number of particles
    from numpy import zeros, sqrt, where, pi, mean, arange, histogram, absolute
    
    r = np.zeros([256,3])  #create a matrix for x,y,z one entry for each of the 256 molecules
    r[:,0]=x   #set the first column equal to x
    r[:,1]=y   #set the second column equal to y
    r[:,2]=z   #set the third column equal to z
    num_particles  = len(r)  
    rMax  = S/2.0;
    edges          = arange(0., rMax + dr, dr)
    num_increments = len(edges) - 1
    g              = zeros(num_increments)
    radii          = zeros(num_increments)
    numberDensity  = len(r) / S**3

    # Compute pairwise correlation for each particle
    for index in range(num_particles):

        d = 0.0

        for i in range(3):
            dp = absolute(r[index,i] - r[:,i])
            mask = dp>S/2.0
            dp[mask] = S - dp[mask]
            d += dp*dp

        d = sqrt(d)
        d[index] = 2 * rMax

        (result, bins) = histogram(d, bins=edges, normed=False)
        g += result

    g = g/(num_particles * numberDensity)

    # Average g(r) divide by apropriate radii
    for i in range(num_increments):
        radii[i] = (edges[i] + edges[i+1]) / 2.
        rOuter = edges[i + 1]
        rInner = edges[i]
        g[i] = g[i] / (4.0 / 3.0 * pi * (rOuter**3 - rInner**3))
    
    
    return (radii, g)

density= .5
X= 8
Y= 8
Z= 8
Eps=1
Sig=1
M=1
T=100
    
#forces(sig, eps, T, rho, m, dt, box_length, minimum_distance)    
results = forces(1, 1, 1, 0.5, 1, 0.01, 8, 0.8)

Temp=1
Pres=3
V=X*Y*Z 
N=256       #Number of particles
     #Put N particle in box in a random way which each particle have minimum distance with other particles not have collision with each other.
RX= results[2]
RY=results[3]
RZ= results[4]
delvmax=5.0
ULJO=0
ULJNew=0
beta=1.0
from math import exp, log
import random
import math
delmax=.05
MC_cycle=100
ULJN=0
RXn=[0 for col in range(N)]
RYn=[0 for col in range(N)]
RZn=[0 for col in range(N)]
rho_ave=0
V_ave=0
ULJOO=0
sumsix=0
sumtwelve=0
RXnew=0
RYnew=0
RZnew=0
Vir_ave=0
deld=0.05*Sig
Nmax=int((X/(2*deld))+1)
hist=[0 for col in range(Nmax)]

for p in range(MC_cycle):
    
    for j in range(N):
        RXnew=RX[j]+ (random.random()*2-1)*delmax
        RYnew=RY[j]+ (random.random()*2-1)*delmax
        RZnew=RZ[j]+ (random.random()*2-1)*delmax
        for k in range(N):
            if k!=j:
                Dx=RX[k]-RX[j]
                Dx=Dx-round(Dx/X)*X
                Dy=RY[k]-RY[j]
                Dy=Dy-round(Dy/Y)*Y
                Dz=RZ[k]-RZ[j]
                Dz=Dz-round(Dz/Z)*Z
                Dis=Dx**2+Dy**2+Dz**2
                Dis=Dis**0.5
                Sigsix=Sig/Dis
                Sigsixsix=(Sigsix)**6
                Sigtwelve=(Sigsix)**12
                ULJO=ULJO+4*Eps*(Sigtwelve-Sigsixsix)
                Dxx=RX[k]-RXnew
                Dxx=Dxx-round(Dxx/X)*X
                Dyy=RY[k]-RYnew
                Dyy=Dyy-round(Dyy/Y)*Y
                Dzz=RZ[k]-RZnew
                Dzz=Dzz-round(Dzz/Z)*Z
                Diss=Dxx**2+Dyy**2+Dzz**2
                Dist=Diss**0.5
                Sigsixx=Sig/Dist
                Sigsixsixx=(Sigsixx)**6
                Sigtwelvee=(Sigsixx)**12
                ULJNew=ULJNew+4*Eps*(Sigtwelvee-Sigsixsixx)
        if ULJNew<ULJO:
            RX[j]=RXnew
            RY[j]=RYnew
            RZ[j]=RZnew
        else:
            Prob=math.exp(-beta*(ULJNew-ULJO))
            R=random.random()
            if R<Prob:
                RX[j]=RXnew
                RY[j]=RYnew
                RZ[j]=RZnew
                
        
        
        ULJNew=0
        ULJO=0
    if random.random()< 1:
        Vnew=V+random.random()*2*delvmax-delvmax
        delV=Vnew-V
        delU1= results[0]
        vir=results[1]
        virn=results[1]
        
        Prob1=math.exp(-beta*(delU1+Pres*delV-N*(1/beta)*math.log(Vnew/V)))
        ratio=(Vnew/V)
      
        if Prob1>1:
            for i in range(N):
                RX[i]=RX[i]*ratio
                RY[i]=RY[i]*ratio
                RZ[i]=RZ[i]*ratio
            V=Vnew
            X=X*ratio
            Y=Y*ratio
            Z=Z*ratio
            vir=virn
        else:
            R=random.random()
            if R<Prob1:
                for i in range(N):
                    RX[i]=RX[i]*ratio
                    RY[i]=RY[i]*ratio
                    RZ[i]=RZ[i]*ratio
                V=Vnew
                X=X*ratio
                Y=Y*ratio
                Z=Z*ratio
                vir=virn
    if p>50:
        V_ave=V_ave+V
        rho_ave=rho_ave+N/V
        Vir_ave=Vir_ave+vir
        for i in range(N-1):
            for j in range(i+1,N):
                Dx=RX[i]-RX[j]
                Dx=Dx-round(Dx/X)*X
                Dy=RY[i]-RY[j]
                Dy=Dy-round(Dy/Y)*Y
                Dz=RZ[i]-RZ[j]
                Dz=Dz-round(Dz/Z)*Z
                Dis= Dx**2+Dy**2+Dz**2
                Dis=Dis**0.5
                
                if Dis<X/2:

                    KK=int(Dis/deld)+1
                    hist[KK]=hist[KK]+2

    for i in range(N):
        if RX[i]>X:
            RX[i]=RX[i]-X
        if RX[i]<0:
            RX[i]=RX[i]+X
        if RY[i]>Y:
            RY[i]=RY[i]-Y
        if RY[i]<0:
            RY[i]=RY[i]+Y
        if RZ[i]>Z:
            RZ[i]=RZ[i]-Z
        if RZ[i]<0:
            RZ[i]=RZ[i]+Z
            
V_ave1=V_ave/100

rho=N/V_ave1
Vir_ave= Vir_ave/100

Paverage=Pres*((N*Temp)/V_ave1)+(Vir_ave/V_ave1)
print("rho", rho)
print("pressure", Paverage)

r, g = pairCorrelationFunction_3D(RX, RY, RZ, 8, 4, .1)
