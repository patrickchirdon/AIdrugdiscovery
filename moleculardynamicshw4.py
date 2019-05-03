?import numpy as np

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
                        
                        F = 24*eps*( ((sig**12) /(dist**12))-( (sig**6)/(dist**6)) )
                        PE = 4 * eps * ((sig/dist)**12) -((sig/dist)**6)
                        pesum+=PE
                        VIR=VIR + 24 *eps * ((2*sig**12)/(dist**12)-(sig**6/ (dist**6)))*(1/3)
                        
                        FX[i] = (F * rxdx) 
                        FY[i] = (F * rxdy) 
                        FZ[i] = (F * rxdz) 
                    
                   
                            
                    Fx[i]=FX[i]
                    Fy[i]=Fy[i]
                    Fz[i]=FZ[i]
                    
                    T_inst=(2*KE)/(3*(particle_count-1)) 
                    if(T_inst >1):
                        lamda = np.sqrt(1 + .001 *(Tdes / T_inst - 1))
                        vx[i] = np.multiply(lamda,vx[i])
                        vy[i]= np.multiply(lamda,vy[i])
                        vz[i]=np.multiply(lamda,vz[i])
                        KE =KE + 0.5 * np.sum((np.sqrt(vx[i]**2+vy[i]**2+vz[i]**2))**2);
                        
                   
                j += 1      
           
            ax[i] = ((FX[i]))/m;
            ay[i] = ((FY[i]))/m;
            az[i] = ((FZ[i]))/m;
            
                
            
            vx[i] = vx[i] + 0.5 * dt*(ax[i]);
            vy[i] = vy[i] + 0.5 * dt*(ay[i]);
            vz[i] = vz[i] + 0.5 * dt*(az[i]);
                
            
            x =  xx[i] + dt*vx[i]
            y =  yy[i] + dt*vy[i]
            z =  zz[i] + dt*vz[i]
            
            
            xx[i] = x
            yy[i] = y
            zz[i] = z
           
            

            KE +=  0.5 * m * np.sqrt(vx[i]**2+vy[i]**2+vz[i]**2);
            
            i += 1
        
        kinetic[u] += KE
        potential[u]+=pesum
        
        KE_avg=kinetic[u]/T;

        VIR_avg=VIR/T;

        Temp=(2.0*KE_avg)/(3*(particle_count-1));
        temp[u]=Temp

        P=(particle_count*Temp*eps+VIR_avg)/(particle_count *1/rho);
        pressure[u]= P
        thetotal[u]+=kinetic[u] + potential[u]
       
        
        mylist=[u,temp,pressure,kinetic,potential,thetotal, xx, yy, zz]
           
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

#forces(sig, eps, T, rho, m, dt, box_length, minimum_distance)    
results = forces(1, 1, 20, 0.5, 1, 0.01, 8, 0.8)

