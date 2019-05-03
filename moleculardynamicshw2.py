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
                    
                        
                        FX[i] = (F * rxdx) 
                        FY[i] = (F * rxdy) 
                        FZ[i] = (F * rxdz) 
                    
                   
                            
                    Fx[i]=FX[i]
                    Fy[i]=Fy[i]
                    Fz[i]=FZ[i]
                   
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
        
        print(kinetic[u], potential[u])
        
        u += 1
    
    return(kinetic, potential)      

    
#forces(sig, eps, T, rho, m, dt, box_length, minimum_distance)    
results = forces(1, 1, 20, 0.5, 1, 0.01, 8, 0.8)
print('Kinetic energy: ', results[0])
print('Potential energy: ', results[1])
print('total energy:', results[0] + results[1])

plt.subplot(3, 1, 1)
plt.plot(results[0], 'o-')
plt.title('Energy vs. Time')
plt.ylabel('kinetic energy')

plt.subplot(3, 1, 2)
plt.plot(results[1], '.-')
plt.xlabel('time (s)')
plt.ylabel('potential energy')

plt.subplot(3, 1, 3)
plt.plot(results[1] + results[0], '.-')
plt.xlabel('time (s)')
plt.ylabel('total energy')
plt.show()
    
