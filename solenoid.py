from __future__ import division
from visual import *
from visual.graph import *	    #import graphing features


#Computer Resolution
widthPC = 1280
heightPC = 800

#Scene Constants
scene.width = widthPC
scene.height = heightPC


########################
###EDITABLE CONSTANTS###
########################

#Editable Cylinder Constants
cylR = 1                            #radius of cylinder
R = 15                              #radius of solenoid
L = 10*R                            #length of solenoid, in terms of R
N = 50                              #number of rotations
n = 6*4                             #number of cylinders per rotation, in multiples of 6 (for color purposes)

#Editable Vector Constants
vecPlaneTotal = 5                   #number of circle vector planes (>= 2)
v = 5                               #number of vectors per rotation (>= 1)
vecFrac = L/4                       #distance between one end of solenoid and last plane of vectors created, in terms of L

numCircles = 3                      #number of vector circles (R1, R2, ..., RnumCircles)
R1 = R*5/13                         #radius of first circle of vectors, in terms of R
R2 = R*10/13                        #radius of second circle of vectors, in terms of R
R3 = R*15/13                        #radius of three circle of vectors, in terms of R

#Editable Plot Constants
P = 50                              #number of points plotted (>= 2) in graph of Bmag vs. Z

#Editable Scale Factors
kB = .6                             #Bvector axis scale factor
kfPlot = 1000                       #force plot scale factor
kax = 3500                          #force (Farrow axis) scale factor


###############
###CONSTANTS###
###############

#Scale Factors
k = 1e5                             #mu scale factor
kq = 1e16                           #charge scale factor

#Physical Constants
I = 500                             #current in amps
mu = 1e-7                           #equal to mu0/(4*pi)
m = 1	#1.67e-27                   #test charge mass
q = 1.6e-19                         #test charge
deltat = 1                          #time interval

#Cylinder Constants
cylTotal = n*N                  #total number of cylinders in solenoid
RadTotal = N*2*pi-.001          #the "...-.001" subtraction is to accomodate for rounding errors with pi
theta = 0                       #initial theta
dtheta = 2*pi/n                 #theta-increment, in terms of n
z = -L/2                        #initial z-position
dz = L/cylTotal                 #z-increment, in terms of L
phi = 90 - dtheta/2             #angle between the line segment from the center of a rotation to the position of one cylinder
                                    #and the line segment from that cylinder's position to an adjacent cylinder's position
d = R*sin(dtheta)/sin(phi)      #distance between the positions of two adjacent cylinders found using law of sines;
                                    #used to define cylinder axis length
h = sqrt(R**2+d**2/4)           #length involved in compensating for fact that cyl pos is at end instead of center of cyl
alpha = asin(d/(2*h))           #angle involved in compensating for fact that cyl pos is at end instead of center of cyl

#Cylinder Color Constants
dColor = 6/cylTotal
upColorList = arange(0, 1, dColor)
downColorList = upColorList + dColor
upCountcyl = 0
downCountcyl = len(downColorList) - 1

#Vector Constants
vecLim = L/2 + vecFrac                          #z-limit 
Vz = -vecLim                                    #initial z-position of arrows
dVz = 2*vecLim/(vecPlaneTotal-1)                #z-increment for creating/displaying arrows
Pz = -vecLim	                                #initial z-position of plotting points
dPz = 2*vecLim/(P-1)                            #z-increment for plotting points
dvecAngle = 2*pi/v                              #angle-increment
vecTotal = vecPlaneTotal*(1+numCircles*v)       #total number of vectors


###########
###LISTS###
###########

#Cylinder List:
cylList = []                        #creates empty list to be filled with cylinder objects in Cylinder Loop

#Vector Lists:                      #vecList objects contain position information
vecListPlot = []                    #creates empty list to be filled with vectors in VecPlot Loop
BListLine = []                    #creates empty list to be filled with vectors in VecLine Loop
BListR1 = []                      #creates empty list to be filled with R1 vectors in VecCircle Loop
BListR2 = []                      #creates empty list to be filled with R2 vectors in VecCircle Loop
BListR3 = []                      #creates empty list to be filled with R3 vectors in VecCircle Loop

#Particle List:
particleList = []

#Force List:
FarrowList = []

#Momentum List:
pList = []

######################
###PRINT STATEMENTS###
######################

#Cylinder Print Statement
print "cylTotal =", cylTotal

#Vector Print Statement
print "vecTotal =", vecTotal


###############
###FUNCTIONS###
###############

def ResetCounts(Total,downList,U,D):
    if U > Total/6 - 1:    #resets upCount to lowest ColorList index at every 1/6th of the length of the solenoid
        U = 0
        return vector(U,D)
    if D < 0:               #resets downCount to highest ColorList index at every 1/6th of the length of the solenoid
        D = len(downList) - 1
        return vector(U,D)
    return vector(U,D)

#Color RGB Functions:       Each is used to alter colors within Cylinder Loop according to the following R,G,B transitions:
                                #1,0,0 > 1,1,0 > 0,1,0 > 0,1,1 > 0,0,1 > 1,0,1 > 1,0,0
def upGreen(upColorList,R,B,upCount):
    return R, upColorList[upCount], B
def downRed(ColorList,G,B,downCount):
    return ColorList[downCount], G, B
def upBlue(upColorList,R,G,upCount):
    return R, G, upColorList[upCount]
def downGreen(ColorList,R,B,downCount):
    return R, ColorList[downCount], B
def upRed(upColorList,G,B,upCount):
    return upColorList[upCount], G, B
def downBlue(ColorList,R,G,downCount):
    return R, G, ColorList[downCount]

#Color Counter Updater:     Calls appropriate Color RGB Function and updates upCount and downCount variables depending on
                                #theta value. There is one if/elif/else statement for each of 6 Color RGB functions.
def ColorCountUpdater(U,D, Step,Limit1,Limit2,Limit3,Limit4,Limit5,Object,Total,upList,downList,reverse):
    UD = vector(U,D)
    #print "UD =",UD
    if reverse == False:
        if Step < Limit1:
            R = 1
            B = 0
            Object.color = upGreen(upList,R,B,U)
            U+= 1
        elif Step < Limit2:
            G = 1
            B = 0
            Object.color = downRed(downList,G,B,D)
            D-= 1
        elif Step < Limit3:
            R = 0
            G = 1
            Object.color = upBlue(upList,R,G,U)
            U+= 1
        elif Step < Limit4:
            R = 0
            B = 1
            Object.color = downGreen(downList,R,B,D)
            D-= 1
        elif Step < Limit5:
            B = 1
            G = 0
            Object.color = upRed(upList,G,B,U)
            U+= 1
        else:
            R = 1
            G = 0
            Object.color = downBlue(downList,R,G,D)
            D-= 1
            
    else:
        if Step < Limit1:
            R = 1
            G = 0
            Object.color = upBlue(upList,R,G,U)
            U+= 1
        elif Step < Limit2:
            G = 0
            B = 1
            Object.color = downRed(downList,G,B,D)
            D-= 1
        elif Step < Limit3:
            R = 0
            B = 1
            Object.color = upGreen(upList,R,B,U)
            U+= 1
        elif Step < Limit4:
            R = 0
            G = 1
            Object.color = downBlue(downList,R,G,D)
            D-= 1
        elif Step < Limit5:
            B = 0
            G = 1
            Object.color = upRed(upList,G,B,U)
            U+= 1
        else:
            R = 1
            B = 0
            Object.color = downGreen(downList,R,B,D)
            D-= 1

    UD = ResetCounts(Total,downList,U,D)

    return UD

#Vector Line Function                   Creates line of vectors along the solenoid axis     
def vecLine(z,vectorList):

    Vector = vector(0,0,z)        
    vectorList.append(Vector)
       
#Vector Circle Function:                Creates set of vector circles along the solenoid axis
def vecCircle(radius,vectorList):

    vecAngle = 0                        #initial angle
    while vecAngle < 2*pi-.001:
        CircleVector = vector(radius*cos(vecAngle), radius*sin(vecAngle), Vz)
        vectorList.append(CircleVector)
        vecAngle+= dvecAngle
        
#Magnetic Field Function:       Creates arrows in same positions as elements in "List," and with axes found by
                                    #adding magnetic field contributions from every cylinder making up the solenoid,
                                    #according to the Biot-Savart Law. Each arrow axis is then added to a
                                    #list called "Blist," which represents Bfield of solenoid and is returned
def Bfield(List,BList,hide):       
    Counter = 0

    while Counter < len(List):      #cycles through List positions and creates arrows at those positions
        B = vector(0,0,0)
        
        cylCounter = 0
    
        while cylCounter < len(cylList):            #adds up Bfield contributions from each cylinder until total B at given
                                                        #position is found, which is used to define arrow axis
            dl = d*cylList[cylCounter].axis/mag(cylList[cylCounter].axis)   #current element
            r = List[Counter] - cylList[cylCounter].pos                     #distance from cylindar to arrow position
            dB = k*mu*I*cross(dl,norm(r))/mag(r)**2                     #magnetic field contribution from single cylinder
            B = B + dB                                                  #updated total Bfield
            cylCounter+= 1                                               

        if hide == False:                                   #no arrows created if hide==False, but BList still updated
            Arrow = arrow(pos = List[Counter], axis = kB*B, opacity = .5)  #creates arrow object

        BList.append(B)                                                 #adds arrow axis to BList
        
        Counter+= 1
        
    return BList


def Bfunc(particlePos):
    B = vector(0,0,0)

    cylCounter = 0
    
    while cylCounter < len(cylList):            #adds up Bfield contributions from each cylinder until total B at given
                                                        #position is found, which is used to define arrow axis
        dl = d*cylList[cylCounter].axis/mag(cylList[cylCounter].axis)   #current element
        r = particlePos - cylList[cylCounter].pos                       #distance from cylinder to arrow position
        dB = k*mu*I*cross(dl,norm(r))/mag(r)**2                         #magnetic field contribution from single cylinder
        B = B + dB                                                      #updated total Bfield
        cylCounter+= 1

    return B


#################
###WHILE LOOPS###
#################

#Cylinder Loop:                 Generates cylinders that make up solenoid, while also updating their color
while theta < RadTotal:              
    xLocatorA = R*cos(theta-dtheta)     #Each Locator finds the x-, y-, or z-position of an imaginary cylinder before (A)
    xLocatorB = R*cos(theta + dtheta)       #or after (B) the real cylinder (named "cyl") created below. The line between  
    yLocatorA = R*sin(theta-dtheta)         #these two imaginary cylinders defines the direction of cyl's axis.
    yLocatorB = R*sin(theta + dtheta)
    zLocatorA = z - dz
    zLocatorB = z + dz

    cyl = cylinder(pos = (h*cos(theta-alpha), h*sin(theta-alpha), z), opacity = .7,
                   axis =(xLocatorB-xLocatorA, yLocatorB-yLocatorA, zLocatorB-zLocatorA), radius=cylR, length=d)
    cylList.append(cyl)

    cylCount = ColorCountUpdater(upCountcyl,downCountcyl,theta,
                                  1*RadTotal/6,2*RadTotal/6,3*RadTotal/6,4*RadTotal/6,5*RadTotal/6,
                                  cyl,cylTotal,upColorList,downColorList,False)

    upCountcyl = cylCount.x
    downCountcyl = cylCount.y
    
    theta+= dtheta
    z+= dz

#Vector Loop:                   Creates vectors representing the direction and magnitude of magnetic field contributions
                                    #from each cylinder making up the solenoid, according to the (scaled) Biot-Savart Law
while Pz <= vecLim:
    vecLine(Pz,vecListPlot)

    Pz+= dPz

while Vz <= vecLim:
    vecLine(Vz,BListLine)                     #creates line of vectors centered along solenoid axis
    vecCircle(R1,BListR1)                     #creates first set of vector circles along solenoid axis
    vecCircle(R2,BListR2)                     #creates second set of vector circles along solenoid axis
    vecCircle(R3,BListR3)                     #creates third set of vector circles along solenoid axis

    Vz+= dVz

vecListCircle = BListR1 + BListR2 + BListR3       #combines all vector lists (line and circles) into one list


############
###Bfield###
############

#BList objects contain magnetic field vector information

BListLine = []                                      #creates empty Bfield list to be used for a list of vectors
BListLine = Bfield(BListLine, BListLine, False)   #arrows displayed (if False) and magnetic field vectors added to list

BListCircle = []
BListCircle = Bfield(vecListCircle, BListCircle, False)

BList = BListLine + BListCircle

BListPlot = []
BListPlot = Bfield(vecListPlot, BListPlot, True)


#Particle (for list)
particle = sphere(pos = (-R/3,R/3,-L/2), radius = cylR/2, opacity = .5, visible = False)
particle.m = m
particle.p = particle.m * vector(.1,.1,.2)

#Farrow (for list)
Farrow = arrow(axis = (0,0,0), color = color.red, visible = False)
Farrow.pos = particle.pos

#Plot display and curves
particlePlot = gdisplay(width = widthPC, height = widthPC*R3/(L/2), xmin = -L/2, xmax = L/2, ymin = -R3, ymax = R3,
                        foreground = color.black, background = color.white,
                        title = "Position = Green, Force = Red, Momentum = Blue, dot(Force,Pos) = magenta")

particleCurve = gcurve(display = particlePlot, dot = True, dot_color = color.yellow, color = color.green)
forceCurve = gcurve(display = particlePlot, color = Farrow.color)
momentumCurve = gcurve(display = particlePlot, color = color.blue)
FPCurve = gcurve(display = particlePlot, color = color.magenta)
##fancyCurve = gcurve(display = particlePlot)
                 
#Bfield Loop
while mag(particle.pos - (0,0,particle.pos.z)) < R3 and particle.pos.z < vecLim:
    particleCurve.plot(pos = (particle.pos.x, particle.pos.y))          #plot particle position
    
    particleList.append(vector(particle.pos))                           #update particle list

    FarrowList.append(vector(Farrow.axis))                              #update Farrow axis list
    
    F = kq*q*cross((particle.p/particle.m),Bfunc(particle.pos))         #force calculation

    Farrow.axis = F	                                    

    forceCurve.plot(pos = (particle.pos.z,kfPlot*Farrow.axis.y))#,kfPlot*Farrow.axis.y))  #plot force
    particle.p = particle.p + F*deltat	                                #momentum update
    particle.pos = particle.pos + particle.p*deltat/particle.m          #position update

    pList.append(vector(particle.p))
    momentumCurve.plot(pos = (particle.pos.z,50*particle.p.y))#,50*particle.p.y))       #plot momentum
    
    Farrow.pos = particle.pos

    FPCurve.plot(pos = (particle.pos.z, dot(norm(Farrow.axis),norm(particle.p))))

while len(particleList)%6 != 0:     #removes last elements of particleList until length of list is divisible by 6
    del particleList[-1]

#Particle (for display)
particleSphere = sphere(pos = particleList[0], radius = 1.5, opacity = .9)

#Particle Trail
trail = curve(radius = .3)

#Farrow (for display)
FarrowVec = arrow(pos = particleList[0], axis = kax*FarrowList[0], opacity = .9)
ParrowVec = arrow(pos = particleList[0], axis = pList[0], opacity = .9, color = color.green)

#Particle Color Constants
dColorSphere = 6/len(particleList)
upColorListSphere = arange(0, 1, dColorSphere)
downColorListSphere = upColorListSphere + dColorSphere
upCountSphere = upCountVec =  upCountMomentum = 0
downCountSphere = downCountVec = downCountMomentum = len(downColorListSphere) - 1

#Particle and Farrow display loop
for z in arange(0,len(particleList)-1):
    rate(100)
    
    sphereCount = ColorCountUpdater(upCountSphere, downCountSphere, particleSphere.pos.z,
                                 particleList[int(1*len(particleList)/6-1)].z, particleList[int(2*len(particleList)/6-1)].z,
                                 particleList[int(3*len(particleList)/6-1)].z, particleList[int(4*len(particleList)/6-1)].z,
                                 particleList[int(5*len(particleList)/6-1)].z,
                                 particleSphere, len(particleList), upColorListSphere, downColorListSphere, True)
    upCountSphere = sphereCount.x
    downCountSphere = sphereCount.y

    vecCount = ColorCountUpdater(upCountVec, downCountVec, particleSphere.pos.z,
                                 particleList[int(1*len(particleList)/6-1)].z, particleList[int(2*len(particleList)/6-1)].z,
                                 particleList[int(3*len(particleList)/6-1)].z, particleList[int(4*len(particleList)/6-1)].z,
                                 particleList[int(5*len(particleList)/6-1)].z,
                                 FarrowVec, len(particleList), upColorListSphere, downColorListSphere, False)
    upCountVec = vecCount.x
    downCountVec = vecCount.y

    momentumCount = ColorCountUpdater(upCountMomentum, downCountMomentum, particleSphere.pos.z,
                                 particleList[int(1*len(particleList)/6-1)].z, particleList[int(2*len(particleList)/6-1)].z,
                                 particleList[int(3*len(particleList)/6-1)].z, particleList[int(4*len(particleList)/6-1)].z,
                                 particleList[int(5*len(particleList)/6-1)].z,
                                 ParrowVec, len(particleList), upColorListSphere, downColorListSphere, True)
    upCountMomentum = vecCount.x
    downCountMomentum = vecCount.y

##    fancyCurve.plot(pos = (particleSphere.pos.x, particleSphere.pos.y), color = particleSphere.color)
    
    FarrowVec.pos = particleList[z]
    FarrowVec.axis = kax*FarrowList[z]
    
    ParrowVec.pos = particleList[z]
    ParrowVec.axis = 50*pList[z]
    
    particleSphere.pos = particleList[z]
    trail.append(pos = particleSphere.pos, color = particleSphere.color)

###########
###PLOTS###
###########

Bplot = gdisplay(x = 0, y = particlePlot.height, width = widthPC, height = heightPC - particlePlot.height, 
                 title='Bmag vs. Z', foreground = color.black, background = color.white)

Bbars = gvbars(display = Bplot, delta = (9/10)*dPz, color = (211/255,94/255,96/255))
Bcurve = gcurve(display = Bplot, color = (57/255,106/255,177/255))

for Z in arange(0,len(vecListPlot)):
    Bbars.plot(pos = (vecListPlot[Z].z,mag(BListPlot[Z])))
    Bcurve.plot(pos = (vecListPlot[Z].z,mag(BListPlot[Z])))


