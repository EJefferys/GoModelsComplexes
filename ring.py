### INTERACTION PARTICLE STUFF ###
### Parts related to making the interaction particles are indicated with ### IP ### to make them easier to find (and sorry that this code isn't super tidy or pretty atm) ###


import os
import random
import MDAnalysis as mda
import numpy as np
import MDAnalysis.analysis.contacts as contacts
import IMP
import IMP.pmi
import IMP.atom
import IMP.algebra
import IMP.rmf
import IMP.core
import RMF
import IMP.container
import IMP.display
import itertools
import secrets

def convertToFrames(time_ns,stepsize_fs):
    FS_PER_NS=1E6
    time_fs=time_ns*FS_PER_NS
    n_frames_float=(time_fs+0.0)/stepsize_fs
    n_frames=int(round(n_frames_float))
    return max(n_frames,1)
    
def getRandomVector(L):
    """get random vector within box of size L"""

    x1,x2=-L/2,L/2
    y1,y2=-L/2,L/2
    z1,z2=-L/2,L/2
    number=secrets.randbits(128)
    rng1=np.random.default_rng(number)
    v=rng1.uniform([x1,y1,z1],[x2,y2,z2],size=(1,3))
    return list(v[0])

def haversine(v):
    """get arc distance between particles"""

    radius=50
    angle=(np.pi*2)/3
    arc=angle*50
    return arc


### IP vvv ###

def createSimpleParticle(m,name,radius,mass,colour,v):
    """create a simple IMP particle"""

    p=IMP.Particle(m,name)
    xyzr=IMP.core.XYZR.setup_particle(p)
    xyzr.set_coordinates_are_optimized(True)
    xyzr.set_coordinates(v)
    xyzr.set_radius(radius)
    IMP.atom.Mass.setup_particle(p,mass)
    IMP.atom.Hierarchy.setup_particle(p)
    IMP.display.Colored.setup_particle(p,IMP.display.get_display_color(colour))
    IMP.atom.Diffusion.setup_particle(p)
    return p

def createCombinedParticle(m,p1,v1,names):
    """create interaction patches"""

    v2=[v1[0]+50,v1[1],v1[2]]
    v3=[v1[0]-50,v1[1],v1[2]]
    p2=createSimpleParticle(m,names[0],1,1,0,v2)
    p3=createSimpleParticle(m,names[1],1,1,0,v3)
    return [p2,p3]


### IP ^^^ ###


# set up model
m=IMP.Model()
pRoot=IMP.Particle(m,'root')
hRoot=IMP.atom.Hierarchy.setup_particle(pRoot)
rs=[] # restrictions

L=2000

particleList=[]
interactionParticleList=[]

for i in range(10):
    print(i)
    name='p'+str(i)

    ### IP vvv ###

    v=getRandomVector(L)
    p=createSimpleParticle(m,name,50,10,i+1,v) # sets up the host particle
    particleList.append(p)
    hP1=IMP.atom.Hierarchy(p)
    hRoot.add_child(hP1)

    ip=createCombinedParticle(m,p,v,[name+'-A',name+'-B']) # creates two interaction particles for the just created host particle
    interactionParticleList.append(ip[0])
    interactionParticleList.append(ip[1])
    hP2=IMP.atom.Hierarchy(ip[0])
    hP3=IMP.atom.Hierarchy(ip[1])
    hRoot.add_child(hP2)
    hRoot.add_child(hP3)
    distance=haversine(v)
    f=IMP.core.Harmonic(distance,1)
    s=IMP.core.DistancePairScore(f) # fixes the distance between the interaction particles on the surface of the other particle ie the space between them is kept constant
    r=IMP.core.PairRestraint(m,s,(ip[0],ip[1]))
    rs.append(r)
    bb=IMP.core.Harmonic(0,1)
    ss=IMP.core.SphereDistancePairScore(bb) # fixes the interaction particles to be on the surface
    r1=IMP.core.PairRestraint(m,ss,(p,ip[0]))
    r2=IMP.core.PairRestraint(m,ss,(p,ip[1]))
    rs.append(r1)
    rs.append(r2)

    ### IP ^^^ ###

pairList=[]
for i in range(len(interactionParticleList[1::2])-1):
    pairList.append([
        interactionParticleList[::2][i],
        interactionParticleList[1::2][i+1]
        ])

pairList.append([interactionParticleList[1::2][0],interactionParticleList[::2][-1]])

for pair in pairList:
    print(pair)
    p1=pair[0]
    p2=pair[1]
    bb=IMP.core.TruncatedHarmonicBound(0,0.1,30)
    ss=IMP.core.SphereDistancePairScore(bb)
    r5=IMP.core.PairRestraint(m,ss,(p1,p2))
    rs.append(r5)


RMF_FILENAME='toy-model.rmf'
K_BB=0.1
K_EXCLUDED=0.1
BD_STEP_SIZE_SEC=10E-10#10E-10
SIM_TIME_SEC=1.5
bd_step_size=BD_STEP_SIZE_SEC*1E+15
sim_time_ns=SIM_TIME_SEC*1E+9
RMF_DUMP_INTERVAL_NS=sim_time_ns/10000.0
bb=IMP.algebra.BoundingBox3D(IMP.algebra.Vector3D(-L/2,-L/2,-L/2),IMP.algebra.Vector3D(L/2,L/2,L/2))
bb_harmonic=IMP.core.HarmonicUpperBound(0,K_BB)
outer_bbss=IMP.core.BoundingBox3DSingletonScore(bb_harmonic,bb)
rs.append(IMP.container.SingletonsRestraint(outer_bbss,hRoot.get_children()))

# Excluded volume restraint
ev=IMP.core.ExcludedVolumeRestraint(IMP.atom.get_leaves(hRoot),K_EXCLUDED)
ev=IMP.core.ExcludedVolumeRestraint(IMP.atom.get_leaves(hRoot),K_EXCLUDED,10,'EV')
rs.append(ev)
sf=IMP.core.RestraintsScoringFunction(rs,'SF')
bd=IMP.atom.BrownianDynamics(m)
bd.set_log_level(IMP.SILENT)
bd.set_scoring_function(sf)
bd.set_maximum_time_step(bd_step_size/1000)
bd.set_temperature(300)

# run BD simulation
sim_time_frames=convertToFrames(sim_time_ns,bd_step_size)
rmf_dump_interval_frames=convertToFrames(RMF_DUMP_INTERVAL_NS,bd_step_size)
rmf=RMF.create_rmf_file(RMF_FILENAME)
rmf.set_description('some stuff that is not relevant rn\n')
IMP.rmf.add_restraints(rmf,rs)
IMP.rmf.add_hierarchy(rmf,hRoot)
IMP.rmf.add_geometry(rmf,IMP.display.BoundingBoxGeometry(bb))
sos=IMP.rmf.SaveOptimizerState(m,rmf)
sos.set_log_level(IMP.SILENT)
sos.set_simulator(bd)
sos.set_period(rmf_dump_interval_frames)
bd.add_optimizer_state(sos)
sos.update_always('initial conformation')
print('maybe running the simulation?')
print('score before: %s'%(sf.evaluate(True)))
bd.optimize(sim_time_frames)
print('score after: %s'%(sf.evaluate(True)))
