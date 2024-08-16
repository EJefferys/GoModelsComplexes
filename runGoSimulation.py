# IMPORT #

# various functions for the simulations and their analysis
from GoModelFunctions import *

# general
import os
import random
import MDAnalysis
import numpy as np
import MDAnalysis.analysis.contacts as contacts
import argparse
import glob
import itertools
import secrets

# IMP
import IMP
import IMP.pmi
import IMP.atom
import IMP.algebra
import IMP.rmf
import IMP.core
import RMF
import IMP.container
import IMP.display

# plotting
import matplotlib.pyplot as plt
import seaborn as sns


# ARGUMENT PARSER #

######################### want to generalise to any system so should be able to specify a file for the interactions between particles and the number of particles as input or default to a ring structure ###################
########################## also add options for radius and mass of the particles #######################

parser=argparse.ArgumentParser()
parser.add_argument('-f','--filename',default='0',help='filename to be used for the rmf output of IMP and all subsequently generated processing and analysis files')
parser.add_argument('-p','--particles',default=None,help='this should be used only if there are interaction particles, and should be one more than the number of interaction particles per main particle e.g. two interaction sites would be -p 3')
parser.add_argument('-b','--boxsize',default=2000,help='size of the simulation box, default is 2000 A')
parser.add_argument('-n','--nparticles',default=10,help='number of particles in the model, default is 10')
parser.add_argument('-r','--radius',default=50,help='radius in Angstroms of the particles, default is 50')
parser.add_argument('-i','--infile',default=None,help='file with interactions between particles')
args=parser.parse_args()

fname=args.filename
if args.particles!=None:
    sliceStep=int(args.particles)+1
else:
    sliceStep=None
L=int(args.boxsize)
nParticles=int(args.nparticles)
radius=float(args.radius)


if args.infile is None: # if no input file of interactions provided then default to a ring
    sliceStep=3
    ring=True
else:                   # otherwise get the interactions from the input list
    ring=False
    f=open(args.infile,'r')
    pairIndices=f.readlines()
    f.close()

# SET UP THE MODEL #

m=IMP.Model()
pRoot=IMP.Particle(m,'root')
hRoot=IMP.atom.Hierarchy.setup_particle(pRoot)
rs=[]


# CREATE THE PARTICLES #

particleList=[]     # list of the main particles
interactionParticleList=[]      # list of the interaction particles
interactionParticleMap={}       # mapping of the interaction particles to their host main particle
mainParticleMap={}

# for the number of particles given, loop through and create the particles

for i in range(nParticles):

    # main particle made and put at a random position in the box
    name='p'+str(i)
    v=getRandomVector(L)
    p=createSimpleParticle(m,name,radius,10,i+1,v)
    particleList.append(p)
    hP1=IMP.atom.Hierarchy(p)
    hRoot.add_child(hP1)
    

    # if the system has interaction particles then make them next
    if sliceStep!=None:
        ip=createCombinedParticle(m,p,v,[name+'-A',name+'-B'])
        interactionParticleList.append(ip[0])
        interactionParticleList.append(ip[1])
        hP2=IMP.atom.Hierarchy(ip[0])
        hP3=IMP.atom.Hierarchy(ip[1])
        hRoot.add_child(hP2)
        hRoot.add_child(hP3)
        distance=haversine(radius,3)
        f=IMP.core.Harmonic(distance,1)
        s=IMP.core.DistancePairScore(f)
        r=IMP.core.PairRestraint(m,s,(ip[0],ip[1]))
        rs.append(r)
        bb=IMP.core.Harmonic(0,1)
        ss=IMP.core.SphereDistancePairScore(bb)
        r1=IMP.core.PairRestraint(m,ss,(p,ip[0]))
        r2=IMP.core.PairRestraint(m,ss,(p,ip[1]))
        rs.append(r1)
        rs.append(r2)
        interactionParticleMap[ip[0]]=p
        interactionParticleMap[ip[1]]=p
        mainParticleMap[p]=(ip[0],ip[1])


# set up the interactions

pairList=[]

if ring==True:      # if the structure is a ring
    for i in range(len(interactionParticleList[1::2])-1):
        pairList.append([
            interactionParticleList[::2][i],
            interactionParticleList[1::2][i+1]
        ])

    pairList.append([interactionParticleList[1::2][0],interactionParticleList[::2][-1]])

elif ring==False and sliceStep!=None:       # not a ring and interaction aprticles needed
    for pair in pairIndices:
        i0=int(pair.split(',')[0])
        i1=int(pair.split(',')[1])
        p0=particleList[i0]
        p1=particleList[i1]
        pairList.append([mainParticleMap[p0][1],mainParticleMap[p1][0]])

else:       # not a ring and no interaction particles
    for pair in pairIndices:
        p0=int(pair.split(',')[0].strip())
        p1=int(pair.split(',')[1].strip())
        pairList.append([particleList[p0],particleList[p1]])

for pair in pairList:
    p1=pair[0]
    p2=pair[1]
    bb=IMP.core.TruncatedHarmonicBound(0,0.1,30) # truncated harmonic so don't just agglomerate together
    ss=IMP.core.SphereDistancePairScore(bb)
    r5=IMP.core.PairRestraint(m,ss,(p1,p2))
    rs.append(r5)


# SETUP BROWNIAN DYNAMICS #

RMF_FILENAME='%s.rmf'%(fname)
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
ev=IMP.core.ExcludedVolumeRestraint(IMP.atom.get_leaves(hRoot),K_EXCLUDED)
ev=IMP.core.ExcludedVolumeRestraint(IMP.atom.get_leaves(hRoot),K_EXCLUDED,10,'EV')
rs.append(ev)
sf=IMP.core.RestraintsScoringFunction(rs,'SF')
bd=IMP.atom.BrownianDynamics(m)
bd.set_log_level(IMP.SILENT)
bd.set_scoring_function(sf)
bd.set_maximum_time_step(bd_step_size/1000)
bd.set_temperature(300)


# RUN BROWNIAN DYNAMICS #

sim_time_frames=convertToFrames(sim_time_ns,bd_step_size)
rmf_dump_interval_frames=convertToFrames(RMF_DUMP_INTERVAL_NS,bd_step_size)
rmf=RMF.create_rmf_file(RMF_FILENAME)
rmf.set_description('Brownian Dynamics simulation of Go type macromolecular complex formation\n')
IMP.rmf.add_restraints(rmf,rs)
IMP.rmf.add_hierarchy(rmf,hRoot)
IMP.rmf.add_geometry(rmf,IMP.display.BoundingBoxGeometry(bb))
sos=IMP.rmf.SaveOptimizerState(m,rmf)
sos.set_log_level(IMP.SILENT)
sos.set_simulator(bd)
sos.set_period(rmf_dump_interval_frames)
bd.add_optimizer_state(sos)
sos.update_always('Initial conformation')
print('Running the simulation')
print('Score before: %s'%(sf.evaluate(True)))
bd.optimize(sim_time_frames)
print('Score after: %s'%(sf.evaluate(True)))
