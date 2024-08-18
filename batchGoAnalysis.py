### NOTE: Current version just runs the analysis on each of the repeats, and nothing else.

import os
import argparse

parser=argparse.ArgumentParser()
parser.add_argument('-p','--particles',default=None,help='this should be used only if there are interaction particles, and should be one more than the number of interaction particles per main particle e.g interaction sites would be -p 3')
parser.add_argument('-b','--boxsize',default=2000,help='size of the simulation box, default is 2000 A')
parser.add_argument('-n','--nparticles',default=10,help='number of particles in the model, default is 10')
parser.add_argument('-r','--radius',default=50,help='radius in Angstroms of the particles, default is 50')
parser.add_argument('-i','--infile',default=None,help='file with interactions between particles')
parser.add_argument('-s','--repeats',default=10,help='number of repeats to do')
args=parser.parse_args()

interactionParticles=args.particles
nParticles=args.nparticles
inputfile=args.infile
box=args.boxsize
radius=args.radius
nRepeats=int(args.repeats)

cwd=os.getcwd()
for i in range(nRepeats):
    os.system('cp runGoAnalysis.py %s'%(i))
    os.system('cp combine.sh %s'%(i))
    if inputfile!=None:
        os.system('cp %s %s'%(inputfile,i))
    os.chdir('%s/%s'%(cwd,i))
    os.system('rm *.xml *.xtc *.txt *.pdb')
    os.system('python3 runGoAnalysis.py -f %s -n %s -p %s -i %s -b %s -r %s'%(fname,nParticles,interactionParticles,inputfile,box,radius))
    os.chdir(cwd)
