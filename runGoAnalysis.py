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

import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import MDAnalysis
import IMP
import RMF
import os
import glob


def intersection(lst1,lst2):
    lst3=[value for value in lst1 if value in lst2]
    return list(set(lst3))

def show_xml(nh,kcs):
    name = nh.get_name()
    name.replace(" ", "_")
    f.write("<node name=\"" + name + "\" id=\"" + str(nh.get_id().get_index())\
    + "\" type=\"" + str(nh.get_type()) + "\"/>\n")
    if verbose:
        for kc in kcs:
            show_data_xml(nh, kc)
    children = nh.get_children()
    for c in children:
        f.write("<child>\n")
        show_xml(c, kcs)
        f.write("</child>\n")

def show_data_xml(nh,kc):
    rh = nh.get_file()
    keys = rh.get_keys(kc)
    opened = False
    for k in keys:
        v = nh.get_value(k)
        if v is not None:
            if not opened:
                f.write("< "+rh.get_name(kc))
                f.write('\n')
                opened = True
            name = rh.get_name(k)
            name.replace(" ", "_")
            f.write(name+" =\"" + str(v) + "\"\n")
    if opened:
        f.write("/>\n")

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
    sliceStep=int(args.particles)+1 # need to slice with one more than the number of interaction particles
else:
    sliceStep=None
L=int(args.boxsize)
nParticles=int(args.nparticles)
radius=float(args.radius)
cutoff=(radius*2)+10

nativePairList=[]

if args.infile is None: # if no input file of interactions provided then default to a ring
    sliceStep=3
    ring=True
    for i in range(nParticles-1):
        p0='p%s'%(i)
        p1='p%s'%(i+1)
        #nativePairList.append([p0,p1])
        nativePairList.append((p0,p1))
    p0='p0'
    #nativePairList.append([p0,p1])
    nativePairList.append((p0,p1))
else:                   # otherwise get the interactions from the input list
    ring=False
    f=open(args.infile,'r')
    pairIndices=f.readlines()
    f.close()
    for pair in pairIndices:
        p0=pair.split(',')[0]
        p1=pair.split(',')[1]
        p01='p%s'%(p0)
        p11='p%s'%(p1)
        #nativePairList.append([p01,p11])
        nativePairList.append((p01,p11))


################## NEED: file name, number of particles, native interaction pairs, number of interaction particles, radius
# don't really need to rebuild the model, just make the interaction list

### set up for the analysis ###

file_name='%s.rmf'%(fname)
myPath=os.getcwd()
os.system('rm \#* *.xtc *.pdb *.xml dat.xyz list-xtc.txt times-pdb.txt')
verbose=True

rh=RMF.open_rmf_file_read_only(file_name)
rh.set_current_frame(RMF.FrameID(0))
for i in range(len(rh.get_frames())):
    rh.set_current_frame(RMF.FrameID(i))
    f=open('%s.xml'%(i),'w')
    f.write("<?xml version=\"1.0\"?>\n")
    f.write("<rmf>\n")
    f.write("<path>\n")
    f.write(file_name+'\n')
    f.write("</path>\n")
    f.write("<description>\n")
    f.write(rh.get_description()+'\n')
    f.write("</description>\n")
    f.write("<path>\n")
    f.write('%s\n'%(myPath))
    f.write("</path>\n")
    kcs=rh.get_categories()
    show_xml(rh.get_root_node(),kcs)
    f.write("</rmf>\n")
    f.close()

xmlCounter=len(glob.glob1(myPath,'*.xml'))

frameList=[]
coordinateList=[]
particleList=[]
radiusList=[]

string="< physics"
f=open('dat.xyz','w')
for i in range(xmlCounter):
    frameList.append(i)
    f.write('%s\n'%(nParticles))
    f.write('description\n')
    fl=open('%s.xml'%(i),'r')
    lines=fl.readlines()
    fl.close()
    tmp=[]
    for j in range(len(lines)):
        if string in lines[j]:
            name=lines[j-1].split('"')[1]
            particleList.append(name)
        if 'coordinates =' in lines[j]:
            coords=lines[j].split('[')[1].split(']')[0]
            xyz=coords.split(',')
            xyz=[float(x) for x in xyz]
            f.write('%s %.8f %.8f %.8f\n'%(name,xyz[0],xyz[1],xyz[2]))
            tmp.append([float(xyz[0]),float(xyz[1]),float(xyz[2])])
        if 'radius =' in lines[j]:
            radiusList.append(lines[j].split('"')[1])
    coordinateList.append(tmp)
    f.write('\n')

f.close()

for i in range(xmlCounter):
    f=open('%s.pdb'%(i),'w')
    f.write('CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1           1\n')
    fl=open('%s.xml'%(i),'r')
    lines=fl.readlines()
    fl.close()
    w=0
    for j in range(len(lines)):
        if string in lines[j]:
            if 'name' in lines[j-1]:
                print(lines[j-1])
                tmp=lines[j-1].split('=')[1]
                name=tmp.split('"')[1]
                if '-' in name:
                    name=name.split('-')[1]
                else:
                    name=name
   #         name=name[7:]
                print(name)
        if 'coordinates =' in lines[j] and '"p' in lines[j-5]:
            print(j)
            coords=lines[j].split('[')[1].split(']')[0]
            xyz=coords.split(',')
            xyz=[float(x) for x in xyz]
            f.write('{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}\n'.format('ATOM',w,'O',' ',name,'A',w,' ',xyz[0],xyz[1],xyz[2],1,1,'O',str(0)))
            w+=1
    f.close()

### make the trajectory ###

os.system('bash combine.sh')
os.system('cp 0.pdb traj-%s.pdb'%(fname))
os.system('gmx trjcat -f `cat list-xtc.txt` -cat -o traj-%s.xtc'%(fname))

### run analyses ###

u=MDAnalysis.Universe('traj-%s.pdb'%(fname),'traj-%s.xtc'%(fname))
residues=u.select_atoms('resname p*')

# set up to collect data from the trajectory

i=0
qList=[]
scoreList=[]
nContactList=[]
contactsPerFrame=[]
frameList=[]
seenGroups=[]
groupList=[]
scoreList=[]
datList=[]
frameList=[]
contactGroups=[]
dat=[]
seenComplex=[]
passageTime=[]
completeFrames=[]


# loop through the trajectory and collect the data

for ts in u.trajectory:

    # frames
    print(i)
    frameList.append(i)

    # contacts
    nContacts=0     # count the total number of contacts in each frame
    contactPairs=[]     # get the pairs of particles that are in contact in each frame
    seen=[]
    for res1 in u.residues[::sliceStep]:        # go through pairs of particles and avoid duplicates
        for res2 in u.residues[::sliceStep]:
            if res1==res2:
                next
            elif (res1,res2) in seen:
                next
            elif (res2,res1) in seen:
                next
            else:
                seen.append((res1,res2))
                seen.append((res2,res1))
                a1=res1.atoms.positions
                a2=res2.atoms.positions
                iscontact=checkContact(a1,a2,cutoff)       # checks contact based on a cutoff that is double the particle radii
                if iscontact==True:
                    contactPairs.append([res1.resname,res2.resname])    # add the pair to the list of pairs in contact for that frame
                    nContacts+=1        # add the contact to the total number of particles in contact
            contactsPerFrame.append(contactPairs)       # this is all the contacts in that frame in a list ie get a list of lists that will be used in getting the pathway
    nContactList.append(nContacts)      # number of contact
    Q=len(intersection(contactPairs,nativePairList))        # fraction of native contacts goes into qList
    qList.append(Q/len(nativePairList))


    # score
    frame=i
    score=getScore(i)
    scoreList.append(score)


    # pathway and frame at which each subcomplex first forms
    datList=[]

    G=toGraph(contactPairs)
    datList=[tuple(c) for c in sorted(connected_components(G),key=len,reverse=True)]        # this turns the contact pairs into groups

    for rname in u.residues[::sliceStep].resnames:
        if rname not in datList:
            datList.append((rname,))    # adds single particles into the previous because they would be missed in the contact pairs


    # remove duplicates in a couple of steps
    newList=[]
    for resname in u.residues[::sliceStep].resnames:
        tmpList=[j for j in datList if resname in j]
        longest=max(tmpList,key=len)
        newList.append(longest)
    datList=newList
    for d in datList:
        if len(d)==nParticles:
            completeFrames.append(i)
    newList=[]
    for ele in datList:
        newList.append(sorted(ele))
    datList=removeDuplicates(datList)

    if datList not in seenComplex:
        seenComplex.append(str(datList)) # this should give a list of all the complexes in the order in which they appear
        passageTime.append(frame) # this should give the frames at which each of the complexes appear
    groupList.append(datList)

    i+=1        # for the frame counter


# SAVE STUFF #

# pathway
pathDict={'frame':passageTime,'complex':seenComplex}        # each subcomplex seen and the time at which it forms
df=pd.DataFrame(pathDict)
df.to_pickle('passagetime-%s.pkl'%(fname))

endFrame=sorted(completeFrames)[0]

if ring==True: # repeat but account for the fact that all particles are equivalent
    newTimes=[]
    newComplexes=[]
    for i in range(len(seenComplex)):
        time=passageTime[i]
        c1=seenComplex[i]
        c2=[]
        for x in c1: # should be subcomplexes within stuff
            tmp=[1 for y in x]
            c2.append(tmp)
#            c2.append(len(x))
        c2=sorted(c2)
        newComplexes.append(c2)
        newTimes.append(time)
    pathDict={'frame':newTimes,'complex':newComplexes}
    df=pd.DataFrame(pathDict)
    df.to_pickle('passagetime_ones-%s.pkl'%(fname))

# per frame data
datDict={'scores':scoreList,'frames':frameList,'groups':groupList,'contacts':nContactList,'Q':qList}      # score, Q, group, number of contacts per frame
df=pd.DataFrame(datDict)
df.to_pickle('perframedata-%s.pkl'%(fname))

#print(len(seenGroups)/len(groupOptions)) # what fraction of options are explored


# PLOTS #

# setup style
sns.set_style()
sns.set(font_scale=0.5)

# per frame data
fig,(ax1,ax2,ax3)=plt.subplots(3,sharex=True,figsize=(20,10))
fig.suptitle('Complex assembly %s'%(fname),fontsize=20)
ax1.plot(list(range(len(qList))),qList)
ax1.set_ylabel('Q',fontsize=18)
ax1.tick_params(axis='y',labelsize=16)
ax2.scatter(list(range(len(scoreList[1:]))),scoreList[1:])
ax2.set_ylabel('E',fontsize=18)
ax2.tick_params(axis='y',labelsize=16)
ax3.scatter(list(range(len(nContactList))),nContactList)
ax3.set_ylabel('N',fontsize=18)
ax3.tick_params(axis='y',labelsize=16)
plt.xlabel('frame',fontsize=18)
plt.xticks(fontsize=16)
plt.savefig('plots-%s.png'%(fname),dpi=300)
plt.close()


# lineplot
particles=list(u.residues.resnames)[::sliceStep]

startPositions=[]
for i in range(len(u.residues[::sliceStep])):
    startPositions.append(i)
newPositions=startPositions

i=0
dat=[]
frames=[]
w=0
currentFrame=0

for i in range(len(u.residues[::sliceStep])):
    dat.append([])


for ts in u.trajectory:
    seen=[]

    # get the residue pair
    for i in range(len(u.residues[::sliceStep])):
        res1=u.residues[::sliceStep][i]
        for j in range(len(u.residues[::sliceStep])-1):
            res2=u.residues[::sliceStep][j]

            # get the positions in the plot of the residues
            pos1=newPositions[i]
            pos2=newPositions[j]

            # exclude identical residues or pairs that have already been seen in this frame
            if i==j or (i,j) in seen or (j,i) in seen:
                next

            else:

                # if new pair add to seen so don't repeat
                seen.append((i,j))
                seen.append((j,i))

                # get positions of the residues and the distance between them
                a1=res1.atoms.positions
                a2=res2.atoms.positions
                distance=np.linalg.norm(a1-a2)
                print(w)

                frames.append(w)
                w+=1
                # are they in contact?
                if distance<=cutoff:

                    # if in contact get a new position for them
                    combined=np.mean([pos1,pos2])
                    newPositions[i]=combined-0.05
                    newPositions[j]=combined+0.05

                    # add position to that part of the data lists
                    dat[i].append(combined-0.05)
                    dat[j].append(combined+0.05)

                # if not in contact add the previous position again
                else:
                    dat[i].append(pos1)
                    dat[j].append(pos2)

    diffList=[]
    for pos1 in newPositions:
        for pos2 in newPositions:
            if pos1==pos2:
                next
            else:
                diff=np.linalg.norm(pos1-pos2)
                diffList.append(diff)
    if max(diffList)-min(diffList)<=nParticles*0.05:
        endFrame=currentFrame
    currentFrame+=1

for subdat in dat:
    i=dat.index(subdat)
    plt.plot(subdat,list(range(len(subdat))),label=particles[i])
plt.xticks(np.linspace(0,nParticles,len(particles)),particles)
plt.legend(loc='lower left',bbox_to_anchor=(1,0.5))
plt.gca().invert_yaxis()
plt.savefig('%s-lineplot.png'%(fname),dpi=300)
plt.close()

#print(endFrame)
#print(len(subdat)/5)
#print(endFrame)
truncate=endFrame*5#:w*5#len(subdat)-int(endFrame)

for subdat in dat:
    i=dat.index(subdat)
    plt.plot(subdat[:truncate],list(range(len(subdat)))[:truncate],label=particles[i])
plt.xticks(np.linspace(0,nParticles,len(particles)),particles)
plt.legend(loc='lower left',bbox_to_anchor=(1,0.5))
plt.gca().invert_yaxis()
plt.savefig('%s-lineplot-truncated.png'%(fname),dpi=300)
plt.close()
