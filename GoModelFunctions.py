### set of functions for use in Go model protein complex assembly simulations and their analysis ###

# IMPORTS #

import os
import glob
import random
import MDAnalysis as mda
import matplotlib.pyplot as plt
import seaborn as sns
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

# SETUP #

def createSimpleParticle(m,name,radius,mass,colour,v):
    """create a simple IMP particle with specified mass, radius, colour and coordinates (v)"""

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
    """create interaction patches for particle p1, whose coordinates are v1"""

    v2=[v1[0]+50,v1[1],v1[2]]
    v3=[v1[0]-50,v1[1],v1[2]]
    p2=createSimpleParticle(m,names[0],1,1,0,v2)
    p3=createSimpleParticle(m,names[1],1,1,0,v3)
    return [p2,p3]

def getRandomVector(L):
    """get random vector within box of size L -- by using secrets module this should be closer to truly random"""

    x1,x2=-L/2,L/2
    y1,y2=-L/2,L/2
    z1,z2=-L/2,L/2
    number=secrets.randbits(128)
    rng1=np.random.default_rng(number)
    v=rng1.uniform([x1,y1,z1],[x2,y2,z2],size=(1,3))
    return list(v[0])

def convertToFrames(time_ns,stepsize_fs):
    """converts the BD simulation time and stepsize into the number of frames"""

    FS_PER_NS=1E6
    time_fs=time_ns*FS_PER_NS
    n_frames_float=(time_fs+0.0)/stepsize_fs
    n_frames=int(round(n_frames_float))
    return max(n_frames,1)

def haversine(r,f): #################### needs tidying up and making more generalisable
    """get arc distance between interaction particles on the surface of a particle of radius r"""

    angle=(np.pi*2)/f # choosing the value of f influences the spacing of the particles and needs to be optimised for each system -- I did by trial and error but there is presumably a better way to do it
    arc=angle*r
    return arc


# FORMAT CONVERSION #

def intersection(lst1,lst2):
    lst3=[value for value in lst1 if value in lst2]
    return list(set(lst3))

def make_xml(i,file_name,rh,myPath,verbose=True):#show_xml(nh,kcs):

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
    nh=rh.get_root_node()
    #show_xml(rh.get_root_node(),kcs)
    name = nh.get_name()
    name.replace(" ", "_")
    f.write("<node name=\"" + name + "\" id=\"" + str(nh.get_id().get_index())\
    + "\" type=\"" + str(nh.get_type()) + "\"/>\n")
    if verbose:
        for kc in kcs:
            show_data_xml(nh, kc,f)
    children = nh.get_children()
    for c in children:
        f.write("<child>\n")
        name=c.get_name()
        name.replace(" ", "_")
        f.write("<node name=\"" + name + "\" id=\"" + str(c.get_id().get_index())\
            + "\" type=\"" + str(c.get_type()) + "\"/>\n")
        for kc in kcs:
            show_data_xml(c,kc,f)
#        show_xml(c, kcs)
        f.write("</child>\n")
    f.write("</rmf>\n")
    f.close()

def show_data_xml(nh,kc,f):
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


def xmlToPDB(fname):
    """convert xml file fname.xml to PDB file because PDB is more handy for people with MD backgrounds"""

    f=open('%s.pdb'%(fname),'w')
    f.write('CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1           1\n')
    fl=open('%s.xml'%(fname),'r')
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
        if 'coordinates =' in lines[j] and '"p' in lines[j-5]:
            coords=lines[j].split('[')[1].split(']')[0]
            xyz=coords.split(',')
            xyz=[float(x) for x in xyz]
            f.write('{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}\n'.format('ATOM',w,'O',' ',name,'A',w,' ',xyz[0],xyz[1],xyz[2],1,1,'O',str(0)))
            w+=1
    f.close()



# ANALYSIS #

#def Manhattan(tup1,tup2):
#   return np.linalg.norm(np.asarray(tup1)-np.asarray(tup2))

def toGraph(l):
    """convert list l into a networkx graph"""

    G=networkx.Graph()
    for part in l:
        # each sublist is a bunch of nodes
        G.add_nodes_from(part)
        # it also imlies a number of edges:
        G.add_edges_from(toEdges(part))
    return G

def toEdges(l):
    """
        treat `l` as a Graph and returns it's edges
        toEdges(['a','b','c','d']) -> [(a,b), (b,c),(c,d)]
    """
    it=iter(l)
    last=next(it)

    for current in it:
        yield last,current
        last=current

def getTrajectoryContacts(fname,cutoff,skip=0):
    """ get a list of the contacts present in each frame of a simulation"""

    u=MDAnalysis.Universe('%s.pdb'%(fname),'%s.xtc'%(fname))

    contactsPerFrame=[]
    seen=[]

    if skip==0:
        for ts in u.trajectory:
            contactPairs=[]

            for res1 in u.residues:
                for res2 in u.residues:
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
                        iscontact=checkContact(a1,a2,cutoff)
                        if iscontact==True:
                            contactPairs.append([res1.resname,res2.resname])
            contactsPerFrame.append(contactPairs)
#            frameList.append(i)
#            nContacts=0
#            pairList=[]
#            for j in range(len(residues)):
#                for k in range(len(residues)):
#                    if j==k:
#                        next
#                    else:
#                        coords0=residues[j].position
#                        coords1=residues[k].position
#                    distance=np.linalg.norm(coords0-coords1)
#                    if distance<=110:
#                        pairList.append((residues.resnames[j],residues.resnames[k]))
#                        pairList.append((residues.resnames[k],residues.resnames[j]))
#                        pairList=list(set(pairList))
#                        nContacts+=1
#        count=len(intersection(pairList,nativePairList))
#        qList.append(count/len(nativePairList))
#        nContactList.append(nContacts)
#        i+=1
    return contactsPerFrame # list of all the contacts in a frame, which can be processed separately

def checkContact(a1,a2,cutoff):
    """given coordinates of two particles a1 and a2, check if they are in contact given the cutoff"""

    distance=np.linalg.norm(a1-a2)
    if distance<=cutoff:
        return True
    else:
        return False

# mean first passage time - number of frames before a specific subcomplex forms
import numpy as np
import pandas as pd
import MDAnalysis
from itertools import groupby,product
import seaborn as sns
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from sklearn.cluster import DBSCAN
import networkx
from networkx.algorithms.components.connected import connected_components
import scipy.special

def getScore(i):
    f=open('%s.xml'%(i),'r')
    lines=f.readlines()
    f.close()
    score=0
    for l in lines:
        if 'score' in l:
            score+=float(l.split('"')[1])
    return score

def removeDuplicates(lst):
    return [t for t in (set(tuple(i) for i in lst))]

def groupCompare(g1,g2):
    list1=[]
    list2=[]
    for tup1 in g1:
        list1.append(sorted(tup1))
    for tup2 in g2:
        list2.append(sorted(tup2))
    if len(list1)==len(list2):
        print(list1)
        print(list2)
    if sorted(list1)==sorted(list2):
        return 1
    else:
        return 0
def sorted_k_partitions(seq,k):
    """Returns a list of all unique k-partitions of `seq`.

    Each partition is a list of parts, and each part is a tuple.

    The parts in each individual partition will be sorted in shortlex
    order (i.e., by length first, then lexicographically).

    The overall list of partitions will then be sorted by the length
    of their first part, the length of their second part, ...,
    the length of their last part, and then lexicographically.
    """
    n=len(seq)
    groups=[]  # a list of lists, currently empty

    def generate_partitions(i):
        if i>=n:
            yield list(map(tuple,groups))
        else:
            if n-i>k-len(groups):
                for group in groups:
                    group.append(seq[i])
                    yield from generate_partitions(i+1)
                    group.pop()

            if len(groups)<k:
                groups.append([seq[i]])
                yield from generate_partitions(i+1)
                groups.pop()

    result=generate_partitions(0)

    # Sort the parts in each partition in shortlex order
    result=[sorted(ps,key=lambda p: (len(p),p)) for ps in result]
    # Sort partitions by the length of each part, then lexicographically.
    result=sorted(result,key=lambda ps: (*map(len,ps), ps))

    return result
def Manhattan(tup1,tup2):
    return np.linalg.norm(np.asarray(tup1)-np.asarray(tup2))

def toGraph(l):
    G=networkx.Graph()
    for part in l:
        # each sublist is a bunch of nodes
        G.add_nodes_from(part)
        # it also imlies a number of edges:
        G.add_edges_from(toEdges(part))
    return G

def toEdges(l):
    """ 
        treat `l` as a Graph and returns it's edges 
        to_edges(['a','b','c','d']) -> [(a,b), (b,c),(c,d)]
    """
    it=iter(l)
    last=next(it)

    for current in it:
        yield last,current
        last=current


# OTHER #
# these are functions or code snippets that were being used in efforts to do certain analyses and which might be useful, but haven't been fitted successfully into things yet
