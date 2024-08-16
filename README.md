# OVERVIEW #
Set of Python codes for using IMP to set up and run Go-type models of macromolecular complex assembly, and carrying out a preliminary analysis on these


# GO SIMULATION #

`python runGoSimulation.py
     -f [output filename]
     -p [number of interaction particles per particle]
     -n [number of main particles]
     -i [input interaction file]
     -b [box dimensions in Angstroms]
     -r [particle radius in Angstroms]`


- This is a generalised script designed to allow someone to make a Go-type model of macromolecular complex assembly of any arbitrary configuration they want
- Each component is represented using a particle with or without interaction particles on its surface
-  The assembly is simulated using Brownian Dynamics in IMP
-  Input interaction file
     - This should just have the interacting pairs of particles given as indices, commas between the particles and a newline between each pair
     - If no interaction file is given, the script will make a ring structure by default
- The interaction particles are small particles on the surface of the main one, and at a fixed distance away from each other and the centre of the main particle
     - Where present, it will be these, and not the main particles, that participate in the direct interactions
     - Not necessary for more compact complexes, but are required to form stable ring or fibril structures
     - The distance of the particles relative to each other has been chosen to suit the ring and fibril structures I was working on, and might need modifying for other complex types; ideally, would have a way to determine this mathematically before the simulation as part of the system set up
- The script sets up a model with particles of equal radius and mass, defines truncated harmonic interactions between specific pairs, and runs a Brownian Dynamics simulation during which the particles assemble into the desired structure
     - Truncated harmonics were required to prevent complexes from immediately coalescing, and means that you get stepwise formation of subcomplexes


# SIMULATION ANALYSIS #

`python runGoAnalysis.py 
     -f [output filename (also name of input rmf file)]
     -p [number of interaction particles per particle]
     -n [number of main particles]
     -i [input interaction file]
     -r [particle radius in Angstroms]`


- This analysis script does some format conversions and basic analyses that attempt to get information about energetics and pathways of assembly, and produce a few simple plots, as well as storing Panda dataframes of the data in pickle format
- First converts the rmf files into something compatible with analysis using things like MDAnalysis and visualisation using Pymol or VMD – creates a trajectory (.xtc format)
     - The way this is done is quite clunky, first making xml files that are then turned into pdb files, which are combined into a trajectory - messy, but it works
- Then loop through the trajectory and collect various bits of data
- Then, for each frame
     - finds those particles that are close enough together to be considered in physical contact with each other, based on double the radii of the particles, getting the contacts and fraction of native contacts
     - calculates the score -- proxy for energy
     - stores a list of the unique complexes i.e. all the different subcomplexes that are seen, as well as the frame at which they first appear -- this is to get mean passage time when compared across simulations
- Plots score score, fraction of native contacts and number of contacts against frame
- Saves the scores, number of contacts, and fraction of native contacts as a Pandas dataframe, and the pathway (route between subcomplexes) as a separate dataframe
- Creates a lineplot that follows the assembly process visually – each particle is shown by a single coloured line that descends in the y-direction as the simulation goes on, and moves in the x-direction to be adjacent to another line when those particles interact
     - The lineplots that are produced are what I have found most visually useful to demonstrating that different repeats follow different paths for complex formation – I wanted to find a numerical way to demonstrate this, but have not managed to figure out how to do so; it is a useful result because experiments on in vitro ribosome assembly have shown the use of multiple different pathways


# BATCH PROCESSES #

`python batchGoSimuluations.py
     -f [output filename]
     -p [number of interaction particles per particle]
     -n [number of main particles]
     -i [input interaction file]
     -b [box dimensions in Angstroms]
     -r [particle radius in Angstroms]
     -s [number of simultions to run]`

- This runs Go simulations with the given input parameters for s number of repeats, putting the results from each in a separate directory, named according to the number of the simulation

`python batchGoAnalyses.py
     -f [output filename]
     -p [number of interaction particles per particle]
     -n [number of main particles]
     -i [input interaction file]
     -b [box dimensions in Angstroms]
     -r [particle radius in Angstroms]
     -s [number of simultions to run]`

- This runs analyses of the repeats and then combines some of the data to give mean passage time for each subcomplex, number of subcomplexes seen overall and what fraction this is for all possible subcomplexes (for number of particles up to ten), and a sort of energy landscape for the transition
- All of these need work to be more useful and informative
