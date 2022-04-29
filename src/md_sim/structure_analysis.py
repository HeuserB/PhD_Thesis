# Import OVITO modules.
from ovito.io import import_file
from ovito.modifiers import IdentifyDiamondModifier, CoordinationAnalysisModifier
from ovito.data import *

import numpy as np
NA = np.newaxis

def split_onions(r, lattice, cell, bins, r_max):
    r_array             = np.array(r)
    diamond_idx         = np.where(lattice == 1)[0]
    # find the center of the diamond structure
    if len(diamond_idx) != 0:
        r_0 = np.mean(r_array[diamond_idx],axis=0)
        #r_max = np.min(np.sqrt(np.sum((cell - np.repeat(r_0[:,NA],4,axis=1))**2,axis=1)))
    else:
        r_0 = np.mean(r_array,axis=0)
        #r_max = np.min(np.sqrt(np.sum((cell - np.repeat(r_0[:,NA],4,axis=1))**2,axis=1)))
    
    r_inner             = np.arange(bins) * r_max/bins
    r_outer             = np.arange(bins) * r_max/bins + r_max/bins
    
    lattice_onions = np.zeros((bins,7),dtype=float)
    
    for r_id in range(bins):
        # look for atoms in onion
        #int("r_Min: %.3f ; r_Max: %.3f" %(r_inner[r_id], r_outer[r_id]))
        r_rel           = np.abs(np.sqrt(np.sum((r-r_0)**2,axis=1)))
        atoms_in        = np.where((r_inner[r_id]<=r_rel) & (r_rel<=r_outer[r_id]))[0]
        
        for lat in range(7):
            # count the four different types in the onion shell
            lattice_onions[r_id][lat]   = np.sum(lattice[atoms_in] == lat)
        
        # normalise the counts to the number of atoms in the onion shell
        lattice_onions[r_id]            = lattice_onions[r_id][:]/np.sum(lattice_onions[r_id][:])

    return lattice_onions, r_outer

def largest_diamond_onion(structure_type, lattice_onions_all, r_outer):
    num_frames          = structure_type.shape[0]
    largest_diamond_r   = np.zeros(num_frames)
    for frame in range(num_frames):
        r_id_max                    = np.where(lattice_onions_all[frame][:,3] > 0.)[0][-1]
        largest_diamond_r[frame]    = r_outer[r_id_max]
    return largest_diamond_r

def diamond_structure(file, bins =200):
    ### for detailed description see https://ovito.org/manual/reference/pipelines/modifiers/identify_diamond.html#particles-modifiers-identify-diamond-structure ###
    
    pipeline = import_file(file)
    pipeline.modifiers.append(IdentifyDiamondModifier())
    pipeline.modifiers.append(CoordinationAnalysisModifier(cutoff = 8.0, number_of_bins = bins))
    
    data = pipeline.compute(0)
    
    num_frames = pipeline.source.num_frames
    num_atoms = data.particles.count

    diamond_lattice = np.zeros((num_frames,7),dtype=int)
    structure_type = np.zeros((num_frames, num_atoms),dtype=int)
    rdf_total = np.zeros((num_frames,bins,2),dtype=float)
    lattice_onions_all = np.zeros((num_frames,bins,7),dtype=float)

    for frame in range(num_frames):
        #print('Processing frame %02d/%02d' %(frame,num_frames))
        data = pipeline.compute(frame)
        rdf_total[frame] = pipeline.compute(frame).tables['coordination-rdf'].xy()
        #rdf = pipeline.compute(frame).tables['coordination-rdf'].xy()
        #print(rdf)
        
        ### if atoms are lost -> leave the rest as zero
        structure_type[frame][:len(data.particles.structure_type.array)] = np.array(data.particles.structure_type.array)
        
        particles = np.array(data.particles.position[:])
        cell = np.array(data.cell[:])
        
        for i in range(7):
            count = len(np.where(structure_type[frame] == i)[0])
            diamond_lattice[frame,i] = count

        if (diamond_lattice[frame].sum() != len(structure_type[frame])):
            print("Atom numbers do not add up!")
        
        lattice_onions_all[frame], r_outer = split_onions(particles, structure_type[frame], cell, bins, r_max=60.0)
    
    where_are_NaNs = np.isnan(lattice_onions_all)
    lattice_onions_all[where_are_NaNs] = 0.

    #largest_diamond_r = largest_diamond_onion(structure_type, lattice_onions_all, r_outer)

    return structure_type, diamond_lattice, rdf_total, lattice_onions_all, r_outer, #largest_diamond_r