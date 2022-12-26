from ovito.io import import_file
from ovito.modifiers import IdentifyDiamondModifier, CoordinationAnalysisModifier, CreateBondsModifier
from ovito.data import *
from tqdm import tqdm
import h5py as h5
import re
import pandas as pd
import logging
import os
import sys
import numpy as np
import time
from load_lammps_log import load_thermo

logger                  = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def load_files(basename, outfile, bins=200, hdf_dir='./'):
    dumpfile, logfile = basename + '.dump', basename + '.log'
    try:
        pipeline          = import_file(dumpfile)
    except:
        logger.error(f"Dumpfile {dumpfile} not found!")
        return
    try:
        dfs, dt, n_atoms, run_steps = load_thermo(logfile)
    except:
        logger.error(f"Logfile {logfile} not found!")
        return

    logger.info(f"Dumpfile {dumpfile} loaded.")

    hdf_file                = hdf_dir + re.split('/',basename)[-1] 
    hdf_file                = outfile

    pipeline.modifiers.append(IdentifyDiamondModifier())
    pipeline.modifiers.append(CoordinationAnalysisModifier(cutoff = 8.0, number_of_bins = bins))
    pipeline.modifiers.append(CreateBondsModifier(cutoff = 1.704))

    data                    = pipeline.compute(0)
    num_frames              = pipeline.source.num_frames
    num_atoms               = data.particles.count
    print(f'Found {num_frames} frames with {num_atoms} atoms each!')

    diamond_lattice         = np.zeros((num_frames,7),dtype=int)
    structure_type          = np.zeros((num_frames, num_atoms),dtype=int)
    rdf_total               = np.zeros((num_frames,bins,2),dtype=float)
    #diamond_center          = np.zeros((num_frames,3),dtype=np.double)
    #diamond_radius          = np.zeros((num_frames),dtype=np.double)
    #diamond_onion           = np.zeros((num_frames,bins,7),dtype=float)
    #diamond_onion_outer_r   = np.zeros((num_frames,bins),dtype=float)
    neighbours              = np.zeros((num_frames,num_atoms),dtype=int)
    timesteps               = np.zeros(num_frames,dtype=int)

    if not os.path.exists(hdf_file):
        f                           = h5.File(hdf_file, "w")

        equilibration_thermo        = f.create_group('ThermoData/equilibration',track_order=True)
        expansion_thermo            = f.create_group('ThermoData/expansion',track_order=True)
        structureAnalysis           = f.create_group('StructureAnalysis',track_order=True)
        structureAnalysis.attrs.create('n_frames',data=num_frames)
        structureAnalysis.attrs.create('processed_frames',data=0)
        
        sets                        = [equilibration_thermo, expansion_thermo]

        for purpose in [0,1]:
            sets[purpose].attrs.create('n_atoms', data=n_atoms)
            sets[purpose].attrs.create('n_steps', data=run_steps[purpose])

            for key in dfs[purpose].columns:
                sets[purpose].create_dataset(key, data=dfs[purpose][key])

        structureAnalysis.create_dataset('StructureType', data=structure_type)
        structureAnalysis.create_dataset('NeighbourCount',data=neighbours)
        structureAnalysis.create_dataset('Timesteps',data=timesteps)
        structureAnalysis.create_dataset('RadialDistributionFunction',data=rdf_total)

        f.close()


    tic             = time.perf_counter()
    start_frame     = 0
    stepsize        = 500

    for _ in range(num_frames//stepsize + 1):
    #for _ in range(1):
        f                           = h5.File(hdf_file, "a")
        start_frame     = f['StructureAnalysis'].attrs['processed_frames']
        end_frame       = start_frame + stepsize
        #for frame in tqdm(range(f['StructureAnalysis'].attrs['processed_frames'] ,f['StructureAnalysis'].attrs['processed_frames'] + stepsize)):
        for frame in tqdm(range(start_frame,end_frame)):
            if frame >= num_frames:
                print(f"{frame} Done")
                break
            diamond_lattice         = np.zeros((7),dtype=int)
            structure_type          = np.zeros((num_atoms),dtype=int)
            rdf_total               = np.zeros((bins,2),dtype=float)
            neighbours              = np.zeros((num_atoms),dtype=int)

            data                    = pipeline.compute(frame)
            f['StructureAnalysis/RadialDistributionFunction'][frame]        = data.tables['coordination-rdf'].xy()
            f['StructureAnalysis/Timesteps'][frame]                         = data.attributes['Timestep']

            ### if atoms are lost -> leave the rest as zero
            structure_type[:len(data.particles.structure_type.array)] = np.array(data.particles.structure_type.array)

            bonds_enum              = BondsEnumerator(data.particles.bonds)
            
            for i in range(7):
                count                       = len(np.where(structure_type == i)[0])
                diamond_lattice[i]    = count

            for particle_index in range(num_atoms):
                for _ in bonds_enum.bonds_of_particle(particle_index):
                    neighbours[particle_index] += 1

            f['StructureAnalysis/StructureType'][frame]     = structure_type
            f['StructureAnalysis/NeighbourCount'][frame]    = neighbours
            #print(f'Frame {frame} structure types are {np.array(data.particles.structure_type.array)[:20]}') 
            ### structure type 0: unordered, 1: diamond, 2: diamond first neighbour, 3: diamond second neihbour ###
            #if diamond_lattice[frame,1:4].sum(for frame in tqdm(range(start_frame,max_frame)):
            
            ### structure type 0: unordered, 1: diamond, 2: diamond first neighbour, 3: diamond second neihbour ###
            #if diamond_lattice[frame,1:4].sum() == 0:
            #    diamond_center[frame]           = np.array([0.,0.,0.])
            #    diamond_radius[frame]           = 0.

        f['StructureAnalysis'].attrs['processed_frames']         =   f['StructureAnalysis'].attrs['processed_frames'] + stepsize
        #f['StructureAnalysis'].attrs['processed_frames']         =   end_frame
        #processedFrames         = f['StructureAnalysis'].attrs['processed_frames']
        processedFrames         = end_frame - start_frame
        toc = time.perf_counter()
        logger.info(f'Calculated all neighbours and structure type for frames {start_frame} to {end_frame} in {toc - tic:0.4f} seconds.\n{(toc - tic)/(processedFrames):0.4f} seconds per frame.')
        f.close()

if __name__ == "__main__":
    logger.setLevel(logging.INFO)
    logging.basicConfig(level=logging.INFO)
    if len(sys.argv) == 0:
        logger.error('Please provide a basefile path and name')
    load_files(sys.argv[1],sys.argv[2])