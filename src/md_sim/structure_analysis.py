# Import OVITO modules.
#from decimal import ROUND_05UP
#from tokenize import group
#from unicodedata import name
#from grpc import stream_stream_rpc_method_handler
#from turtle import position

from ovito.io import import_file
from ovito.modifiers import IdentifyDiamondModifier, CoordinationAnalysisModifier
from ovito.data import *
#from sklearn.neighbors import radius_neighbors_graph
from miniball import miniball
from tqdm import tqdm
import h5py as h5
import re
import pandas as pd
import logging
import os


import time
from load_lammps_log import load_all_runs

from multiprocessing import Pool

import numpy as np
NA = np.newaxis

# Create a custom logger
logger                  = logging.getLogger(__name__)
stream_handler          = logging.StreamHandler()

f_handler               = logging.FileHandler('file.log')
stream_handler.setLevel(logging.DEBUG)
f_handler.setLevel(logging.DEBUG)

stream_format           = logging.Formatter('%(name)s - %(levelname)s - %(message)s')
f_format                = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

stream_handler.setFormatter(stream_format)
f_handler.setFormatter(f_format)

logger.addHandler(stream_handler)
logger.addHandler(f_handler)

logger.setLevel(logging.DEBUG)


def split_onions(positions : np.array, structure_type : np.array , diamond_lattice : np.array , bins : int, r_max : float = 20.0):
    mask                            = (structure_type < 4) & (0 < structure_type)
    lattice_onions                  = np.zeros((bins,7),dtype=float)
    r_inner             = np.arange(bins) * r_max/bins
    r_outer             = np.arange(bins) * r_max/bins + r_max/bins

    # find the center of the diamond structure if there are atoms in a diamond lattice
    if diamond_lattice[1:4].sum() == 0:
        return lattice_onions, r_outer
    else:
        r_0 = np.mean(positions[mask],axis=0)
    
    for shell_id in range(bins):
        # look for atoms in onion
        #int("r_Min: %.3f ; r_Max: %.3f" %(r_inner[r_id], r_outer[r_id]))
        r_rel           = np.abs(np.sqrt(np.sum((positions-r_0)**2,axis=1)))
        atoms_in_shell  = np.where((r_inner[shell_id]<=r_rel) & (r_rel<=r_outer[shell_id]))[0]
        
        for lat in range(7):
            # count the four different types in the onion shell
            lattice_onions[shell_id][lat]   = np.sum(structure_type[atoms_in_shell] == lat)
        
        # normalise the counts to the number of atoms in the onion shell
        lattice_onions[shell_id]            = lattice_onions[shell_id][:]/np.sum(lattice_onions[shell_id][:])

    return lattice_onions, r_outer

def largest_diamond_onion(structure_type, lattice_onions_all, r_outer):
    num_frames          = structure_type.shape[0]
    largest_diamond_r   = np.zeros(num_frames)
    for frame in range(num_frames):
        r_id_max                    = np.where(lattice_onions_all[frame][:,3] > 0.)[0][-1]
        largest_diamond_r[frame]    = r_outer[r_id_max]
    return largest_diamond_r

def fit_spheres(points):
    res             = miniball(points)
    center          = res['center']
    radius          = res['radius']
    return [center, radius]

def diamond_structure(file, bins=200):
    ### for detailed description see https://ovito.org/manual/reference/pipelines/modifiers/identify_diamond.html#particles-modifiers-identify-diamond-structure ###
    
    pipeline                = import_file(file)
    pipeline.modifiers.append(IdentifyDiamondModifier())
    pipeline.modifiers.append(CoordinationAnalysisModifier(cutoff = 8.0, number_of_bins = bins))
    
    data                    = pipeline.compute(0)
    
    num_frames              = pipeline.source.num_frames
    num_atoms               = data.particles.count

    diamond_lattice         = np.zeros((num_frames,7),dtype=int)
    structure_type          = np.zeros((num_frames, num_atoms),dtype=int)
    rdf_total               = np.zeros((num_frames,bins,2),dtype=float)
    diamond_center          = np.zeros((num_frames,3),dtype=np.double)
    diamond_radius          = np.zeros((num_frames),dtype=np.double)
    diamond_onion           = np.zeros((num_frames,bins,7),dtype=np.float)
    diamond_onion_outer_r   = np.zeros((num_frames,bins),dtype=np.float)
    timesteps               = np.zeros(num_frames,dtype=int)

    for frame in tqdm(np.arange(num_frames)):
        #logger.debug(f'Processing frame {frame}/{num_frames}')
        data                = pipeline.compute(frame)
        rdf_total[frame]    = pipeline.compute(frame).tables['coordination-rdf'].xy()
        timesteps[frame]    = pipeline.compute().attributes['Timestep']
        
        ### if atoms are lost -> leave the rest as zero
        structure_type[frame][:len(data.particles.structure_type.array)] = np.array(data.particles.structure_type.array)
        
        particles           = np.array(data.particles.position[:])
        if particles.shape[0] != num_atoms:
            logger.debug(f'We seem to have lost some atoms in the simulation, putting the rest far away')
            new_particles                           = np.ones((num_atoms,3)) * 200.
            new_particles[:particles.shape[0],:]    = particles
            particles                               = new_particles

        for i in range(7):
            count                       = len(np.where(structure_type[frame] == i)[0])
            diamond_lattice[frame,i]    = count

        if (diamond_lattice[frame].sum() != len(structure_type[frame])):
            logger.warning("Atom numbers do not add up!")
            print("Atom numbers do not add up!")
        
        ### structure type 0: unordered, 1: diamond, 2: diamond first neighbour, 3: diamond second neihbour ###
        if diamond_lattice[frame,1:4].sum() == 0:
            diamond_center[frame]           = np.array([0.,0.,0.])
            diamond_radius[frame]           = 0.
            continue
        mask                            = (structure_type[frame] < 4) & (0 < structure_type[frame])
        
        #logger.debug(f'Created mask of shape: {np.array(mask).shape}')
        #logger.debug(f'Schape of particles is {particles.shape}')
        #tic                             = time.perf_counter()
        #logger.debug(f'Now fitting sphere dimension of masked particles is {np.array(particles[mask],dtype=np.double).shape}')
        center, radius                  = fit_spheres(np.array(particles[mask],dtype=np.double))
        #toc                             = time.perf_counter()
        #print(f"Fitting sphere took {toc - tic:0.4f} seconds")

        diamond_center[frame]           = center
        diamond_radius[frame]           = radius

        diamond_onion[frame], \
        diamond_onion_outer_r[frame]    = split_onions(particles, structure_type[frame], diamond_lattice[frame],bins)

    return structure_type, diamond_lattice, rdf_total, diamond_center, diamond_radius, diamond_onion, diamond_onion_outer_r, timesteps

#def store_as_hdf(dir : str, filename : str, bins : int, run_group : str = "airebo_m"):

def store_as_hdf(dir : str, filename : str, bins : int, run_group : str = "airebo_m"):
    run_data, thermo_dfs    = load_all_runs(dir)
    print(f"Found {len(run_data)} runs.")
    keys                    = ["structure_type", "diamond_lattice", "rdf_total", \
                                "diamond_center", "diamond_radius", "diamond_onion", \
                                "diamond_onion_outer_r", "timesteps"]
    columns                 = ['lattice_constant','temperature','cells','timestep_eq',\
                                'timestep_ex','drift','n_atoms','dE','dumpfile','run_steps_eq', 'run_steps_ex']
    written_files           = []

    for run_id in range(len(run_data)):
        #if run_data.loc[run_id]["dumpfile"] != '/mnt/TOSHIBA_EXT/Nextcloud/Work/W_PhD/W_PhD_Analysis/Simulations/LAMMPS/free_expansion_airebo_m/grid_3_35/2000K_8x8x8':
        #    print(f"Skipping file: {run_data['dumpfile'][run_id]}.")
        #    continue
        if run_data.loc[run_id]['cells'] != 8:
            logger.debug(f"Skipping file: {run_data['dumpfile'][run_id]} because it does not have 8x8x8 cells")
            print(f"Skipping file because it does not have 8x8x8 cells")
            continue
        tmp                     = re.sub("\.","_",str(run_data['lattice_constant'][run_id])) + "_" +re.split('/',run_data['dumpfile'][run_id])[-1]

        hdf_file                = dir + filename + '_' + tmp
        written_files.append(hdf_file + ".h5")
        f = h5.File(hdf_file + ".h5", "w")
        #main_group              = f.create_group(tmp,track_order=True)
        print(f"Processing file {run_id+1}/{len(run_data)}")
        logger.debug(f"Processing file {run_id+1}/{len(run_data)} with name {run_data['dumpfile'][run_id]}")
        dumpfile                = run_data['dumpfile'][run_id]
        run_data
        data                    = diamond_structure(dumpfile, bins)
        
        run_metadata            = f.create_group(tmp,track_order=True)
        
        main_attrs              = run_metadata.attrs
        for column in columns:
            logger.debug(f'Added column {column} to the main attributes')
            main_attrs.create(column, data=run_data[column].values[run_id])

        for key, set in zip(keys, data):
            hf_data             = run_metadata.create_dataset(key, data=set)

        for incentive_id, run_incentive in enumerate(['equilibration','expansion']):
            ### incentive is equilibration or expansion
            ### we want to write the data for both into a new dataframe
            thermo_df               = thermo_dfs[run_id][incentive_id]
            name                    = "thermo_" + run_incentive
            thermo_group            = run_metadata.create_group(name,track_order=True)
            for key in thermo_df.keys():
                thermo_data         = thermo_group.create_dataset(key, data=thermo_df[key].values)
        f.close()
        logger.debug(f'Written data to file {hdf_file + ".h5"}')
    print(f'Done with all runs. Combining hdf5 files.')
    combine_hdf5(dir + filename + '.h5', written_files, run_group)

def load_hdf(filename : str) -> [pd.DataFrame, dict] :
    """
    :param filename:             name of the *.h5 file with all the processed data
    :return:                     pandas dataframe with all run information, dict with the structure analysis
    """
    with h5.File(filename, "r") as f:
        group_key           = list(f.keys())[0]
        logger.debug(f'Found group key: {group_key}')
        group               = f[group_key]
        run_keys            = list(group.keys())

        #logger.debug(f'Found run keys: {run_keys}')
        
        #logger.debug(f'The metadata keys are: {metadata_keys}')
        #metadata            = np.array([group.attrs.__getitem__(attr) for attr in metadata_keys])
        #metadata_df         = pd.DataFrame(data=metadata.T, columns=metadata_keys)
        run_data_dic        = {}
        for key in run_keys:
            metadata_keys       = list(group[key].attrs.keys())
            metadata            = np.array([[group[key].attrs.__getitem__(attr) for attr in metadata_keys]])
            #logger.debug(f'metadata keys are: {metadata_keys} metadata shape ois {metadata}')
            metadata_df         = pd.DataFrame(data=metadata, columns=metadata_keys)
            # for each run define a subdirectory
            subdict         = {}
            data            = []
            keys            = []
            for subkey in group[key]:
                if type(group[key][subkey]) == h5.Group:
                    #logger.debug(f'Subkey {subkey} is an hdf5 group object!')
                    subsubdict      = {}
                    for subsubkey in group[key][subkey].keys():
                        subsubdict.update({subsubkey : group[key][subkey][subsubkey][:]})
                    subdict.update({subkey : subsubdict})
                    continue
                data.append(np.array(group[key][subkey][:]))
                keys.append(subkey)
                #logger.debug(f'Subkey is {subkey}')
                ## fill it with the values from the structure analysis
                subdict.update({subkey : group[key][subkey][:]})
            #for metakey in metadata_keys:
            #subdict.update({'run_data' : group[key][subkey][:]})
            subdict.update({'metadata' : metadata_df})
            ## and add it to the run directory 
            run_data_dic.update({key : subdict})
        
        return run_data_dic #metadata_df, run_data_dic

def combine_hdf5(outfile : str, h5files : list, run_group : str):
    with h5.File(outfile, "w") as f_dst:
        main_group      = f_dst.create_group(run_group)
        for i, filename in enumerate(h5files):
            print(filename)
            with h5.File(filename) as f_scr:
                logger.debug(f'Now reading structure analysis data for file {filename}')
                groupname       = list(f_scr.keys())[0]
                group           = f_scr[groupname]
                subgroup        = main_group.create_group(groupname)
                for obj in group.keys():
                    group.copy(obj, subgroup, name=obj)
                
                metadata_keys       = list(group.attrs.keys())
                #[print(f"type of attribute {attr} is: {type(group.attrs.__getitem__(attr))}") for attr in metadata_keys]
                #metadata            = np.array([group.attrs.__getitem__(attr) for attr in metadata_keys])
                #print(f'Metadata is: {metadata_keys} : {metadata}')
                main_attrs          = subgroup.attrs
                for id, key in enumerate(metadata_keys):
                    #print(f'data is: {key} : {type(metadata[id])}')
                    main_attrs.create(key, data=group.attrs.__getitem__(key))
    #cleanup(h5files)
    logger.debug(f'All files were written to the new HDF5 file {outfile}')

def cleanup(h5files):
    logger.debug(f'Cleaning up.')
    for file in h5files:
        os.remove(file)
    logger.debug(f'Done!')
