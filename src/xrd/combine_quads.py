#!/usr/bin/env python
# coding: utf-8


import numpy as np
from scipy.interpolate import interp1d
from scipy.signal import find_peaks
import re
import sys
from lmfit import Parameters, fit_report, minimize
import matplotlib.pyplot as plt

def residual(pars, y, data=None, eps=None):
    # unpack parameters: extract .value attribute for each parameter
    parvals = pars.valuesdict()
    scale = parvals['scale']

    model = scale * y  

    if data is None:
        return model
    if eps is None:
        return (model - data)
    return (model-data) / eps

def merge_quads(quad_a, quad_b, plot=False):
    ''' 
    merge two quads and scale teh insensity so that the intensity profile is continuous
    Input:
    quad_a:
    quad_b:
    plot: bool ; should the output be plotted 
    '''
    if np.min(quad_b[0]) < np.min(quad_a[0]):
        tmp = quad_a
        quad_a = quad_b
        quad_b = tmp

    print(f'Quad_A going from : {quad_a[0].min()} to {quad_a[0].max()} ')
    print(f'Quad_B going from : {quad_b[0].min()} to {quad_b[0].max()} ')

    # Crop both arrays because the edges are never nice due to low statistics and saturation 
    quart_distance      = abs(np.max(quad_a[0])-np.min(quad_b[0]))/3
    id_A                = quad_a[0] < np.max(quad_a[0]) - quart_distance
    id_B                = quad_b[0] > np.min(quad_b[0]) + quart_distance
    quad_a              = np.array([quad_a[0][id_A],quad_a[1][id_A]])
    quad_b              = np.array([quad_b[0][id_B],quad_b[1][id_B]])

    # find the range between which we have data for both quads

    low, high = np.min(quad_b[0]), np.max(quad_a[0])
    overlap_a, overlap_b = np.where(quad_a[0] >= low)[0], np.where(quad_b[0] <= high)[0]

    # interpolate the two quads
    f_a, f_b = interp1d(quad_a[0], quad_a[1], kind='linear'), interp1d(quad_b[0], quad_b[1], kind='linear')
    
    # define the two different spacings in theta - not both quads are (nessecarily) integrated with the same theta-spacings - 
    dT_a, dT_b          = np.mean(quad_a[0][1::2] - quad_a[0][0:-(len(quad_a)%2+1):2]), np.mean(quad_b[0][1::2] - quad_b[0][0:-(len(quad_b)%2+1):2])
    pure_a, pure_b      = np.ones_like(quad_a[0],dtype=bool), np.ones_like(quad_b[0],dtype=bool)
    pure_a[overlap_a]   = 0 
    pure_b[overlap_b]   = 0

    # divide the regions where we have both data into theta-steps that are the mean of the two quads
    theta_low, I_low    = quad_a[0][np.where(pure_a)], quad_a[1][np.where(pure_a)]
    theta_mid           = np.arange(low, high, 0.5 *(dT_a+dT_b))
    theta_high, I_high  = quad_b[0][np.where(pure_b)], quad_b[1][np.where(pure_b)]
    
    # Fit the linear scale - this is more complicatedthen nessecary but it is also a nice exercise for minimisation
    fit_params          = Parameters()
    fit_params.add('scale', value=0.9)
    out                 = minimize(residual, fit_params, args=(f_b(theta_mid),), kws={'data': f_a(theta_mid)});
    param_dict          = out.params.valuesdict()
    scale               = param_dict['scale']

    I_high              *= scale
    I_mid               = (f_b(theta_mid) * scale + f_a(theta_mid)) / 2.
    
    quads_merged = np.array([np.concatenate([theta_low,theta_mid,theta_high]), np.concatenate([I_low, I_mid, I_high])] )

    if plot==True:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(1,3)
        ax[0].plot(quad_a[0],quad_a[1],label='Quad A')
        ax[1].plot(quad_b[0],quad_b[1],label='Quad B')
        ax[2].plot(quads_merged[0],quads_merged[1],label='Merged Quads')
        plt.show()
                     
    return quads_merged


def merge_xy_files(a,b,filename=None,plot=True):
    if not filename:
        print("No filename given!")
        filename = re.sub(".xy", "_merged.xy", a)

    quad3 = np.loadtxt(a).T
    quad2 = np.loadtxt(b).T
    quads_merged  = merge_quads(quad2, quad3, plot) 

    np.savetxt(filename, quads_merged.T)

def peaks_from_file(file_peak,theta_min=0.,theta_max=90.):
    peak_data           = np.loadtxt(file_peak,skiprows=1)
    two_theta_vals      = peak_data[:,7]
    peak_positions      = np.array(list(set(two_theta_vals)))
    multiplicity        = np.array([len(np.where(two_theta_vals == theta)[0]) for theta in peak_positions])
    rel_intensity       = np.array([peak_data[np.where(two_theta_vals == i)[0][0],8] for i in peak_positions])
    scaled_int          = 0.01 * (rel_intensity * multiplicity) / multiplicity[np.argmax(rel_intensity)]
    mask                = np.where(np.logical_and(peak_positions>=theta_min, peak_positions<=theta_max))
    return peak_positions[mask], scaled_int[mask]

def combine_quads_cif(files_quad, file_peak):
    num_quads       = len(files_quad)
    # define a list of the quad datas [Theta,I]
    list_quads_tmp      = [[]] * num_quads
    theta_start = np.zeros(num_quads)

    for quad_id, file in enumerate(files_quad):
        list_quads_tmp[quad_id]         = np.loadtxt(file).T
        theta_start[quad_id]            = np.min(list_quads_tmp[quad_id][0])
    
    sort = np.argsort(theta_start)
    list_quads = [list_quads_tmp[i] for i in sort]


    # define the theta coverage for each quad
    theta_range     = np.zeros([num_quads,2],dtype=np.float32)

    I_exp_ref                = 0.       # The highest intensity of all quads is used as the experimental reference intensity
    max_quad_Id                   = int(0)

    for quad_id, quad in enumerate(list_quads):
        theta_range[quad_id,0]      = np.min(quad[0])
        theta_range[quad_id,1]      = np.max(quad[0])
        if np.max(quad[1]) > I_exp_ref:
            I_exp_ref          = np.max(quad[1])
            max_quad_Id             = quad_id
        
    theta_min       = np.min(theta_range[:,0])
    theta_max       = np.max(theta_range[:,1])


    
    peak_positions, scaled_int = peaks_from_file(file_peak, theta_min, theta_max)
    # define a list of the ID of the peaks with the maximal theoretical intesinty in each quad
    max_peak_ids    = np.zeros(num_quads,dtype=int)
    peak_masks = [[]] * 4 # all peaks visible in each quad

    for quad_id, quad in enumerate(list_quads):
        peak_mask                 = np.where(np.logical_and(peak_positions>=theta_range[quad_id,0], peak_positions<=theta_range[quad_id,1])) # subset of peaks in the quad
        peak_masks[quad_id]       = peak_mask
        max_int_peak              = peak_positions[peak_mask][np.argmax(scaled_int[peak_mask])]                                              # the theta value of the peak with the max. theoretical I in the quad
        max_peak_ids[quad_id]     = np.where(peak_positions == max_int_peak)[0]
    
    total_max_peak                = max_peak_ids[np.argmax(scaled_int[max_peak_ids])] # this is the peak with the highest theoretical I in all quad

    # theoretical intensity of the refence peak 
    I_theo_ref =  scaled_int[max_peak_ids[max_quad_Id]]
    print(f"The highest Intensity is in quad {max_quad_Id} for peak {max_peak_ids[max_quad_Id]} ")

    scale = [[]] * 4
    # Go through all quads and set a scale factor so that the intensities match the expectation
    for quad_id, quad in enumerate(list_quads):
        if quad_id == max_quad_Id:
            scale[quad_id]      = 1.
            continue
        I_theo_quad             = scaled_int[max_peak_ids[quad_id]]
        scale[quad_id]          =  (I_exp_ref / np.max(quad[1])) * (I_theo_quad / I_theo_ref)


    tmp_1, tmp_2        = np.copy(list_quads[0]), np.copy(list_quads[1])
    tmp_1[1]            *= scale[0]
    tmp_2[1]            *= scale[1]

    print(f"Combined data: {tmp_1.shape} tmp: {tmp_2.shape}")
    combined_data       = merge_quads(tmp_1,tmp_2)

    for step in range(num_quads - 2):
        tmp             = np.copy(list_quads[step+2])
        tmp[1]             *= scale[step+2]
        print(f"Combined data: {combined_data.shape} tmp: {tmp.shape}")
        combined_data   = merge_quads(combined_data, tmp)    
        

    print(f"The {num_quads} Quads span a theta range from: {theta_min} to {theta_max}")
    print(f"Scaling factors are: {scale}")
    #fig, ax = plt.subplots()
    #for quad_id, quad in enumerate(list_quads):
    #    ax.plot(quad[0],quad[1]*scale[quad_id],label=f"Quad {quad_id}")
    #ax.scatter(peak_positions,scaled_int*I_exp_ref, label="Theoretical peak intensities")
    #ax.plot(combined_data[0],combined_data[1]+70,label=f"Combined data")
    #fig.legend()
    #plt.show()
    return combined_data, scale

def merge_four_quads(quad_data, scale):
    num_quads       = len(quad_data)
    # define a list of the quad datas [Theta,I]
    list_quads_tmp      = [[]] * num_quads
    theta_start = np.zeros(num_quads)

    for quad_id, data in enumerate(quad_data):
        list_quads_tmp[quad_id]         = data
        theta_start[quad_id]            = np.min(data[0])
    
    sort = np.argsort(theta_start)
    list_quads = [list_quads_tmp[i] for i in sort]

    [print(f"Shape of quad {quad_id} is {data.shape}") for quad_id, data in enumerate(list_quads)]

    tmp_1, tmp_2        = np.copy(list_quads[0]), np.copy(list_quads[1])
    tmp_1[1]            *= scale[0]
    tmp_2[1]            *= scale[1]

    print(f"Combined data: {tmp_1.shape} tmp: {tmp_2.shape}")
    combined_data       = merge_quads(tmp_1,tmp_2)

    for step in range(num_quads - 2):
        tmp             = np.copy(list_quads[step+2])
        tmp[1]             *= scale[step+2]
        print(f"Combined data: {combined_data.shape} tmp: {tmp.shape}")
        combined_data   = merge_quads(combined_data, tmp)    
        
    return combined_data

def combine_patterns(patterns):
    x_min = []
    for pattern in patterns:
        x = pattern[0]
        x_min.append(np.min(x))

    sorted_pattern_ind = np.argsort(x_min)

    pattern = patterns[sorted_pattern_ind[0]]
    for ind in sorted_pattern_ind[1:]:
        x1, y1 = pattern[0] , pattern[1] 
        x2, y2 = patterns[ind][0], patterns[ind][0]

        pattern2_interp1d = interp1d(x2, y2, kind='linear')

        overlap_ind_pattern1 = np.where((x1 <= np.max(x2)) & (x1 >= np.min(x2)))[0]
        left_ind_pattern1 = np.where((x1 <= np.min(x2)))[0]
        right_ind_pattern2 = np.where((x2 >= np.max(x1)))[0]

        combined_x1 = x1[left_ind_pattern1]
        combined_y1 = y1[left_ind_pattern1]
        combined_x2 = x1[overlap_ind_pattern1]
        combined_y2 = (y1[overlap_ind_pattern1] + pattern2_interp1d(combined_x2)) / 2
        combined_x3 = x2[right_ind_pattern2]
        combined_y3 = y2[right_ind_pattern2]

        combined_x = np.hstack((combined_x1, combined_x2, combined_x3))
        combined_y = np.hstack((combined_y1, combined_y2, combined_y3))

        pattern = np.vstack([combined_x, combined_y])

    return pattern


if __name__=="__main__":
    files_quad = [sys.argv[1],sys.argv[2], sys.argv[3],sys.argv[4]]
    file_peak = sys.argv[5]
    combined_data, scale  = combine_quads_cif(files_quad, file_peak)
    np.savetxt(sys.argv[6], combined_data.T)
    np.savetxt("quad_scale_LW03.txt", scale)


    #merge_xy_files(sys.argv[1],sys.argv[2], sys.argv[3])
