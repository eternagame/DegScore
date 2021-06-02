import numpy as np
import sys, os

try:
    from arnie.mfe import mfe                   
except ImportError:
    print('Warning: could not find Arnie for DegScore structure prediction.\n\
        Secondary structures must be input in DegScore class.')

from assign_loop_type import write_loop_assignments
DEBUG=False

coeffs_2_1 = [-0.020, -0.027, -0.026, -0.017, 0.005, 0.000, 0.005, -0.006, -0.011, 0.006, 0.031,
 0.021, 0.036, 0.034, 0.026, -0.024, -0.005, 0.028, -0.022, -0.015, -0.043, -0.043, -0.029, -0.026,
  0.016, -0.077, -0.001, -0.016, 0.031, -0.001, 0.064, 0.065, 0.064, 0.069, 0.029, 0.044, -0.003,
  0.012, -0.006, -0.004, -0.072, -0.066, -0.061, -0.065, 0.014, 0.037, 0.051, 0.017, 0.054, 0.037,
 -0.065, -0.068, -0.058, -0.041, -0.014, 0.075, -0.007, 0.005, -0.010, -0.006, 0.009, 0.014, 0.019,
  0.037, 0.005, -0.097, 0.013, -0.005, 0.001, 0.002, -0.026, -0.026, -0.036, -0.008, 0.041, 0.067, 0.017,
   0.007, 0.034, 0.028, -0.077, -0.079, -0.092, -0.064, 0.012, 0.022, 0.041, 0.041, 0.057, 0.038, 0.010,
    0.017, -0.004, 0.050, -0.018, -0.144, -0.001, 0.009, 0.013, -0.017, -0.012, -0.019, -0.047, -0.032,
     0.124, 0.164, 0.089, 0.063, 0.076, 0.055, 0.012, 0.021, -0.038, -0.050, 0.014, -0.044, -0.023,
      -0.015, -0.037, -0.059, -0.027, 0.042, 0.003, -0.017, -0.090, -0.057, -0.140, -0.005, -0.031,
 -0.256, -0.353, -0.178, -0.503, -0.506, -0.090, -0.058, -0.032, -0.035, 0.015, -0.070,
0.048, 0.002, -0.071, 0.005, -0.002, 0.024, 0.000, -0.007, -0.014, -0.029, 0.033, 0.011,
 0.006, -0.039, -0.023, 0.046, 0.009, 0.009, 0.003, -0.002, 0.015, 0.022, 0.007, 0.015,
  -0.003, -0.121, -0.014, -0.021, -0.016, -0.010, -0.046, -0.059, -0.031, -0.013, 0.005,
   0.051, -0.006, 0.003, 0.009, -0.006, 0.001, 0.021, 0.025, 0.017, 0.012, 0.033, -0.003,
-0.030, 0.002, 0.008, -0.003, -0.001, 0.003, 0.007, 0.004, -0.099, -0.003, 0.006,
 -0.018, -0.002, -0.018, -0.023, -0.016, -0.008, -0.013, -0.027, -0.006, 0.005,
  0.011, -0.001, 0.035, 0.024, 0.042, 0.040, -0.010, 0.270, 0.006, -0.040, -0.001,
  -0.012, -0.042, -0.040, -0.015, -0.021, 0.008, -0.294, -0.005, -0.008, 0.001,
 -0.014, -0.037, -0.047, -0.034, -0.028, 0.025, 0.144, 0.016, 0.023, 0.026, 0.015,
  0.014, 0.009, 0.020, 0.031, 0.026, -0.096, 0.002, -0.012, 0.030, 0.000]

k_deg_m, k_deg_b = 0.002170959651184987, 0.05220886935630193

def encode_input(sequence, bprna_string, window_size=12, pad=0, seq=True, struct=True):
    '''Creat input/output for regression model for predicting structure probing data.
    Inputs:
    
    sequence (str): RNA sequence
    bprna_string (str): loop assignment string (HEIMBS)
    window_size: size of window (in one direction). so window_size=1 is a total window size of 3.
    pad: number of nucleotides at start to not include
    seq (bool): include sequence encoding
    struct (bool): include bpRNA structure encoding
    
    Outputs:
    Input array (n_samples x n_features): array of windowed input features
    feature_names (list): feature names
    '''

    assert len(sequence) == len(bprna_string)

    feature_kernel=[]
    if seq:
        feature_kernel.extend(['A','U','G','C'])
    if struct:
        feature_kernel.extend(['H','E','I','M','B','S'])
    
    inpts = []
    labels = []
    
    #for i, row in df.iterrows():
        
    length = len(sequence)
    
    arr = np.zeros([length,len(feature_kernel)])

    for index in range(length):
        ctr=0

        #encode sequence
        if seq:
            for char in ['A','U','G','C']:
                if sequence[index]==char:
                    arr[index,ctr]+=1
                ctr+=1

        if struct:
            for char in ['H','E','I','M','B','S']:
                if bprna_string[index]==char:
                    arr[index,ctr]+=1
                ctr+=1

    # add zero padding to the side

    padded_arr = np.vstack([np.zeros([window_size,len(feature_kernel)]), arr, np.zeros([window_size,len(feature_kernel)])])

    for index in range(length):
        new_index = index+window_size-pad
        tmp = padded_arr[new_index-window_size:new_index+window_size+1]
        inpts.append(tmp.flatten())
            
    return np.array(inpts)

def create_U_mask(seq):
    mask=[]
    for i, x in enumerate(list(seq)):
        if x!='U':
            mask.append(i)
    return np.array(mask)

class DegScore():
    def __init__(self, sequence, structure=None, mask_U=False, package='eternafold',
        start_ind=None, end_ind=None, linear=False, coefficients=None, intercept=None):
        '''Class to calculate DegScore-2.1, a ridge regression model to predict degradation.
        H Wayment-Steele, 2020/2021.

        Inputs:
        sequence (str): RNA sequence
        structure (str): RNA dot-bracket structure. If not provided and Arnie is provided,
        will re-calculate based provided 'package' and 'linear' keywords.
        mask_U (bool, default False): If True, sets U positions to zero to mimic pseudouridine stabilization.
        start_ind: starting position to sum degscore (default 0).
        end_ind: ending position to sum degscore (defuault len(sequence)).

        Structure prediction options:
        package: package to use to calculate MFE secondary structure (example options in Arnie: 'vienna', 'eternafold')
        linear (bool): Use linearfold calculation (must be set up in Arnie.)

        Regression options:
        coefficients (list) (optional): coefficients from Ridge regression. Default is DegScore 2.1 coefficients.
        intercept (float) (optional): intercept from Ridge regression. Default is DegScore 2.1 intercept.

        Attributes:
        loop_assignments (str): loop assignments:
                    H: Hairpin, E: External, S: Stem, I: Internal, B: Bulge, M: Multiloop	

		degscore: DegScore (float), summed across all nucleotides.
        degscore_by_position (vector): DegScore at each position.
        est_k_deg (float): Estimated degradation rate (hrs^-1).
        est_half_life (float): Estimated half life (hrs).

        '''
        self.sequence = sequence
        if structure is not None:
            self.structure = structure
        else:
            self.structure = mfe(sequence, package=package, linear=linear, viterbi=True)

        assert len(self.sequence) == len(self.structure)

        self.loop_assignments = write_loop_assignments(self.structure)

        if coefficients is None:
            self.coefficients_ =  coeffs_2_1
        else:
            self.coefficients_ = coefficients
        if intercept is None:
            self.intercept_ = 1.122
        else:
            self.intercept_ = intercept

        if DEBUG: print(self.bprna_string)

        self.encoding_ = encode_input(self.sequence, self.loop_assignments)

        if DEBUG: print("encoding shape", self.encoding_.shape)

        self.degscore_by_position = np.sum(self.encoding_ * self.coefficients_, axis=1) + self.intercept_

        if mask_U:
            mask_inds = create_U_mask(self.sequence)

            mask = np.ones(self.degscore_by_position.size, dtype=bool)
            mask[mask_inds] = False
            self.degscore_by_position[mask] = 0

        if start_ind is None:
            start_ind = 0
        if end_ind is None:
            end_ind = len(sequence)

        self.degscore = np.sum([self.degscore_by_position[x] for x in range(start_ind, end_ind)])

        self.est_k_deg = k_deg_m*self.degscore + k_deg_b

        self.est_half_life = np.log(2)/self.est_k_deg
