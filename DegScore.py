from arnie.utils import *
import numpy as np
from collections import Counter
package_locs = load_package_locations()
from arnie.mfe import mfe
DEBUG=False

'''
Usage with arnie:
add to arnie file:
bprna: /path/to/bpRNA
bprna is cloned from https://github.com/hendrixlab/bpRNA.
'''

coeffs = [-0.020, -0.027, -0.026, -0.017, 0.005, 0.000, 0.005, -0.006, -0.011, 0.006, 0.031,
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



def encode_input(sequence, bprna_string, window_size=12, pad=0, seq=True, struct=True):
    '''Creat input/output for regression model for predicting structure probing data.
    Inputs:
    
    sequence (in EternaBench RDAT format)
    window_size: size of window (in one direction). so window_size=1 is a total window size of 3
    pad: number of nucleotides at start to not include
    seq (bool): include sequence encoding
    struct (bool): include bpRNA structure encoding
    
    Outputs:
    Input array (n_samples x n_features): array of windowed input features
    feature_names (list): feature names
    
    '''

    assert len(sequence) == len(bprna_string)
    MAX_LEN = 1588

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


class DegScore():
    def __init__(self, sequence, structure=None, package='eternafold', linear=False):
        '''Class to handle DegScore information.
        Initialize by providing sequence and structure.
        Uses bpRNA to handle structure parsing.
        H: Hairpin, E: External, S: Stem, I: Internal, B: Bulge, M: Multiloop
        Attributes:
		counts: dictionary of counts of nucleotides in different secstruct features (hairpin, external, stem, internal, bulge, multiloop)
		weight: dictionary of floats for weight attributed to each secstruct feature
		score: DegScore (float), scores feature counts based on weights
		set_weights(weight_dictionary): update weights and re-score without re-running bpRNA.
        '''
        self.sequence = sequence
        if structure is not None:
            self.structure = structure
        else:
            self.structure = mfe(sequence, package=package, linear=linear, viterbi=True)
            
        assert len(sequence) == len(self.structure)

        fname = "%s.dbn" % filename()

        # old weights from DegScore 1
        #self.weights = {'H':0.7, 'E':1.0, 'S':0.0, 'I':0.2, 'B': 0.8, 'M': 0.6}
        self.coeffs_ =  coeffs #np.loadtxt('DegScore2.1_coeffs.txt',usecols=1)
        self.intercept_ = 1.122

        _ = write([sequence, self.structure],fname=fname)
        LOC=package_locs['bprna']
        command = ['perl', "%s/bpRNA.pl" % LOC, fname]
        if DEBUG:
            print(' '.join(command))
        p = sp.Popen(command, stdout=sp.PIPE, stderr=sp.PIPE)
        stdout, stderr = p.communicate()
        if DEBUG:
            print('stdout')
            print(stdout)
            print('stderr')
            print(stderr)

        if p.returncode:
            raise Exception('bpRNA failed: on %s\n%s' % (sequence, stderr))
        with open(os.path.basename(fname).replace('.dbn','.st')) as f:
            self.bprna_string = f.readlines()[5].strip()

        self.bprna_string = self.bprna_string.replace('X','E')

        if DEBUG: print(self.bprna_string)

        self.encoding_ = encode_input(self.sequence, self.bprna_string)

        if DEBUG: print("encoding shape", self.encoding_.shape)

        self.degscore_by_position = np.sum(self.encoding_ * self.coeffs_, axis=1) + self.intercept_
        self.degscore = np.sum(self.degscore_by_position)

        os.remove(fname)
        os.remove(os.path.basename(fname).replace('dbn','st'))

        # counter = Counter(bprna_string)
        # self.bprna_string = bprna_string
        # self.counts = {}
        # for k in self.weights.keys():
        #     if k in counter.keys():
        #         self.counts[k] = counter[k]
        #     else:
        #         self.counts[k] = 0
        # self.score = np.sum([self.counts[k]*self.weights[k] for k in self.weights.keys()])


    # def set_weights(self, weights_dict):
    #     if not isinstance(weights_dict, dict):
    #         raise RuntimeError('Provided weights_dict must be a dictionary instance.')
    #     for k in self.weights.keys():
    #         if k not in weights_dict.keys():
    #             raise RuntimeError("%s not in weights_dict keys." % k)
    #     self.weights = weights_dict
    #     self.score = np.sum([self.counts[k]*self.weights[k] for k in self.weights.keys()])

