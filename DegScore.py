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
    def __init__(self, sequence, package='eternafold', linear=True):
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
        self.structure = mfe(sequence, package=package, linear=linear)
        assert len(sequence) == len(self.structure)

        fname = "%s.dbn" % filename()

        # old weights from DegScore 1
        #self.weights = {'H':0.7, 'E':1.0, 'S':0.0, 'I':0.2, 'B': 0.8, 'M': 0.6}
        self.coeffs_ = np.loadtxt('DegScore2.1_coeffs.txt',usecols=1)
        self.intercept_ = 1.1220380163671115

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
        print(self.bprna_string)

        self.encoding_ = encode_input(self.sequence, self.bprna_string)

        print(self.encoding_.shape)

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

