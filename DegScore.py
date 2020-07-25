from arnie.utils import *
import numpy as np
from collections import Counter
package_locs = load_package_locations()
DEBUG=False
'''
Usage with arnie:
add to arnie file:
bprna: /path/to/bpRNA
bprna is cloned from https://github.com/hendrixlab/bpRNA.
'''
class DegScore():
    def __init__(self, sequence, structure):
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
        assert len(sequence) == len(structure)
        fname = "%s.dbn" % filename()
        self.weights = {'H':0.7, 'E':1.0, 'S':0.0, 'I':0.2, 'B': 0.8, 'M': 0.6}
        _ = write([sequence, structure],fname=fname)
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
        with open(os.path.basename(fname).replace('dbn','st')) as f:
            bprna_string = f.readlines()[5].strip()
        counter = Counter(bprna_string)
        self.bprna_string = bprna_string
        self.counts = {}
        for k in self.weights.keys():
            if k in counter.keys():
                self.counts[k] = counter[k]
            else:
                self.counts[k] = 0
        self.score = np.sum([self.counts[k]*self.weights[k] for k in self.weights.keys()])
        os.remove(fname)
        os.remove(os.path.basename(fname).replace('dbn','st'))
    def set_weights(self, weights_dict):
        if not isinstance(weights_dict, dict):
            raise RuntimeError('Provided weights_dict must be a dictionary instance.')
        for k in self.weights.keys():
            if k not in weights_dict.keys():
                raise RuntimeError("%s not in weights_dict keys." % k)
        self.weights = weights_dict
        self.score = np.sum([self.counts[k]*self.weights[k] for k in self.weights.keys()])