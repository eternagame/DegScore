
## Setup:

1. Clone `bpRNA` from https://github.com/hendrixlab/bpRNA.

2. Add to arnie file:

  `bprna: /path/to/bpRNA`

3. May need to install perl package `Graph`: if not installed, will get an error like `Can't locate Graph.pm in @INC`.

  Sherlock: Run `ml perl; cpanm Graph`.
  Otherwise, something like `sudo cpan install Graph`.


## Usage:
Update Nov 10:

No longer using LinearFold-E, uses just plain EternaFold now.

\[Update for Nov 2020 OpenVaccine work\]:

The DegScore class uses Arnie to calculate an MFE structure, which is then used with bpRNA to compute the loop type string (`bprna_string`).

The DEFAULT setting in this class is LinearFold-E (LinearFold with EternaFold) to mimic DegScore 2.1 calculation. This requires having successfully set up arnie with LinearFold-E.

Example syntax:

```
from DegScore import DegScore

sequence='GGGGAAACCCC'

mdl = DegScore(sequence)
print('bpRNA string:')
print(mdl.bprna_string)
print('Predicted degscore per nucleotide:')
print(mdl.degscore_by_position)
print('Total degscore:')
print(mdl.degscore)
```

Output:
```
bpRNA string:
SSSSHHHSSSS
Predicted degscore per nucleotide:
[0.173 0.129 0.277 0.388 0.414 0.434 0.445 0.233 0.104 0.122 0.702]
Total degscore:
3.421000000000001
```

To change from these default LinearFold-E settings for MFE calculation to other Arnie settings:

ViennaRNA:
```
mdl = DegScore(sequence, package='vienna', linear=False)
```

LinearFold-V:
```
mdl = DegScore(sequence, package='vienna', linear=True)

```
Eternafold:
```
mdl = DegScore(sequence, package='eternafold', linear=False)
```
