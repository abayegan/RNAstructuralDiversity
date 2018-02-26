Amir Bayegan, Peter Clote -- Boston College


A number of measures of structural diversity have been introduced in the literature, ranging from positional entropy to ensemble defect. Here, we unified  such measures into a C++ program. This program can efficiently compute 17 measures for RNA structural analysis including some newly defined ones all of which depend on base pairing probabilities.  


USAGE: ./strDiversity -s sequence -f fasta  [-d 0|2] [-c tragetStructure] [-t temperature] [-q tragetSequence] [-u turner99|turner04|andronescu07] [-x epsilon] [-o output] [-n] [-v]
try ./strDiversity -h for help


compute various structural diversity measures:

	-s sequence: The input RNA sequence

	-f fasta: The input fasta file. Target structure for each rna can be define after its sequence. Otherwise MFE structure will be used. 

	-d 0|2: set dangling end contributions to 0 or 2 (default is 2)

	-c targetStructure: target structure for computation of expected base pair distance,ensemble defect,expected number of native contacts
	and expected ensemble distance. if it is not defined the minimum free energy structure is used to calculate these measures

	-q targetSequence: target sequence for computation of ensemble diversity between two sequences

	-t temperature: fold sequences at the given temperature

	-u turner99|turner04|andronescu07 : energy parameter used to fold the sequences

	-n : normalize measures by length

	-x : epsilon for calculation of epsilon-near expected bp distance and ensemble defect (default is 2)

	-v : verbose mode

	-h : print help

	-o : measures to be computed separated by comma. Each measure is defined by an integer value:
	1: expected positional entropy(if -n is not used computes the sum over all poistions)
	2:Morgan Higgs diversity
	3: Vienna diversity
	4: expected base pair distance
	5: ensemble defect
	6: expected number of base pairs
	7: expected proportion of target contacts
	8: expected 2D contact order(if -n is not used computes the sum over all poistions)
	9: epsilon-near BP distance
	10: epsilon-near ensemble defect
	11: robustness
	12: ensemble Hamming distance
	13: ensemble base pair distance
	14: positional entropy distance
	15: expected height
	16: expected modified height
	17: For each position the probability of bein basepaired to left,right and being unpaired is reported respectively.	Example: using -o 1,2,5,7 will only compute the corresppnding measures. 
