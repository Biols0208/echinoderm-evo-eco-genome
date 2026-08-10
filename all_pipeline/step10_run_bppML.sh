#!/bin/bash
bpp2.sif bppml param=bppML_V3.bpp
bpp2.sif mapnh param=mapNH_V3.bpp

<<!!bppML_V3.bpp
DATA = rename.keep.sp.id.aln

# Sequences:

# The alphabet to use:
# DNA, RNA or Protein
alphabet = Codon(letter=DNA)
genetic_code = Standard

# The sequence file to use (sequences must be aligned!)
input.data1=alignment(file = ./rename.keep.sp.id.aln, format=Fasta, sites_to_use = all, max_gap_allowed=50%, remove_stop_codons = yes,remove_saturated_sites = yes)

# ----------------------------------------------------------------------------------------
#                                     Input tree file
# ----------------------------------------------------------------------------------------

# user or random
input.tree1=user(file = ./rename.tree.pruned)

# ----------------------------------------------------------------------------------------
#                                     Process specification
# ----------------------------------------------------------------------------------------
# See the manual for a description of the syntax and available options.
#

model1 = YN98(kappa=1, omega=1, frequencies=F3X4(), initFreqs=observed, data=1)

rate_distribution1 = Constant()

process1 = Homogeneous(model=1, tree=1, rate=1)

phylo1 = Single(process=1, data=1)

# ----------------------------------------------------------------------------------------
#                                     Optimization
# ----------------------------------------------------------------------------------------

# Should we reestimate likelihood parameters? Tree topology will not be optimized.
# (recommanded)
# Method to use for optimizing numerical parameters:
# - None, no optimization performed
# - DB derivatives for branch lengths + Brent for other parameters
# - FullD derivatives for all parameters, using numerical derivatives for non-branch lengths parameters.

optimization = FullD(derivatives=Newton)

# Parameters to ignore (for instance equilibrium frequencies)
optimization.ignore_parameters =

# Maximum number of likelihood evaluations:
optimization.max_number_f_eval = 10000

# Precision to reach:
optimization.tolerance = 0.000001

# idem for error or warning messages:
optimization.message_handler = $(DATA).ml_h.messages

# A file where to dump optimization steps (a file path or std for standard output)
optimization.profiler = $(DATA).ml_h.profile


# Should we write the resulting tree? none or file name.
output.tree.file = $(DATA).ml_h.dnd
output.tree.format = Newick

# Alignment information log file (site specific rates, etc):
output.infos = $(DATA).ml_h.infos

# Write numerical parameter estimated values:
output.estimates = $(DATA).ml_h.params.bpp
!!bppML_V3.bpp

<<!!mapNH_V3.bpp
DATA = rename.keep.sp.id.aln
TYPE = DnDs

alphabet = Codon(letter=DNA)
genetic_code = Standard

input.data1=alignment(file = ./$(DATA), format=Fasta, sites_to_use = all, max_gap_allowed=50%, remove_stop_codons = yes,remove_saturated_sites = yes)


#Set imput tree file. Original branch lengths will be used for mapping.
input.tree1=user(file = $(DATA).ml_h.dnd_1)

### File of the modeling

param = $(DATA).ml_h.params.bpp

### normalization with the same model in which omega=1

nullProcessParams = YN98.omega*=1

### Type of event counted

map.type = $(TYPE)

count.max = 3

# Split counts from model opportunities
output.counts = PerBranch(prefix=CM_splitcounts_, splitNorm=true, perWordSize=true)
!!mapNH_V3.bpp
