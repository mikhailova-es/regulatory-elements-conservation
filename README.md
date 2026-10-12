# regulatory-elements-conservation

This repository contains scripts for analyzing the conservation of bacterial regulatory elements and intergenic regions in *Escherichia coli* and *Bacillus subtilis*. 

## Analysis workflow

All scripts should be run from the root directory.

### 1. Prepare input data

Download genome sequences in FASTA format and genome annotations in GFF format from NCBI. Place the files in `data/<species>/genomes/`, with a separate directory for each genome containing its FASTA and GFF files.

Orthologous groups must also be prepared as described in the Materials and Methods section of the associated manuscript. Place the resulting tables at the following locations:

- *Bacillus subtilis*: `data/B_subtilis/orthologous_groups/proteinortho.tsv`
- *Escherichia coli*: `data/E_coli/orthologous_groups/gene_presence_absence.csv`

### 2. Extract intergenic regions

Run `scripts/<species>/prepare_intergenic_regions/extract_igr.py`.

The script combines genome annotations with orthologous group assignments and saves the resulting annotations in `intermediate_data/<species>/ortho_annotation/`.

It then extracts intergenic regions (IGRs) located between adjacent CDSs belonging to the same pair of orthologous groups across genomes. IGRs are grouped according to the orientation of the flanking CDSs.

The resulting FASTA files are saved in `intermediate_data/<species>/intergenic_regions/`, in the following directories:

- `ortologs_sinvergons-/`
- `ortologs_sinvergons+/`
- `ortologs_convergons/`
- `ortologs_divergons/`

### 3. Filter and align sequences

Run `scripts/<species>/prepare_intergenic_regions/filter_and_align.py`.

The script filters IGRs sequences according to the criteria described in the Materials and Methods section. The filtered FASTA files are saved in directories with the `_filter` suffix.

The remaining sequences are aligned using MAFFT. The alignments are saved in corresponding directories with the `_align` suffix.

MAFFT must be installed, and its executable path must be specified in the script before running it.

### 4. Map regulatory elements

The next step is to identify regulatory elements in the IGR alignments using scripts in `scripts/<species>/blast_elements/`. BLAST must be installed before running these scripts. Regulatory element sequences are stored in `data/<species>/regulatory_elements/`.

The scripts create subdirectories in the alignment directories to store the resulting alignments of regulatory elements and their flanking regions.

- `blast_prom.py` maps promoter elements and generates alignments of proximal and distal promoter elements with their flanking regions.
- `tf_blast.py` maps transcription factor binding sites (TFBSs) and generates alignments of the elements with their flanking regions.
- `blast_rna.py` maps non-coding RNAs and produces alignments with both flanks, without flanks, and with either the left or right flank.
- `blast_term.py` maps terminators and produces alignments with both flanks, without flanks, and with the right flank.


### 5. Calculate polymorphism rates

The scripts in `scripts/<species>/polymorphism_count/` use the alignments generated in the previous step, regulatory element annotations, and dot-bracket representations of RNA and terminator secondary structures stored in `data/<species>/regulatory_elements/`.

The scripts perform the following analyses:

- `count_tf.py` calculates polymorphism rates in transcription factor binding sites (TFBSs) and their flanking regions. 
- `count_tf_mean.py` calculates mean polymorphism rates across all binding sites and their flanking regions associated with each transcription factor
- `count_prom.py` calculates polymorphism rates in proximal and distal promoter elements and their flanking regions.
- `position_prom.py` analyzes polymorphism rates at individual positions of -10 and -35 promoter elements.
- `count_rna.py` calculates polymorphism rates in RNA elements and their flanking regions.
- `count_rna_structure.py` compares polymorphism rates at paired and unpaired positions in RNA secondary structures.
- `count_term.py` calculates polymorphism rates in terminator elements and their flanking regions.
- `count_term_structure.py` analyzes polymorphism rates across different regions of terminator secondary structures.
- `compens_rna.py` analyzes nucleotide substitutions at complementary base pairs to investigate compensatory substitutions in RNA secondary structures.
- `compens_term.py` analyzes nucleotide substitutions at complementary base pairs to investigate compensatory substitutions in terminator secondary structures.
Scripts that analyze elements and flanks also include a randomization procedure described in the Materials and Methods section.
The resulting tables are saved in `results/<species>/`.

The terminator structure analysis (`count_term_structure.py`) and compensatory substitution analysis in terminators (`compens_term.py`) are performed only for *Bacillus subtilis*.
### 6. Statistical analysis and visualization

The R Markdown file `figures/figure.Rmd` contains the code for statistical tests and figure generation. The resulting figures are also stored in the `figures/` directory.
