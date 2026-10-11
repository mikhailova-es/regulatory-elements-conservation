import os
from Bio import SeqIO

# Input folders and output file
folder_paths = ['intermediate_data/E_coli/intergenic_regions/ortologs_divergons_align/rnas_alignments/without_flanks/', 
                'intermediate_data/E_coli/intergenic_regions/ortologs_sinvergons+_align/rnas_alignments/without_flanks/', 
                'intermediate_data/E_coli/intergenic_regions/ortologs_sinvergons-_align/rnas_alignments/without_flanks/', 
                'intermediate_data/E_coli/intergenic_regions/ortologs_convergons_align/rnas_alignments/without_flanks/']

output_file = 'results/E_coli/rnas_structure.tsv'

# Read RNA secondary structures and identify paired and unpaired positions
rnas_dot = {}
with open('data/E_coli/regulatory_elements/rnas_dot.txt', 'r') as dot:
    line = dot.readline()
    while line:
        if line.startswith('>'):
            name_rna = line.strip()[1:]
            line = dot.readline()
            paired = []
            unpaired = []

            for i in range(len(line.strip())):
                if line[i] == ',' or line[i] == '.' or line[i] == '_' or line[i] == '-':
                    unpaired.append(i)
                else:
                    paired.append(i)
            rnas_dot[name_rna] = [unpaired, paired]
        line = dot.readline()


# Read RNA annotations
rnas = {}
with open('data/E_coli/regulatory_elements/rnas.tsv', 'r') as inp:
    for line in inp:
        rdb = line.split('\t')[0].strip()
        rnas[rdb] = line.strip()


# Count polymorphic positions separately for paired and unpaired regions
def count_pol(result, start, end, name, rnas_dot):
    unpaired = 0
    paired = 0

    for i in range(start, end):
        if result[i] > 1 and i in rnas_dot[name][0]:
            unpaired += 1
        if result[i] > 1 and i in rnas_dot[name][1]:
            paired += 1

    return unpaired, paired

with open(output_file, 'w') as out:
    print("riID	regulatorName	class	seq	paired	unpaired	total	size	ort	mezhgen	orientation", file = out)
    for folder_path in folder_paths:
        path = folder_path
        dir_names = [f for f in os.listdir(path)]

        for file in dir_names:
            mezh_gen = '_'.join(file.split('_')[0:-3])
            orientation = file.split('_')[-1]
            out_file = path + '/' + file
            sequence_length = int(file.split('_')[-2].strip())
            sequences = []

            # Read sequences of the expected length and orient them consistently
            try:
                for record in SeqIO.parse(out_file, "fasta"):
                    if len(record.seq) == sequence_length:
                        if orientation == '-1':
                            sequences.append(str(record.seq.reverse_complement()))
                        else:
                            sequences.append(str(record.seq))
            except:
                continue

            if not sequences:
                continue

            name = file.split('_')[-3].strip()
            name_leng = name + '_' + str(sequence_length)
            size = len(sequences)

            # Count nucleotides and gaps at each alignment position
            position_counts = {}
            for i in range(sequence_length):
                position_counts[i] = {'a': 0, 't': 0, 'g': 0, 'c': 0, '-': 0}

            for seq in sequences:
                for i in range(sequence_length):
                    if seq[i] in position_counts[i]:
                        position_counts[i][seq[i]] += 1

            # Calculate the number of sequences differing from the major allele
            result = {}
            for key, counts in position_counts.items():
                most_common_letter = max(counts, key=counts.get)
                total_count = sum(value for letter, value in counts.items() if letter != most_common_letter)
                result[key] = total_count

            # Count polymorphic paired and unpaired positions
            paired_count = count_pol(result, 0, sequence_length, name_leng, rnas_dot)[1]
            unpaired_count = count_pol(result, 0, sequence_length, name_leng, rnas_dot)[0]

            paired_col = len(rnas_dot[name_leng][1])
            unpaired_col = len(rnas_dot[name_leng][0])

            # Calculate polymorphism rates for paired and unpaired positions
            paired = paired_count / paired_col
            unpaired = unpaired_count / unpaired_col
            total = (paired_count + unpaired_count) / sequence_length

            # Save structure-specific polymorphism rates
            if len(name.split('-')) == 2:
                name = name.split('-')[0]
            print(rnas[name], paired, unpaired, total, size, folder_path, mezh_gen, orientation, sep='\t', file=out)