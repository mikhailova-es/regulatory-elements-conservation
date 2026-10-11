import os
from Bio import SeqIO


# Input folders and output file
folder_paths = ['intermediate_data/B_subtilis/intergenic_regions/ortologs_divergons_align/terminators_alignments/without_flanks/', 
                'intermediate_data/B_subtilis/intergenic_regions/ortologs_sinvergons+_align/terminators_alignments/without_flanks/', 
                'intermediate_data/B_subtilis/intergenic_regions/ortologs_sinvergons-_align/terminators_alignments/without_flanks/', 
                'intermediate_data/B_subtilis/intergenic_regions/ortologs_convergons_align/terminators_alignments/without_flanks/']

output_file = 'results/B_subtilis/term_structure.tsv'


# Calculate the fraction of polymorphic positions
def count_pol(part_size, result, start):
    end = start + part_size
    part_count = 0
    for i in range(start, end):
        if result[i] > 1:
            part_count += 1
    return part_count / part_size


# Read terminator annotations
term = {}
with open('data/B_subtilis/regulatory_elements/terminators.tsv', 'r') as inp:
    for line in inp:
        idd = line.split('\t')[0].strip()
        term[idd] = line.strip()


# Read coordinates of terminator structural elements
term_cord = {}
with open('data/B_subtilis/regulatory_elements/terminator_cord.txt', 'r') as inp2:
    for linee in inp2:
        linee = linee.strip().split('\t')
        idi = linee[0]
        h1s, h1e, h2s, h2e = linee[1], linee[2], linee[3], linee[4]
        term_cord[idi] = [h1s, h1e, h2s, h2e]

with open(output_file, 'w') as out:
    print("tmId	Genes	Operonname	Terminatorposition	DeltaG[kcal/mole]	Terminator	Evidence	References	in_tail	paired_hairpin	loop	size	ort	mezhgen	orientation", file = out)
    for folder_path in folder_paths:
        dir_names = [f for f in os.listdir(folder_path) if not f.endswith('align')]

        for file in dir_names:
            empty_el = 0
            mezh_gen = '_'.join(file.split('_')[1:-3])
            orientation = file.split('_')[-1]
            input_file = folder_path + '/' + file
            true_length = int(file.split('_')[-2].strip())

            # Read aligned terminator sequences
            try:
                sequences = [str(record.seq) for record in SeqIO.parse(input_file, "fasta")]
                sequence_length = len(sequences[0])
            except:
                continue

            if not sequences:
                print(input_file)
                continue

            size = len(sequences)

            # Count nucleotides and gaps at each position
            position_counts = {}
            for i in range(sequence_length):
                position_counts[i] = {'a': 0, 't': 0, 'g': 0, 'c': 0, '-': 0}

            for seq in sequences:
                count = seq.count("-")
                if count * 1.5 >= sequence_length:
                    empty_el += 1
                    continue

                for i in range(sequence_length):
                    if seq[i] in position_counts[i]:
                        position_counts[i][seq[i]] += 1

            # Calculate the number of sequences differing from the major nucleotide
            result = {}
            pos = 0
            for key, counts in position_counts.items():
                most_common_letter = max(counts, key=counts.get)
                total_count = sum(value for letter, value in counts.items() if letter != most_common_letter and letter != '-')
                result[pos] = total_count
                pos += 1

            # Check that the alignment length matches the expected terminator length
            if len(result) != true_length:
                continue

            part_seq_size = int(true_length)
            name = file.split('_')[-3].strip()

            if empty_el > 1:
                print(empty_el, term[name], folder_path, mezh_gen, orientation)
            if empty_el > 40:
                continue

            # Get structural coordinates
            h1s, h1e, h2s, h2e = [int(term_cord[name][i]) for i in range(4)]
            tail_length = part_seq_size - h2e

            # Calculate polymorphism for the tail and individual structural components
            if orientation == '1':
                in_tail = count_pol(tail_length, result, h2e)
                hairpin_1 = count_pol(h1e - h1s, result, h1s)
                hairpin_2 = count_pol(h2e - h2s, result, h2s)
                loop = count_pol(h2s - h1e, result, h1e)

            if orientation == '-1':
                in_tail = count_pol(tail_length, result, 0)
                hairpin_1 = count_pol(h1e - h1s, result, part_seq_size - h1e)
                hairpin_2 = count_pol(h2e - h2s, result, tail_length)
                loop = count_pol(h2s - h1e, result, part_seq_size - h2s)

            # Save polymorphism rates for terminator structural components
                print(term[name], in_tail, (hairpin_1 + hairpin_2) / 2, loop, size, folder_path, mezh_gen, orientation, sep='\t', file=out)