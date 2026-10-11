import os
from Bio import SeqIO


# Input folders and output files for each terminator analysis
analysis_paths = {
    'full': {
        'folders': ['intermediate_data/E_coli/intergenic_regions/ortologs_divergons_align/terminators_alignments/both_flanks/',
                    'intermediate_data/E_coli/intergenic_regions/ortologs_sinvergons+_align/terminators_alignments/both_flanks/',
                    'intermediate_data/E_coli/intergenic_regions/ortologs_sinvergons-_align/terminators_alignments/both_flanks/',
                    'intermediate_data/E_coli/intergenic_regions/ortologs_convergons_align/terminators_alignments/both_flanks/'],
        'output': 'results/E_coli/term.tsv'
    },
    'right': {
        'folders': ['intermediate_data/E_coli/intergenic_regions/ortologs_divergons_align/terminators_alignments/right_flank/',
                    'intermediate_data/E_coli/intergenic_regions/ortologs_sinvergons+_align/terminators_alignments/right_flank/',
                    'intermediate_data/E_coli/intergenic_regions/ortologs_sinvergons-_align/terminators_alignments/right_flank/',
                    'intermediate_data/E_coli/intergenic_regions/ortologs_convergons_align/terminators_alignments/right_flank/'],
        'output': 'results/E_coli/term_right.tsv'
    }
}


# Count polymorphic positions
def count_pol(result, start, end):
    part_count = 0
    for i in range(start, end):
        if result[i] > 1:
            part_count += 1
    return part_count


# Read terminator annotations
term = {}
with open('data/E_coli/regulatory_elements/terminators.tsv', 'r') as inp:
    for line in inp:
        idd = line.split('\t')[0].strip()
        term[idd] = line.strip()


for analysis_type, analysis in analysis_paths.items():
    folder_paths = analysis['folders']
    output_file = analysis['output']

    with open(output_file, 'w') as out:
        if analysis_type == 'full':
            print("tmId	seq	type	name	in_site	left_site	right_site	size	ort	mezhgen	orientation", file=out)
        else:
            print("tmId	seq	type	name	in_site	right_site	size	ort	mezhgen	orientation", file=out)


    for folder_path in folder_paths:
        dir_names = os.listdir(folder_path)

        for file in dir_names:
            empty_el = 0
            mezh_gen = '_'.join(file.split('_')[0:-3])
            orientation = file.split('_')[-1]
            out_file = folder_path + '/' + file

            try:
                sequences = [str(record.seq) for record in SeqIO.parse(out_file, "fasta")]
                sequence_length = len(sequences[0])
            except:
                continue

            if not sequences:
                continue

            size = len(sequences)

            # Split the sequence into equal regions 
            if analysis_type == 'full':
                part_size = int(sequence_length / 3)
            else:
                part_size = int(sequence_length / 2)

            position_counts = {}

            # Initialize nucleotide and gap counts at each position
            for i in range(sequence_length):
                position_counts[i] = {'a': 0, 't': 0, 'g': 0, 'c': 0, '-': 0}

            # Exclude sequences with excessive gaps in the analyzed region
            for seq in sequences:
                if analysis_type == 'right':
                    count = seq.count("-")
                    if count * 2 >= sequence_length:
                        empty_el += 1
                        continue

                for i in range(sequence_length):
                    if seq[i] in position_counts[i]:
                        position_counts[i][seq[i]] += 1

            name = file.split('_')[-3].strip()

            if analysis_type == 'right':
                if empty_el > 1:
                    print(empty_el, term[name], folder_path, mezh_gen, orientation)
                if empty_el > 400:
                    continue

            # Calculate the number of sequences differing from the major allele
            result = {}
            pos = 0
            for key, counts in position_counts.items():
                most_common_letter = max(counts, key=counts.get)
                total_count = 0
                for letter, value in counts.items():
                    if letter != most_common_letter:
                        total_count += value

                result[pos] = total_count
                pos += 1

            if analysis_type == 'full':
                # Calculate polymorphism in the left flank, terminator, and right flank
                first_part_count = count_pol(result, 0, part_size)
                second_part_count = count_pol(result, part_size, 2 * part_size)
                third_part_count = count_pol(result, 2 * part_size, 3 * part_size)

                in_site = second_part_count / part_size

                if orientation == '1':
                    left_site = first_part_count / part_size
                    right_site = third_part_count / part_size
                else:
                    right_site = first_part_count / part_size
                    left_site = third_part_count / part_size

                # Save full terminator results with both flanks
                print(term[name], in_site, left_site, right_site, size, folder_path, mezh_gen, orientation, sep='\t', file=out)

            else:
                # Calculate polymorphism in the terminator and downstream flank
                first_part_count = count_pol(result, 0, part_size)
                second_part_count = count_pol(result, part_size, 2 * part_size)

                if orientation == '1':
                    in_site = first_part_count / part_size
                    right_site = second_part_count / part_size
                else:
                    right_site = first_part_count / part_size
                    in_site = second_part_count / part_size

                print(term[name], in_site, right_site, size, folder_path, mezh_gen, orientation, sep='\t', file=out)