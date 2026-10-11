import os
from Bio import SeqIO
import numpy as np


# Input folders and output files for each RNA region analysis
analysis_paths = {
    'full': {
        'folders': ['intermediate_data/B_subtilis/intergenic_regions/ortologs_divergons_align/rnas_alignments/both_flanks/', 
                'intermediate_data/B_subtilis/intergenic_regions/ortologs_sinvergons+_align/rnas_alignments/both_flanks/', 
                'intermediate_data/B_subtilis/intergenic_regions/ortologs_sinvergons-_align/rnas_alignments/both_flanks/', 
                'intermediate_data/B_subtilis/intergenic_regions/ortologs_convergons_align/rnas_alignments/both_flanks/'],
        'output': 'results/B_subtilis/rna.tsv'
    },
    'left': {
        'folders': ['intermediate_data/B_subtilis/intergenic_regions/ortologs_divergons_align/rnas_alignments/left_flank/', 
                'intermediate_data/B_subtilis/intergenic_regions/ortologs_sinvergons+_align/rnas_alignments/left_flank/', 
                'intermediate_data/B_subtilis/intergenic_regions/ortologs_sinvergons-_align/rnas_alignments/left_flank/', 
                'intermediate_data/B_subtilis/intergenic_regions/ortologs_convergons_align/rnas_alignments/left_flank/'],
        'output': 'results/B_subtilis/rna_left.tsv'
    },
    'right': {
        'folders':['intermediate_data/B_subtilis/intergenic_regions/ortologs_divergons_align/rnas_alignments/right_flank/', 
                'intermediate_data/B_subtilis/intergenic_regions/ortologs_sinvergons+_align/rnas_alignments/right_flank/', 
                'intermediate_data/B_subtilis/intergenic_regions/ortologs_sinvergons-_align/rnas_alignments/right_flank/', 
                'intermediate_data/B_subtilis/intergenic_regions/ortologs_convergons_align/rnas_alignments/right_flank/'],
        'output': 'results/B_subtilis/rna_right.tsv'
    }
}


# Count polymorphic positions
def count_pol(result, start, end):
    part_count = 0
    for i in range(start, end):
        if result[i] > 1:
            part_count += 1
    return part_count


# Read RNA annotations
rna = {}
with open('data/B_subtilis/regulatory_elements/rnas.tsv', 'r') as inp:
    for line in inp:
        idd = line.split('\t')[0].strip()
        rna[idd] = line.strip()


for analysis_type, analysis in analysis_paths.items():
    folder_paths = analysis['folders']
    output_file = analysis['output']

    with open(output_file, 'w') as out:
        if analysis_type == 'full':
            print("riId	regulatorName	class	confidenceLevel	in_site	left_site	right_site	size	ort	mezhgen	orientation	mean_x	mean_y	sd_low	sd_up", file=out)
        else:
            print("riId	regulatorName	class	confidenceLevel	in_site	 out_site	size	ort	mezhgen	orientation	mean_x	mean_y	sd_low	sd_up", file=out)

        for folder_path in folder_paths:
            dir_names = os.listdir(folder_path)

            for file in dir_names:
                mezh_gen = '_'.join(file.split('_')[0:-3])
                orientation = file.split('_')[-1]
                out_file = folder_path + '/' + file

                # Read RNA sequences
                try:
                    sequences = [str(record.seq) for record in SeqIO.parse(out_file, "fasta")]
                    sequence_length = len(sequences[0])
                except:
                    continue

                if not sequences:
                    continue

                size = len(sequences)

                # Count nucleotides and gaps at each position
                position_counts = {}
                for i in range(sequence_length):
                    position_counts[i] = {'a': 0, 't': 0, 'g': 0, 'c': 0, '-': 0}

                for seq in sequences:
                    for i in range(sequence_length):
                        if seq[i] in position_counts[i]:
                            position_counts[i][seq[i]] += 1

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

                # Split the sequence into equal regions for comparison
                if analysis_type == 'full':
                    part_size = int(sequence_length / 3)
                    first_part_count = count_pol(result, 0, part_size)
                    second_part_count = count_pol(result, part_size, 2 * part_size)
                    third_part_count = count_pol(result, 2 * part_size, 3 * part_size)
                    total_pol = first_part_count + second_part_count + third_part_count
                else:
                    part_size = int(sequence_length / 2)
                    first_part_count = count_pol(result, 0, part_size)
                    second_part_count = count_pol(result, part_size, 2 * part_size)
                    total_pol = first_part_count + second_part_count

                x_random_list = []
                y_random_list = []
                deviations = []

                # Randomize polymorphic positions to define the neutral zone
                for i in range(5000):
                    random_positions = np.random.choice(sequence_length, size=total_pol, replace=False)

                    if analysis_type == 'full':
                        poly_in_left = sum(1 for pos in random_positions if pos < part_size)
                        poly_in_site = sum(1 for pos in random_positions if part_size <= pos < part_size * 2)
                        poly_in_right = sum(1 for pos in random_positions if pos >= part_size * 2)

                        x_random = poly_in_site / part_size
                        y_random = (poly_in_left + poly_in_right) / (2 * part_size)
                    else:
                        poly_in_left = sum(1 for pos in random_positions if pos < part_size)
                        poly_in_right = sum(1 for pos in random_positions if part_size <= pos < part_size * 2)

                        x_random = poly_in_left / part_size
                        y_random = poly_in_right / part_size

                    x_random_list.append(x_random)
                    y_random_list.append(y_random)

                    # Calculate the distance from the y = x line
                    d = (y_random - x_random) / np.sqrt(2)
                    deviations.append(d)

                x_random_list = np.array(x_random_list)
                y_random_list = np.array(y_random_list)
                deviations = np.array(deviations)

                mean_x = np.mean(x_random_list)
                mean_y = np.mean(y_random_list)
                mean_d = np.mean(deviations)
                sd_d = np.std(deviations)

                lower_sd = mean_d - sd_d
                upper_sd = mean_d + sd_d

                # Calculate polymorphism rates and orient them relative to the RNA
                if analysis_type == 'full':
                    in_site = second_part_count / part_size

                    if orientation == '1':
                        left_site = first_part_count / part_size
                        right_site = third_part_count / part_size
                    else:
                        right_site = first_part_count / part_size
                        left_site = third_part_count / part_size

                elif analysis_type == 'left':
                    if orientation == '1':
                        left_site = first_part_count / part_size
                        in_site = second_part_count / part_size
                    else:
                        in_site = first_part_count / part_size
                        left_site = second_part_count / part_size

                else:
                    if orientation == '1':
                        in_site = first_part_count / part_size
                        right_site = second_part_count / part_size
                    else:
                        right_site = first_part_count / part_size
                        in_site = second_part_count / part_size

                name = file.split('_')[-3].strip()

                if len(name.split('-')) == 2:
                    name = name.split('-')[0]

                # Save results
                if analysis_type == 'full':
                    print(rna[name], name, in_site, left_site, right_site, size, folder_path, mezh_gen, orientation, mean_x, mean_y, lower_sd, upper_sd, sep='\t', file=out)
                else:
                    out_site = left_site if analysis_type == 'left' else right_site
                    print(rna[name], name, in_site, out_site, size, folder_path, mezh_gen, orientation, mean_x, mean_y, lower_sd, upper_sd, sep='\t', file=out)