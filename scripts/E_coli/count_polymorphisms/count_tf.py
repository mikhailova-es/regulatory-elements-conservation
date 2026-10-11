import os
from Bio import SeqIO
import numpy as np


# Input folders
folder_paths = ['intermediate_data/E_coli/intergenic_regions/ortologs_divergons_align/tf_factors_alignments/', 
                'intermediate_data/E_coli/intergenic_regions/ortologs_sinvergons+_align/tf_factors_alignments/', 
                'intermediate_data/E_coli/intergenic_regions/ortologs_sinvergons-_align/tf_factors_alignments/', 
                'intermediate_data/E_coli/intergenic_regions/ortologs_convergons_align/tf_factors_alignments/']

# Polymorphism threshold: change to 1, 5, 10, or 100
polymorphism_threshold = 100


# Count positions with polymorphism above the selected threshold
def count_pol(result, start, end):
    part_count = 0
    for i in range(start, end):
        if result[i] > polymorphism_threshold:
            part_count += 1
    return part_count


tf = {}

# Read transcription factor annotations
with open('data/E_coli/regulatory_elements/tf-factors.tsv', 'r') as inp:
    for line in inp:
        idd = line.split('\t')[0].strip()
        tf[idd] = line.strip()


output_file = f'results/E_coli/tf_ecoli_{polymorphism_threshold}.tsv'

with open(output_file, 'w') as out1:
    print("riId	regulatorName	seq	class	gene	sigma	reg_gene	confidenceLevel	in_site	left_site	right_site	size	ort	mezhgen	orientation	mean_x	mean_y	sd_low	sd_up", file=out1)

    for folder_path in folder_paths:
        dir_names = os.listdir(folder_path)

        for file in dir_names:
            empty_el = 0
            mezh_gen = '_'.join(file.split('_')[0:-3])
            orientation = file.split('_')[-1]
            out_file = folder_path + '/' + file
            name = file.split('_')[-3].strip()

            try:
                sequences = [str(record.seq) for record in SeqIO.parse(out_file, "fasta")]
                sequence_length = len(sequences[0])
            except:
                continue

            if not sequences:
                continue

            size = len(sequences)
            part_size = int(sequence_length / 3)
            position_counts = {}

            # Initialize nucleotide and gap counts at each position
            for i in range(sequence_length):
                position_counts[i] = {'a': 0, 't': 0, 'g': 0, 'c': 0, '-': 0}

            # Exclude sequences with excessive gaps in the TF binding site
            for seq in sequences:
                count = seq[part_size:2 * part_size].count("-")
                if count * 1.5 > part_size:
                    empty_el += 1
                    continue

                for i in range(sequence_length):
                    if seq[i] in position_counts[i]:
                        position_counts[i][seq[i]] += 1

            if empty_el > 1 and name in tf:
                print(empty_el, tf[name], folder_path, mezh_gen, orientation)

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

            # Count polymorphic positions in the left flank, TF binding site, and right flank
            first_part_count = count_pol(result, 0, part_size)
            second_part_count = count_pol(result, part_size, 2 * part_size)
            third_part_count = count_pol(result, 2 * part_size, 3 * part_size)
            total_pol = first_part_count + second_part_count + third_part_count

            x_random_list = []
            y_random_list = []
            deviations = []

            # Randomize polymorphic positions to define the neutral zone
            for i in range(5000):
                random_positions = np.random.choice(sequence_length, size=total_pol, replace=False)
                poly_in_left = sum(1 for pos in random_positions if pos < part_size)
                poly_in_site = sum(1 for pos in random_positions if part_size <= pos < part_size * 2)
                poly_in_right = sum(1 for pos in random_positions if pos >= part_size * 2)

                x_random = poly_in_site / part_size
                y_random = (poly_in_left + poly_in_right) / (2 * part_size)
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

            in_site = second_part_count / part_size

            if orientation == '1':
                left_site = first_part_count / part_size
                right_site = third_part_count / part_size
            else:
                right_site = first_part_count / part_size
                left_site = third_part_count / part_size
                
            # Save results
            if name in tf:
                print(tf[name], in_site, left_site, right_site, size, folder_path, mezh_gen, orientation, mean_x, mean_y, lower_sd, upper_sd, sep='\t', file=out1)