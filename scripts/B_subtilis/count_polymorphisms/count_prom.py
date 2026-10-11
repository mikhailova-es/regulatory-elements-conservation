import os
from Bio import SeqIO
import numpy as np

folder_paths = ['intermediate_data/B_subtilis/intergenic_regions/ortologs_divergons_align/promoters_alignments/', 
                'intermediate_data/B_subtilis/intergenic_regions/ortologs_sinvergons+_align/promoters_alignments/', 
                'intermediate_data/B_subtilis/intergenic_regions/ortologs_sinvergons-_align/promoters_alignments/', 
                'intermediate_data/B_subtilis/intergenic_regions/ortologs_convergons_align/promoters_alignments/']

prox_out = 'results/B_subtilis/prox.tsv'
dist_out = 'results/B_subtilis/dist.tsv'

# Count polymorphic positions 
def count_pol(result, start, end):
    part_count = 0
    for i in range(start, end):
        if result[i] > 1:
            part_count += 1
    return part_count


prom = {}

# Read promoter annotations
with open('data/B_subtilis/regulatory_elements/promoters.tsv', 'r') as inp:
    for line in inp:
        idd = line.split('\t')[0].strip()
        prom[idd] = line.strip()


with open(prox_out, 'w') as out_prox, open(dist_out, 'w') as out_dist:
    print("pmId	sigmaFactor	pmName	gene	coordinates	position	pmSequence	confidenceLevel	site_start_prox	site_end_prox	site_start_dist	site_end_dist	in_site	left_site	right_site	size	ort	mezhgen	orientation	mean_x	mean_y	sd_low	sd_up", file = out_prox)
    print("pmId	sigmaFactor	pmName	gene	coordinates	position	pmSequence	confidenceLevel	site_start_prox	site_end_prox	site_start_dist	site_end_dist	in_site	left_site	right_site	size	ort	mezhgen	orientation	mean_x	mean_y	sd_low	sd_up", file=out_dist)
    for folder_path in folder_paths:
        dir_names = os.listdir(folder_path)

        for file in dir_names:
            type = file.split('_')[0]
            mezh_gen = '_'.join(file.split('_')[1:-7])
            orientation = file.split('_')[-1]
            out_file = folder_path + '/' + file

            # Extract promoter coordinates from the file name
            try:
                sequences = [str(record.seq) for record in SeqIO.parse(out_file, "fasta")]
                sequence_length = len(sequences[0])
            except:
                continue

            if not sequences:
                continue

            size = len(sequences)

            # Count nucleotide and gap at every alignment position
            position_counts = {}
            for i in range(sequence_length):
                position_counts[i] = {'a': 0, 't': 0, 'g': 0, 'c': 0, '-': 0}

            for seq in sequences:
                for i in range(sequence_length):
                    if seq[i] in position_counts[i]:
                        position_counts[i][seq[i]] += 1

            # Count sequences that differ from the major allele at each position
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

            # Divide the sequence into left flank, promoter, and right flank
            part_size = int(sequence_length / 3)

            first_part_count = count_pol(result, 0, part_size)
            second_part_count = count_pol(result, part_size, 2 * part_size)
            third_part_count = count_pol(result, 2 * part_size, 3 * part_size)
            total_pol = first_part_count + second_part_count + third_part_count

            # Generate the expected distribution of polymorphisms by randomization
            x_random_list = []
            y_random_list = []
            deviations = []

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

            # Assign the two flanks according to the strand orientation
            if orientation == '1':
                left_site = first_part_count / part_size
                right_site = third_part_count / part_size
            else:
                right_site = first_part_count / part_size
                left_site = third_part_count / part_size

            name = file.split('_')[-7].strip()

            # Save results for proximal and distal promoter elements
            if type == 'prox':
                print(prom[name], in_site, left_site, right_site, size, folder_path, mezh_gen, orientation, mean_x, mean_y, lower_sd, upper_sd, sep='\t', file=out_prox)

            if type == 'dist':
                print(prom[name], in_site, left_site, right_site, size, folder_path, mezh_gen, orientation, mean_x, mean_y, lower_sd, upper_sd, sep='\t', file=out_dist)