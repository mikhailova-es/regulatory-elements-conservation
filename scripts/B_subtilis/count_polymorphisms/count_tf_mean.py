import os
from Bio import SeqIO
import numpy as np

# Input folders and output file
folder_paths = ['intermediate_data/B_subtilis/intergenic_regions/ortologs_divergons_align/tf_factors_alignments/', 
                'intermediate_data/B_subtilis/intergenic_regions/ortologs_sinvergons+_align/tf_factors_alignments/', 
                'intermediate_data/B_subtilis/intergenic_regions/ortologs_sinvergons-_align/tf_factors_alignments/', 
                'intermediate_data/B_subtilis/intergenic_regions/ortologs_convergons_align/tf_factors_alignments/']

output_file = 'results/B_subtilis/tf_mean.tsv'

# Count polymorphic positions
def count_pol(result, start, end):
    part_count = 0
    for i in range(start, end):
        if result[i] > 1:
            part_count += 1
    return part_count

# Read transcription factor annotations
tf_name = {}
tf_stat = {}

with open('data/B_subtilis/regulatory_elements/tf-factors.tsv', 'r') as inp:
    for line in inp:
        idd = line.split('\t')[0].strip()
        tf_name[idd] = line.strip().split('\t')[1]

with open(output_file, 'w') as out1:
    print('name	mean_in_site	mean_out_site	mean_x	mean_y	sd_low	sd_up', file=out1)

    for folder_path in folder_paths:
        dir_names = os.listdir(folder_path)

        for file in dir_names:
            id = file.split('_')[-5].strip()
            name = tf_name[id]
            mezh_gen = '_'.join(file.split('_')[0:-5])
            orientation = file.split('_')[-1]
            out_file = folder_path + '/' + file
            print(file)
            tf_start = int(file.split('_')[-3])
            tf_end = int(file.split('_')[-2])
            true_length = (tf_end - tf_start) * 3

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

            # Calculate number of sequences differing from the major allele
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

            # Split the sequence into the site and two flanking regions

            part_size = int(sequence_length / 3)

            first_part_count = count_pol(result, 0, part_size)
            second_part_count = count_pol(result, part_size, 2 * part_size)
            third_part_count = count_pol(result, 2 * part_size, 3 * part_size)
            total_pol = first_part_count + second_part_count + third_part_count

            in_site = second_part_count / part_size
            if orientation == '1':
                left_site = first_part_count / part_size
                right_site = third_part_count / part_size
            else:
                right_site = first_part_count / part_size
                left_site = third_part_count / part_size
            out_site = (left_site + right_site) / 2

            if name not in tf_stat:
                tf_stat[name] = [[in_site], [out_site], [part_size], [total_pol], [sequence_length]]
            else:
                tf_stat[name][0].append(in_site)
                tf_stat[name][1].append(out_site)
                tf_stat[name][2].append(part_size)
                tf_stat[name][3].append(total_pol)
                tf_stat[name][4].append(sequence_length)

    for key, value in tf_stat.items():
        sequence_lengths = value[4]
        total_pols = value[3]
        part_sizes = value[2]
        out_sites = value[1]
        in_sites = value[0]

        # Randomize polymorphic positions to define the neutral zone
        x_random_list = []
        y_random_list = []
        deviations = []
        boundaries = np.cumsum([0] + sequence_lengths)

        print(total_pols, sequence_lengths)

        for i in range(5000):
            random_positions = np.random.choice(sum(sequence_lengths), size=sum(total_pols), replace=False)
            x_random_this_round = []
            y_random_this_round = []

            for j in range(len(sequence_lengths)):
                start = boundaries[j]
                end = boundaries[j + 1]
                part_size = part_sizes[j]

                positions_in_seq = [pos - start for pos in random_positions if start <= pos < end]

                poly_in_left = sum(pos < part_size for pos in positions_in_seq)
                poly_in_site = sum(part_size <= pos < 2 * part_size for pos in positions_in_seq)
                poly_in_right = sum(pos >= 2 * part_size for pos in positions_in_seq)

                x_random = poly_in_site / part_size
                y_random = (poly_in_left + poly_in_right) / (2 * part_size)

                x_random_this_round.append(x_random)
                y_random_this_round.append(y_random)

            x_factor_random = np.mean(x_random_this_round)
            y_factor_random = np.mean(y_random_this_round)
            x_random_list.append(x_factor_random)
            y_random_list.append(y_factor_random)

            # Calculate the distance from the y = x line
            d_factor_random = (y_factor_random - x_factor_random) / np.sqrt(2)
            deviations.append(d_factor_random)

        x_random_list = np.array(x_random_list)
        y_random_list = np.array(y_random_list)
        deviations = np.array(deviations)

        # Calculate the mean randomized point and deviation range
        mean_x = np.mean(x_random_list)
        mean_y = np.mean(y_random_list)
        mean_d = np.mean(deviations)
        sd_d = np.std(deviations)
        lower_sd = mean_d - sd_d
        upper_sd = mean_d + sd_d

        # Observed factor-level point
        mean_in_site = np.mean(in_sites)
        mean_out_site = np.mean(out_sites)

        print(key, mean_in_site, mean_out_site, mean_x, mean_y, lower_sd, upper_sd, sep='\t', file=out1)