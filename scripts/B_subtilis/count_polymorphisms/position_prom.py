import os
from Bio import SeqIO

# Input folders
folder_paths = ['intermediate_data/B_subtilis/intergenic_regions/ortologs_divergons_align/promoters_alignments/', 
                'intermediate_data/B_subtilis/intergenic_regions/ortologs_sinvergons+_align/promoters_alignments/', 
                'intermediate_data/B_subtilis/intergenic_regions/ortologs_sinvergons-_align/promoters_alignments/', 
                'intermediate_data/B_subtilis/intergenic_regions/ortologs_convergons_align/promoters_alignments/']

# Read promoter annotations
promoter = {}
with open('data/B_subtilis/regulatory_elements/promoters.tsv', 'r') as inp:
    for line in inp:
        idd = line.split('\t')[0].strip()
        promoter[idd] = line.strip()

prox_count = 0
dist_count = 0
total_count_d_prox = {i: 0 for i in range(18)}
total_count_d_dist = {i: 0 for i in range(18)}

for folder_path in folder_paths:
    dir_names = os.listdir(folder_path)

    for file in dir_names:
        name = file.split('_')[-7].strip()
        sigma = promoter[name].split('\t')[1]
        type = file.split('_')[0]
        mezh_gen = '_'.join(file.split('_')[1:-7])
        orientation = file.split('_')[-1]
        out_file = folder_path + '/' + file

# Use only SigA promoters
        if sigma != 'SigA':
            continue

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

        # Count promoters with polymorphic positions in proximal element
        if type == 'prox' and orientation == '1' and len(result) == 18:
            prox_count += 1
            for i in range(18):
                if result[i] > 1:
                    total_count_d_prox[i] += 1

        if type == 'prox' and orientation == '-1' and len(result) == 18:
            prox_count += 1
            for i in range(17, -1, -1):
                if result[i] > 1:
                    total_count_d_prox[17 - i] += 1

        # Count promoters with polymorphic positions in distal element
        if type == 'dist' and orientation == '1' and len(result) == 18:
            dist_count += 1
            for i in range(18):
                if result[i] > 1:
                    total_count_d_dist[i] += 1

        if type == 'dist' and orientation == '-1' and len(result) == 18:
            dist_count += 1
            for i in range(17, -1, -1):
                if result[i] > 1:
                    total_count_d_dist[17 - i] += 1

        print(prox_count)
        print(dist_count)
        print('prox', total_count_d_prox)
        print('dist', total_count_d_dist)