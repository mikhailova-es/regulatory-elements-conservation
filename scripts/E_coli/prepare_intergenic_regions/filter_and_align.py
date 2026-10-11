import os
from statistics import median
from Bio import SeqIO
import subprocess
from Bio.Seq import Seq

# Input folders
sin_minus_dir = 'intermediate_data/E_coli/intergenic_regions/ortologs_sinvergons-/'
sin_plus_dir = 'intermediate_data/E_coli/intergenic_regions/ortologs_sinvergons+/'
divergent_dir = 'intermediate_data/E_coli/intergenic_regions/ortologs_divergons/'
convergent_dir = 'intermediate_data/E_coli/intergenic_regions/ortologs_convergons/'
folder_paths = [sin_minus_dir, sin_plus_dir, divergent_dir, convergent_dir]
filter_folder_paths = [path + '_filter/' for path in folder_paths]

mafft_path = '/path/to/mafft'

# Filter sequences by length and ambiguous nucleotides
for folder in folder_paths:
    data = os.listdir(folder)

    for filename in data:
        input_file = folder + '/' + filename

        if not filename.endswith('.fna'):
            continue

        records = {}

        for record in SeqIO.parse(input_file, 'fasta'):
            if record.id not in records:
                records[record.id] = record.seq

        if len(records) != 0:
            lengths = [len(seq) for seq in records.values()]
            median_length = median(lengths)
        else:
            median_length = 0

        to_delete = []

        for seq_id, seq in records.items():
            seq_len = len(seq)

            if seq_len == 0:
                to_delete.append(seq_id)
                continue

            n_count = str(seq).count('N')

            if seq_len >= median_length * 2 or seq_len <= median_length / 2 or n_count / seq_len > 0.20:
                to_delete.append(seq_id)

        for seq_id in to_delete:
            del records[seq_id]

        if len(records) != 0:
            lengths = [len(seq) for seq in records.values()]
            median_length = median(lengths)
        else:
            median_length = 0

        if median_length > 30 and len(records) > 2000:
            with open(folder + '_filter/' + filename, 'w') as out:
                for seq_id, seq in records.items():
                    print(f">{seq_id}", file=out)
                    print(seq, file=out)


# Align filtered sequences with MAFFT
for folder_path in filter_folder_paths:
    files = os.listdir(folder_path)

    for file in files:
        if not file.endswith('.fna'):
            continue

        inp_file = folder_path + '/' + file
        out_file = folder_path[:-7] + '_align/' + file

        # Run MAFFT alignment
        cmd = mafft_path + ' --maxiterate 2 --quiet --thread 5 ' + inp_file + ' > ' + out_file
        subprocess.run(cmd, shell=True)

        align_records = list(SeqIO.parse(out_file, 'fasta'))
        align_sequences = [str(align_record.seq) for align_record in align_records]
        sequence_length = len(align_sequences[0])

        position_counts = {}

        for i in range(sequence_length):
            position_counts[i] = {'a': 0, 't': 0, 'g': 0, 'c': 0, '-': 0}

        for al_seq in align_sequences:
            for i in range(sequence_length):
                if al_seq[i] in position_counts[i]:
                    position_counts[i][al_seq[i]] += 1

        keep_position = []

        for pos, counts in position_counts.items():
            most_common_letter = max(counts, key=counts.get)

            if most_common_letter == "-":
                letters_count = sum(value for letter, value in counts.items() if letter != "-")
                if letters_count <= 3:
                    continue

            keep_position.append(pos)

        # Write the filtered alignment
        new_records = []

        for recordd, seqq in zip(align_records, align_sequences):
            new_seq = "".join(seqq[pos] for pos in keep_position)
            recordd.seq = Seq(new_seq)
            new_records.append(recordd)

        SeqIO.write(new_records, out_file, "fasta")