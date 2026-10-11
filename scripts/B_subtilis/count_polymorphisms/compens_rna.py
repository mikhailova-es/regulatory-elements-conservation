import os
from Bio import SeqIO

rnas_dot = {}

paths = [
    'intermediate_data/B_subtilis/intergenic_regions/ortologs_divergons_align/rnas_alignments/without_flanks/',
    'intermediate_data/B_subtilis/intergenic_regions/ortologs_sinvergons+_align/rnas_alignments/without_flanks/',
    'intermediate_data/B_subtilis/intergenic_regions/ortologs_sinvergons-_align/rnas_alignments/without_flanks/',
    'intermediate_data/B_subtilis/intergenic_regions/ortologs_convergons_align/rnas_alignments/without_flanks/'
]

# Read RNA secondary structures from dot-bracket notation
with open('data/B_subtilis/regulatory_elements/rnas_dot.txt', 'r') as dot:
    line = dot.readline()
    while line:
        if line.startswith('>'):
            name_rna = line.strip().split('_')[0][1:]
            line = dot.readline()
            rnas_dot[name_rna] = line.strip()
        line = dot.readline()

# Analyze RNA base-pair variation
for path in paths:
    for file in os.listdir(path):
        stat = {}
        input_file = os.path.join(path, file)
        name = file.split('_')[-3]
        records = SeqIO.parse(input_file, "fasta")
        structure = rnas_dot[name]

        # Count observed base pairs at each paired position
        for rec in records:
            sequence = rec.seq
            stack = []
            left_paired = ['[', '<', '(']
            right_paired = [']', '>', ')']

            for i in range(len(structure)):
                try:
                    if structure[i] in left_paired:
                        stack.append(i)

                    if structure[i] in right_paired:
                        vyshel = stack.pop()

                        base_1 = sequence[vyshel].upper()
                        base_2 = sequence[i].upper()

                        if base_1 not in 'ATGC' or base_2 not in 'ATGC':
                            continue

                        pair = (vyshel, i)
                        base_pair = base_1 + base_2

                        if pair not in stat:
                            stat[pair] = {}

                        if base_pair not in stat[pair]:
                            stat[pair][base_pair] = 0

                        stat[pair][base_pair] += 1
                except:
                    break

        # Keep positions with at least three non-singleton base-pair variants
        with open('results/B_subtilis/rna_stat.txt', 'a') as out:
            for key, dictt in stat.items():
                dictt = {pair: count for pair, count in dictt.items() if count != 1}

                if len(dictt) > 2:
                    for pair, count in sorted(dictt.items(), key=lambda x: x[1], reverse=True):
                        print(name, key, pair, count, sep='\t', file=out)
