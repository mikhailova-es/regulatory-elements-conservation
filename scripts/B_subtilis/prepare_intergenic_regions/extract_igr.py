import os
import gffutils

# Input and output paths
genomes_dir = 'data/B_subtilis/genomes/'
proteinortho_file = 'data/B_subtilis/orthologous_groups/proteinortho.tsv'
orthologous_annotation_dir = 'intermediate_data/B_subtilis/ortho_annotation/'
sin_minus_dir = 'intermediate_data/B_subtilis/intergenic_regions/ortologs_sinvergons-/'
sin_plus_dir = 'intermediate_data/B_subtilis/intergenic_regions/ortologs_sinvergons+/'
divergent_dir = 'intermediate_data/B_subtilis/intergenic_regions/ortologs_divergons/'
convergent_dir = 'intermediate_data/B_subtilis/intergenic_regions/ortologs_convergons/'

desired_types = {'CDS', 'tRNA', 'rRNA', 'tmRNA', 'RNase_P_RNA', 'SRP_RNA', 'pseudogene'}


# Collect GFF annotation files from genome directories
genome_dirs = os.listdir(genomes_dir)
annotation_files = []

for genome_dir in genome_dirs:
    if genome_dir.startswith('GCF'):
        genome_path = os.path.join(genomes_dir, genome_dir)
        for filename in os.listdir(genome_path):
            if filename.endswith('gff'):
                annotation_files.append(os.path.join(genome_path, filename))


# Parse genome annotations and store feature coordinates and strand information
annotation_dict = {}

for annotation_file in annotation_files:
    genome_name = annotation_file.split('/')[2]
    annotation_dict[genome_name] = {}

    db = gffutils.create_db(annotation_file, dbfn=':memory:', force=True, merge_strategy='create_unique')
    main_seqid = next(db.all_features()).seqid
    all_features = db.all_features()

    for feature in all_features:
        if feature.featuretype in desired_types:
            if feature.seqid != main_seqid:
                continue

            cds = feature['ID'][0][4:]
            genome_cds = main_seqid + '-' + cds
            annotation_dict[genome_name][genome_cds] = [feature.featuretype, feature.start, feature.end, feature.strand]


# Add orthologous group names from ProteinOrtho results

def find_name(row):
    for j in range(len(row)):
        if row[j] != '*':
            return j
        
with open(proteinortho_file, 'r') as inp:
    header = inp.readline().strip().split('\t')[3:]
    line = inp.readline().strip().split('\t')[3:]

    while line:
        for i in range(len(line)):
            if line[i] != '*':
                if len(line[i].split(',')) != 1:
                    cds = line[i].split(',')[0]
                else:
                    cds = line[i]

                gcf = header[i]
                ortholog_position = find_name(line)
                ortholog_name = line[ortholog_position]

                if len(ortholog_name.split(',')) != 1:
                    ortholog_name = ortholog_name.split(',')[0]

                try:
                    annotation_dict[gcf][cds].append(ortholog_name)
                except:
                    continue

        line = inp.readline().strip().split('\t')[3:]


# Write genome annotations with orthologous group names
for genome_ac, annotations in annotation_dict.items():
    with open(orthologous_annotation_dir + genome_ac, 'w') as out:
        for cds, value in annotations.items():
            if len(value) == 4:
                print(cds, cds, value[0], value[1], value[2], value[3], sep='\t', file=out)
            else:
                print(value[4], cds, value[0], value[1], value[2], value[3], sep='\t', file=out)

#Get genome sequence
def get_genome_seq(genome):
    genome_files = [ os.path.join(genomes_dir + genome, filename) for filename in os.listdir(genomes_dir + genome) if filename.startswith('GCF')]

    with open(genome_files[0], 'r') as fna:
        sequence = ''
        line_fna = fna.readline().strip()

        while True:
            line_fna = fna.readline().strip()
            if not line_fna or line_fna.startswith('>'):
                break
            sequence += line_fna

    return sequence


# Extract intergenic sequences between consecutive CDSs
orthologous_annotation_files = os.listdir(orthologous_annotation_dir)

for genome_acc in orthologous_annotation_files:
    sequence = get_genome_seq(genome_acc)

    with open(orthologous_annotation_dir + genome_acc, mode='r') as orthologs_file:
        orthologs_list = []
        for line in orthologs_file:
            orthologs_list.append(line.strip())

    for j in range(len(orthologs_list) - 1):
        strand_1 = orthologs_list[j].split('\t')[5]
        strand_2 = orthologs_list[j + 1].split('\t')[5]
        coordinate_start = int(orthologs_list[j].split('\t')[4])
        coordinate_end = int(orthologs_list[j + 1].split('\t')[3])
        ortholog_name_1 = orthologs_list[j].split('\t')[0]
        ortholog_name_2 = orthologs_list[j + 1].split('\t')[0]
        feature_type_1 = orthologs_list[j].split('\t')[2]
        feature_type_2 = orthologs_list[j + 1].split('\t')[2]

        if feature_type_1 == 'CDS' and feature_type_2 == 'CDS':

            if strand_1 == '-' and strand_2 == '-':
                with open(sin_minus_dir + ortholog_name_1 + '-' + ortholog_name_2 + '.fna', 'a') as sin1:
                    sin1.write('>' + genome_acc + '\n')
                    sin1.write(sequence[coordinate_start:coordinate_end - 1] + '\n')

            if strand_1 == '+' and strand_2 == '+':
                with open(sin_plus_dir + ortholog_name_1 + '-' + ortholog_name_2 + '.fna', 'a') as sin2:
                    sin2.write('>' + genome_acc + '\n')
                    sin2.write(sequence[coordinate_start:coordinate_end - 1] + '\n')

            if strand_1 == '-' and strand_2 == '+':
                with open(divergent_dir + ortholog_name_1 + '-' + ortholog_name_2 + '.fna', 'a') as div:
                    div.write('>' + genome_acc + '\n')
                    div.write(sequence[coordinate_start:coordinate_end - 1] + '\n')

            if strand_1 == '+' and strand_2 == '-':
                with open(convergent_dir + ortholog_name_1 + '-' + ortholog_name_2 + '.fna', 'a') as con:
                    con.write('>' + genome_acc + '\n')
                    con.write(sequence[coordinate_start:coordinate_end - 1] + '\n')
