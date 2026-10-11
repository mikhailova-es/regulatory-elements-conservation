import os
import gffutils

# Input and output paths
genomes_dir = 'data/E_coli/genomes/'
panaroo_file = 'data/E_coli/orthologous_groups/gene_presence_absence.csv'
annotation_dir = 'intermediate_data/E_coli/ortho_annotation/'
igr_dir = 'intermediate_data/E_coli/intergenic_regions/'

# Create output directories

for folder in [
    'ortologs_sinvergons-',
    'ortologs_sinvergons+',
    'ortologs_divergons',
    'ortologs_convergons'
]:
    os.makedirs(os.path.join(igr_dir, folder), exist_ok=True)
    
# Find GFF annotation files for all genomes
spis = os.listdir(genomes_dir)
annotation_files = []
desired_types = {'CDS', 'tRNA', 'rRNA', 'tmRNA', 'RNase_P_RNA', 'SRP_RNA', 'pseudogene'}

for i in spis:
    if i.startswith('GCF'):
        path1 = genomes_dir + i
        local_spis = os.listdir(path1)
        for j in local_spis:
            if j.endswith('gff'):
                annotation_files.append(path1 + '/' + j)

# Parse GFF annotations and store feature coordinates and strand
ann_dict = {}
gcf_gff = {}

for file in annotation_files:
    name = file.split('/')[3]
    ann_dict[name] = {}
    db = gffutils.create_db(file, dbfn=':memory:', force=True, merge_strategy='create_unique')
    all_features = db.all_features()
    main_seqid = next(db.all_features()).seqid
    gcf_gff[name] = main_seqid

    for seq in all_features:
        if seq.featuretype in desired_types:
            if seq.seqid != main_seqid:
                continue
            cds = seq['ID'][0][4:]
            genome_cds = main_seqid + '-' + cds
            ann_dict[name][genome_cds] = [seq.featuretype, seq.start, seq.end, seq.strand]

# Add Panaroo orthologous group names to the corresponding genome annotations
with open(panaroo_file, 'r') as inp:
    header = inp.readline().strip().split(',')[3:]
    linee = inp.readline().strip().split(',')
    name = linee[0]
    line = linee[3:]

    while line:
        for el in range(len(line)):
            if line[el] != '':
                if line[el].split('-')[0] != 'cds':
                    continue
                cds_id = line[el].split('-')[1]
                gcf = header[el]
                cds_genome = gcf_gff[gcf]
                cds = cds_genome + '-' + cds_id
                orto_name = name
                try:
                    ann_dict[gcf][cds].append(orto_name)
                except:
                    continue

        linee = inp.readline().strip().split(',')
        name = linee[0]
        line = linee[3:]

# Save the annotations of each genome
for key, inner_dict in ann_dict.items():
    with open(annotation_dir + key, 'w') as out:
        for cds, value in inner_dict.items():
            if len(value) == 4:
                print(cds, cds, value[0], value[1], value[2], value[3], sep='\t', file=out)
            else:
                print(value[4], cds, value[0], value[1], value[2], value[3], sep='\t', file=out)

# Read the сhromosome sequence from a genome FASTA file
def get_genome_seq(file):
    genome_file = [os.path.join(genomes_dir + file, k) for k in os.listdir(genomes_dir + file) if k.startswith('GCF')]
    with open(genome_file[0], 'r') as fna:
        sequence = ''
        line = fna.readline().strip()
        while True:
            line = fna.readline().strip()
            if not line or line.startswith('>'):
                break
            sequence += line
    return sequence

# Extract intergenic regions between adjacent CDSs and group them by the orientation of the two CDSs
ortho_ann_result = os.listdir(annotation_dir)

for file in ortho_ann_result:
    sequence = get_genome_seq(file)

    with open(annotation_dir + file, mode='r') as ortologs_file:
        ortologs_spisok = []
        for line_ort in ortologs_file:
            ortologs_spisok.append(line_ort.strip())

    for j in range(len(ortologs_spisok) - 1):
        cep_1 = ortologs_spisok[j].split('\t')[5]
        cep_2 = ortologs_spisok[j + 1].split('\t')[5]
        cord_start = int(ortologs_spisok[j].split('\t')[4])
        cord_end = int(ortologs_spisok[j + 1].split('\t')[3])
        name_ort_1 = ortologs_spisok[j].split('\t')[0]
        name_ort_2 = ortologs_spisok[j + 1].split('\t')[0]
        typee1 = ortologs_spisok[j].split('\t')[2]
        typee2 = ortologs_spisok[j + 1].split('\t')[2]

        # Only extract regions located between two CDSs
        if typee1 == 'CDS' and typee2 == 'CDS':
            if cep_1 == '-' and cep_2 == '-':
                with open(igr_dir + 'ortologs_sinvergons-/' + name_ort_1 + '-' + name_ort_2 + '.fna', 'a') as sin1:
                    sin1.write('>' + file + '\n')
                    sin1.write(sequence[cord_start:cord_end - 1] + '\n')

            if cep_1 == '+' and cep_2 == '+':
                with open(igr_dir + 'ortologs_sinvergons+/' + name_ort_1 + '-' + name_ort_2 + '.fna', 'a') as sin2:
                    sin2.write('>' + file + '\n')
                    sin2.write(sequence[cord_start:cord_end - 1] + '\n')

            if cep_1 == '-' and cep_2 == '+':
                with open(igr_dir + 'ortologs_divergons/' + name_ort_1 + '-' + name_ort_2 + '.fna', 'a') as div:
                    div.write('>' + file + '\n')
                    div.write(sequence[cord_start:cord_end - 1] + '\n')

            if cep_1 == '+' and cep_2 == '-':
                with open(igr_dir + 'ortologs_convergons/' + name_ort_1 + '-' + name_ort_2 + '.fna', 'a') as con:
                    con.write('>' + file + '\n')
                    con.write(sequence[cord_start:cord_end - 1] + '\n')
