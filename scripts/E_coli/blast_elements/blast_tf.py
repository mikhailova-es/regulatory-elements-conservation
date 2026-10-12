import subprocess
import os
from Bio.Blast import NCBIXML
from Bio import SeqIO


def filter_blast_results(blast_record):
    return any(alignment.hsps for alignment in blast_record.alignments)


# Input paths
tf_fna_path = 'data/E_coli/regulatory_elements/tf-factors.fna'
sin_minus_dir = 'intermediate_data/E_coli/intergenic_regions/ortologs_sinvergons-_align/'
sin_plus_dir = 'intermediate_data/E_coli/intergenic_regions/ortologs_sinvergons+_align/'
divergent_dir = 'intermediate_data/E_coli/intergenic_regions/ortologs_divergons_align/'
convergent_dir = 'intermediate_data/E_coli/intergenic_regions/ortologs_convergons_align/'

blast_path = 'path_to_blast'
folder_paths = [sin_minus_dir, sin_plus_dir, divergent_dir, convergent_dir]


for folder_path in folder_paths:
    os.makedirs( folder_path + '/tf_factors_alignments/', exist_ok=True)
    for filename in os.listdir(folder_path):
        if filename[-3:] != 'fna':
            continue

        input_file = folder_path + '/' + filename

        with open(input_file, 'r') as file:
            records = {}
            identifiers = set()

            for record in SeqIO.parse(file, 'fasta'):
                if record.id not in identifiers:
                    identifiers.add(record.id)
                    records[record.id] = record.seq

        # Create a BLAST database from the intergenic region sequences
        cmd = blast_path + '/makeblastdb -in ' + input_file + ' -dbtype nucl'
        subprocess.run(cmd, shell=True)

        # Search for transcription factor binding sites
        cmd = blast_path + '/blastn -query ' + tf_fna_path + ' -db ' + input_file + ' -out blast_out -word_size 7 -evalue 1e-6 -outfmt 5 -num_threads 3 -max_target_seqs 5000'
        subprocess.run(cmd, shell=True)

        # Parse BLAST results and select the best hit for each query
        with open('blast_out') as result_handle:
            blast_records = NCBIXML.parse(result_handle)
            results = {}

            try:
                for blast_record in filter(filter_blast_results, blast_records):
                    query = blast_record.query.split(' ')[0]
                    blast_record.alignments.sort(key=lambda align: max(hsp.bits for hsp in align.hsps), reverse=True)
                    best_alignment = blast_record.alignments[0]
                    best_hsp = best_alignment.hsps[0]

                    identity = best_hsp.identities / best_hsp.align_length
                    left_mismatch = best_hsp.query_start - 1
                    element_length = int(query.split('_')[-1])
                    element_start = 10
                    element_end = element_length - 10
                    right_mismatch = element_length - best_hsp.query_end
                    gaps_left = 0
                    gaps_right = 0

                    for index, letter in enumerate(best_hsp.query):
                        if letter == '-':
                            if index < element_start:
                                gaps_left += 1
                            if index >= element_end:
                                gaps_right += 1

                    if left_mismatch < 2 and right_mismatch < 2 and identity > 0.95:
                        if best_hsp.frame[1] == 1:
                            results[query + '_1'] = [best_hsp.sbjct_start, best_hsp.sbjct_end, left_mismatch, right_mismatch, gaps_left, gaps_right]
                        if best_hsp.frame[1] == -1:
                            results[query + '_-1'] = [best_hsp.sbjct_start, best_hsp.sbjct_end, left_mismatch, right_mismatch, gaps_left, gaps_right]

                # Extract the mapped binding site and its flanking regions
                for query_name, result in results.items():
                    element_length = int(query_name.split('_')[-2])
                    element_start = 10
                    element_end = element_length - 10
                    reverse_start = 10

                    with open(folder_path + '/tf_factors_alignments/' + filename[0:-4] + '_' + query_name, 'w') as site:
                        if result:
                            start_mismatch = int(result[2])
                            end_mismatch = int(result[3])
                            gap_left = int(result[4])
                            gap_right = int(result[5])
                            left = int(result[0])
                            right = int(result[1])

                            if left < right:
                                left_coordinate = left - start_mismatch + element_start + gap_left - 1
                                right_coordinate = right + end_mismatch - reverse_start - gap_right
                            else:
                                left_coordinate = right - end_mismatch + reverse_start + gap_right - 1
                                right_coordinate = left + start_mismatch - element_start - gap_left

                            element_length_in_genome = right_coordinate - left_coordinate

                            first_record_id = list(records.keys())[0]
                            genome_length = len(records[first_record_id])

                            print(records[first_record_id][left_coordinate:right_coordinate])

                            if left_coordinate - element_length_in_genome >= 0 and right_coordinate + element_length_in_genome <= genome_length:
                                for record_id, sequence in records.items():
                                    print('>' + record_id, file=site)
                                    print(sequence[left_coordinate - element_length_in_genome:right_coordinate + element_length_in_genome], file=site)

            except:
                print('xml file empty')
