import subprocess
import os
from Bio.Blast import NCBIXML
from Bio import SeqIO


def filter_blast_results(blast_record):
    return any(alignment.hsps for alignment in blast_record.alignments)


# Input paths
term_fna_path = 'data/E_coli/regulatory_elements/terminators.fna'
sin_minus_dir = 'intermediate_data/E_coli/intergenic_regions/ortologs_sinvergons-_align/'
sin_plus_dir = 'intermediate_data/E_coli/intergenic_regions/ortologs_sinvergons+_align/'
divergent_dir = 'intermediate_data/E_coli/intergenic_regions/ortologs_divergons_align/'
convergent_dir = 'intermediate_data/E_coli/intergenic_regions/ortologs_convergons_align/'

blast_path = 'path_to_blast'

folder_paths = [sin_minus_dir, sin_plus_dir, divergent_dir, convergent_dir]

for folder_path in folder_paths:
    os.makedirs(folder_path + '/terminators_alignments/', exist_ok=True)
    os.makedirs(folder_path + '/terminators_alignments/without_flanks', exist_ok=True)
    os.makedirs(folder_path + '/terminators_alignments/both_flanks', exist_ok=True)
    os.makedirs(folder_path + '/terminators_alignments/right_flank', exist_ok=True)
    for filename in os.listdir(folder_path):
        if filename[-3:] != 'fna':
            continue

        input_file = folder_path + '/' + filename

        # Read alignment sequuences
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

        # Search for terminator sequences in the intergenic regions
        cmd = blast_path + '/blastn -query ' + term_fna_path + ' -db ' + input_file + ' -out out_blast -word_size 7 -evalue 1e-6 -outfmt 5 -num_threads 3 -max_target_seqs 5000'
        subprocess.run(cmd, shell=True)

        # Parse BLAST results and select the best alignment for each query
        with open('out_blast') as result_handle:
            blast_records = NCBIXML.parse(result_handle)
            results = {}

            try:
                for blast_record in filter(filter_blast_results, blast_records):
                    query = blast_record.query.split(' ')[0]
                    blast_record.alignments.sort(key=lambda align: max(hsp.bits for hsp in align.hsps), reverse=True)
                    best_alignment = blast_record.alignments[0]
                    best_hsp = best_alignment.hsps[0]

                    # Check sequence identity and the number of unaligned nucleotides at each end
                    identity = best_hsp.identities / best_hsp.align_length
                    left_mismatch = best_hsp.query_start - 1
                    element_length = int(query.split('_')[-1])
                    element_end = element_length - 10
                    element_start = 10
                    right_mismatch = element_length - best_hsp.query_end
                    gaps_left = 0
                    gaps_right = 0

                    # Count gaps in the 10-bp flanking regions included in the query
                    for index, letter in enumerate(best_hsp.query):
                        if letter == '-':
                            if index < element_start:
                                gaps_left += 1
                            if index >= element_end:
                                gaps_right += 1

                    # Retain hits with high identity and no more than one unaligned nucleotide at either end
                    if left_mismatch < 2 and right_mismatch < 2 and identity > 0.95:
                        if best_hsp.frame[1] == 1:
                            results[query + '_1'] = [best_hsp.sbjct_start, best_hsp.sbjct_end, left_mismatch, right_mismatch, gaps_left, gaps_right]
                        if best_hsp.frame[1] == -1:
                            results[query + '_-1'] = [best_hsp.sbjct_start, best_hsp.sbjct_end, left_mismatch, right_mismatch, gaps_left, gaps_right]

                # Calculate genomic coordinates of the mapped terminator
                for query_name, result in results.items():
                    element_length = int(query_name.split('_')[-2])
                    element_start = 10
                    element_end = element_length - 10
                    reverse_start = 10

                    # Prepare separate files for the terminator, terminator with flanks, and downstream region
                    with open(folder_path + '/terminators_alignments/both_flanks/' + filename[0:-4] + '_' + query_name, 'w') as site, open(folder_path + '/terminators_alignments/without_flanks/' + filename[0:-4] + '_' + query_name, 'w') as site2, open(folder_path + '/terminators_alignments/right_flank/' + filename[0:-4] + '_' + query_name, 'w') as site3:
                        if result:
                            start_mismatch = int(result[2])
                            end_mismatch = int(result[3])
                            gap_left = int(result[4])
                            gap_right = int(result[5])
                            left = int(result[0])
                            right = int(result[1])
                            first_record_id = list(records.keys())[0]
                            genome_length = len(records[first_record_id])

                            # Calculate coordinates for hits on the forward strand
                            if left < right:
                                left_coordinate = left - start_mismatch + element_start + gap_left - 1
                                right_coordinate = right + end_mismatch - reverse_start - gap_right
                                element_length_in_genome = right_coordinate - left_coordinate
                                print(element_length_in_genome, element_end - element_start)
                                print(records[first_record_id][left_coordinate:right_coordinate])

                                # Save the terminator and its downstream region
                                if left_coordinate >= 0 and right_coordinate + element_length_in_genome <= genome_length:
                                    for record_id, sequence in records.items():
                                        print('>' + record_id, file=site3)
                                        print(sequence[left_coordinate:right_coordinate + element_length_in_genome], file=site3)

                            # Calculate coordinates for hits on the reverse strand
                            else:
                                left_coordinate = right - end_mismatch + reverse_start + gap_right - 1
                                right_coordinate = left + start_mismatch - element_start - gap_left
                                element_length_in_genome = right_coordinate - left_coordinate
                                print(records[first_record_id][left_coordinate:right_coordinate])
                                print(element_length_in_genome, element_end - element_start)

                                # Save the terminator and its downstream region 
                                if left_coordinate - element_length_in_genome >= 0 and right_coordinate <= genome_length:
                                    for record_id, sequence in records.items():
                                        print('>' + record_id, file=site3)
                                        print(sequence[left_coordinate - element_length_in_genome:right_coordinate], file=site3)

                            # Save the terminator with equal-length flanking regions
                            if left_coordinate - element_length_in_genome >= 0 and right_coordinate + element_length_in_genome <= genome_length:
                                for record_id, sequence in records.items():
                                    print('>' + record_id, file=site)
                                    print(sequence[left_coordinate - element_length_in_genome:right_coordinate + element_length_in_genome], file=site)

                            # Save the terminator sequence without flanking regions
                            if left_coordinate >= 0 and right_coordinate <= genome_length:
                                for record_id, sequence in records.items():
                                    print('>' + record_id, file=site2)
                                    print(sequence[left_coordinate:right_coordinate], file=site2)

            except:
                print('xml file empty')
