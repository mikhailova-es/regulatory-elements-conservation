import subprocess
import os
from Bio.Blast import NCBIXML
from Bio import SeqIO

# Input paths
rnas_fna_path = 'data/E_coli/regulatory_elements/rnas.fna'

sin_minus_dir = 'intermediate_data/E_coli/intergenic_regions/ortologs_sinvergons-_align/'
sin_plus_dir = 'intermediate_data/E_coli/intergenic_regions/ortologs_sinvergons+_align/'
divergent_dir = 'intermediate_data/E_coli/intergenic_regions/ortologs_divergons_align/'
convergent_dir = 'intermediate_data/E_coli/intergenic_regions/ortologs_convergons_align/'

blast_path = 'path_to_blast'

folder_paths = [sin_minus_dir, sin_plus_dir, divergent_dir, convergent_dir]

def filter_blast_results(blast_record):
    return any(alignment.hsps for alignment in blast_record.alignments)


for folder_path in folder_paths:
    for filename in os.listdir(folder_path):
        if filename[-3:] == 'fna':

            # Read alignment sequences
            input_file = os.path.join(folder_path, filename)
            with open(input_file, 'r') as file:
                records = {}
                identifiers = set()
                for record in SeqIO.parse(file, 'fasta'):
                    if record.id not in identifiers:
                        identifiers.add(record.id)
                        records[record.id] = record.seq

            # Create a BLAST database from the intergenic region sequences
            cmd = blast_path + "makeblastdb -in " + input_file + " -dbtype nucl"
            subprocess.run(cmd, shell=True)

            # Search for Rfam RNA sequences in the intergenic regions
            cmd = blast_path + "blastn -query " + rnas_fna_path + " -db " + input_file + " -out out_blast -word_size 7 -evalue 1e-6 -outfmt 5 -num_threads 3 -max_target_seqs 5000"
            subprocess.run(cmd, shell=True)

            # Parse BLAST results and select the best alignment for each query
            with open("out_blast") as result_handle:
                blast_records = NCBIXML.parse(result_handle)
                hits = {}
                try:
                    for blast_record in filter(filter_blast_results, blast_records):
                        query = blast_record.query.split(' ')[0]
                        blast_record.alignments.sort(key=lambda align: max(hsp.bits for hsp in align.hsps), reverse=True)
                        best_alignment = blast_record.alignments[0]
                        best_hsp = best_alignment.hsps[0]

                        identity = best_hsp.identities / best_hsp.align_length
                        query_start = best_hsp.query_start - 1
                        query_length = int(query.split('_')[-1])
                        query_end = query_length - best_hsp.query_end

                        # Retain hits with >95% identity and at most one unaligned nucleotide at each end
                        if (query_start < 2) and (query_end < 2) and identity > 0.95:
                            if best_hsp.frame[1] == 1:
                                hits[query + '_1'] = [best_hsp.sbjct_start, best_hsp.sbjct_end, query_start, query_end]
                            if best_hsp.frame[1] == -1:
                                hits[query + '_-1'] = [best_hsp.sbjct_start, best_hsp.sbjct_end, query_start, query_end]

                    # Calculate genomic coordinates and extract the RNA regions
                    for key, hit in hits.items():
                        file_name = key

                        # Open separate outputs for RNA only, both flanks, and each individual flank
                        with open(folder_path + '/rnas_alignments/both_flanks/' + filename[0:-4] + '_' + file_name, 'w') as both_flanks, open(folder_path + '/rnas_alignments/without_flanks/' + filename[0:-4] + '_' + file_name, 'w') as rna_only, open(folder_path + '/rnas_alignments/right_flank/' + filename[0:-4] + '_' + file_name, 'w') as right_flank, open(folder_path + '/rnas_alignments/left_flank/' + filename[0:-4] + '_' + file_name, 'w') as left_flank:
                            if hit:
                                start_query = int(hit[2])
                                end_query = int(hit[3])
                                left = int(hit[0])
                                right = int(hit[1])

                                first_record_id = list(records.keys())[0]
                                genome_length = len(records[first_record_id])

                                # Convert BLAST coordinates to genomic coordinates for forward and reverse matches
                                if left < right:
                                    left_coordinate = left - start_query - 1
                                    right_coordinate = right + end_query
                                    element_length = right_coordinate - left_coordinate
                                else:
                                    left_coordinate = right - end_query - 1
                                    right_coordinate = left + start_query
                                    element_length = right_coordinate - left_coordinate

                                # Save the RNA region without flanking sequences
                                if (left_coordinate >= 0) and (right_coordinate <= genome_length):
                                    for record_id, sequence in records.items():
                                        print('>' + record_id, file=rna_only)
                                        print(sequence[left_coordinate:right_coordinate], file=rna_only)

                                # Save the RNA region with equal-length flanks on both sides
                                if (left_coordinate - element_length >= 0) and (right_coordinate + element_length <= genome_length):
                                    for record_id, sequence in records.items():
                                        print('>' + record_id, file=both_flanks)
                                        print(sequence[left_coordinate - element_length:right_coordinate + element_length], file=both_flanks)
                                    continue

                                # If only the right flank fits, save the RNA region with the right flank
                                if (left_coordinate >= 0) and (right_coordinate + element_length <= genome_length):
                                    for record_id, sequence in records.items():
                                        print('>' + record_id, file=right_flank)
                                        print(sequence[left_coordinate:right_coordinate + element_length], file=right_flank)

                                # If only the left flank fits, save the RNA region with the left flank
                                if (left_coordinate - element_length >= 0) and (right_coordinate <= genome_length):
                                    for record_id, sequence in records.items():
                                        print('>' + record_id, file=left_flank)
                                        print(sequence[left_coordinate - element_length:right_coordinate], file=left_flank)

                except:
                    print('xml file empty')