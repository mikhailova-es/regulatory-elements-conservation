import subprocess
import os
from Bio.Blast import NCBIXML
from Bio import SeqIO

# Input paths
promoters_fna_path = 'data/E_coli/regulatory_elements/promoters.fna'

sin_minus_dir = 'intermediate_data/E_coli/intergenic_regions/ortologs_sinvergons-_align/'
sin_plus_dir = 'intermediate_data/E_coli/intergenic_regions/ortologs_sinvergons+_align/'
divergent_dir = 'intermediate_data/E_coli/intergenic_regions/ortologs_divergons_align/'
convergent_dir = 'intermediate_data/E_coli/intergenic_regions/ortologs_convergons_align/'

blast_path = 'path_to_blast'

folder_paths = [sin_minus_dir, sin_plus_dir, divergent_dir, convergent_dir]

def filter_blast_results(blast_record):
    return any(alignment.hsps for alignment in blast_record.alignments)

for folder_path in folder_paths:
    os.makedirs(folder_path + '/promoters_alignments/', exist_ok=True)
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


            # Create a BLAST database and search for promoter sequences
            cmd = blast_path + "makeblastdb -in " + input_file + " -dbtype nucl"
            subprocess.run(cmd, shell=True)

            cmd = blast_path + "blastn -query " + promoters_fna_path + " -db " + input_file + " -out blast_out -word_size 7 -evalue 1e-6 -outfmt 5 -num_threads 3 -max_target_seqs 5000"
            subprocess.run(cmd, shell=True)

            # Parse BLAST results and select the best alignment 
            with open("blast_out") as result_handle:
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
                        length = int(query.split('_')[-5])

                        # Coordinates of the proximal and distal promoter regions within the query
                        prox_end = int(query.split('_')[-1])
                        prox_start = int(query.split('_')[-2])
                        dist_start = int(query.split('_')[-4])
                        dist_end = int(query.split('_')[-3])
                        query_end = length - best_hsp.query_end

                        # Count gaps 
                        prox_gaps_left = 0
                        prox_gaps_right = 0
                        dist_gaps_left = 0
                        dist_gaps_right = 0

                        for i, letter in enumerate(best_hsp.query):
                            if letter == '-':
                                if i < prox_start:
                                    prox_gaps_left += 1
                                if i >= prox_end:
                                    prox_gaps_right += 1
                                if i < dist_start:
                                    dist_gaps_left += 1
                                if i >= dist_end:
                                    dist_gaps_right += 1

                        # Retain hits with >95% identity and at most one unaligned nucleotide at each end
                        if (query_start < 2) and (query_end < 2) and identity > 0.95:
                            if best_hsp.frame[1] == 1:
                                hits[query + '_1'] = [best_hsp.sbjct_start, best_hsp.sbjct_end, query_start, query_end, prox_gaps_left, prox_gaps_right, dist_gaps_left, dist_gaps_right]
                            if best_hsp.frame[1] == -1:
                                hits[query + '_-1'] = [best_hsp.sbjct_start, best_hsp.sbjct_end, query_start, query_end, prox_gaps_left, prox_gaps_right, dist_gaps_left, dist_gaps_right]

                    # Convert BLAST coordinates to genomic coordinates for proximal and distal promoter regions
                    for key, hit in hits.items():
                        file_name = key

                        prox_end = int(key.split('_')[-2])
                        prox_start = int(key.split('_')[-3])
                        dist_start = int(key.split('_')[-5])
                        dist_end = int(key.split('_')[-4])
                        query_length = int(key.split('_')[-6])

                        reverse_start_prox = query_length - prox_end
                        reverse_start_dist = query_length - dist_end

                        proximal_length = prox_end - prox_start
                        distal_length = dist_end - dist_start

                        # Save proximal and distal promoter regions with equal-length flanking sequences
                        with open(folder_path + '/promoters_alignments/' + 'prox_' + filename[0:-4] + '_' + file_name, 'w') as prox_site, open(folder_path + '/promoters_alignments/' + 'dist_' + filename[0:-4] + '_' + file_name, 'w') as dist_site:
                            if hit:
                                start_query = int(hit[2])
                                end_query = int(hit[3])

                                gap_left_prox = int(hit[4])
                                gap_right_prox = int(hit[5])
                                gap_left_dist = int(hit[6])
                                gap_right_dist = int(hit[7])

                                left = int(hit[0])
                                right = int(hit[1])

                                # Calculate coordinates separately for forward and reverse strand matches
                                if left < right:
                                    left_coordinate_prox = left - start_query + prox_start + gap_left_prox - 1
                                    right_coordinate_prox = right + end_query - reverse_start_prox - gap_right_prox
                                    left_coordinate_dist = left - start_query + dist_start + gap_left_dist - 1
                                    right_coordinate_dist = right + end_query - reverse_start_dist - gap_right_dist
                                else:
                                    left_coordinate_prox = right - end_query + reverse_start_prox + gap_right_prox - 1
                                    right_coordinate_prox = left + start_query - prox_start - gap_left_prox
                                    left_coordinate_dist = right - end_query + reverse_start_dist + gap_right_dist - 1
                                    right_coordinate_dist = left + start_query - dist_start - gap_left_dist

                                proximal_sequence_length = right_coordinate_prox - left_coordinate_prox
                                distal_sequence_length = right_coordinate_dist - left_coordinate_dist

                                first_record_id = list(records.keys())[0]
                                genome_length = len(records[first_record_id])


                                # Save proximal promoter sequences with flanks of the same length
                                if (left_coordinate_prox - proximal_sequence_length >= 0) and (right_coordinate_prox + proximal_sequence_length <= genome_length):
                                    for record_id, sequence in records.items():
                                        print('>' + record_id, file=prox_site)
                                        print(sequence[left_coordinate_prox - proximal_sequence_length:right_coordinate_prox + proximal_sequence_length], file=prox_site)

                                # Save distal promoter sequences with flanks of the same length
                                if (left_coordinate_dist - distal_sequence_length >= 0) and (right_coordinate_dist + distal_sequence_length <= genome_length):
                                    for record_id, sequence in records.items():
                                        print('>' + record_id, file=dist_site)
                                        print(sequence[left_coordinate_dist - distal_sequence_length:right_coordinate_dist + distal_sequence_length], file=dist_site)

                except:
                    print('xml file empty')
