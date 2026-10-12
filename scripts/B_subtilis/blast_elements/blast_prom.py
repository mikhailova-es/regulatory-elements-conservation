import subprocess
import os
from Bio.Blast import NCBIXML
from Bio import SeqIO


def filter_blast_results(blast_record):
    return any(alignment.hsps for alignment in blast_record.alignments)


# Input factor and alignment folders
promoters_fna_path = 'data/B_subtilis/regulatory_elements/promoters_cord.fna'
sin_minus_dir = 'intermediate_data/B_subtilis/intergenic_regions/ortologs_sinvergons-_align/'
sin_plus_dir = 'intermediate_data/B_subtilis/intergenic_regions/ortologs_sinvergons+_align/'
divergent_dir = 'intermediate_data/B_subtilis/intergenic_regions/ortologs_divergons_align/'
convergent_dir = 'intermediate_data/B_subtilis/intergenic_regions/ortologs_convergons_align/'
blast_path = 'path_to_blast'
folder_paths = [sin_minus_dir, sin_plus_dir, divergent_dir, convergent_dir]

for folder_path in folder_paths:
    os.makedirs(folder_path + '/promoters_alignments/', exist_ok=True)
    for filename in os.listdir(folder_path):
        if filename[-3:] == 'fna':

            # Read alignment sequences 
            with open(folder_path + '/' + filename, 'r') as file:
                records = {}
                identifiers = set()
                for record in SeqIO.parse(file, 'fasta'):
                    if record.id not in identifiers:
                        identifiers.add(record.id)
                        records[record.id] = record.seq


            # Create a BLAST database and search for promoter sequences

            cmd = blast_path + "makeblastdb -in " + folder_path + '/' + filename + " -dbtype nucl"
            subprocess.run(cmd, shell=True)
            cmd = blast_path + "blastn -query " + promoters_fna_path + " -db "  + folder_path + '/' + filename + " -out out_blast -word_size 7 -evalue 1e-6 -outfmt 5 -num_threads 3 -max_target_seqs 500"
            subprocess.run(cmd, shell=True)

            # Parse BLAST results and select the best alignment
            with open(filename + "out") as result_handle:
                blast_records = NCBIXML.parse(result_handle)
                d = {}
                try:
                    for blast_record in filter(filter_blast_results, blast_records):
                        query = blast_record.query.split(' ')[0]
                        blast_record.alignments.sort(key=lambda align: max(hsp.bits for hsp in align.hsps), reverse=True)
                        best_alignment = blast_record.alignments[0]
                        best_hsp = best_alignment.hsps[0]

                        identity = best_hsp.identities / best_hsp.align_length
                        sravn1 = best_hsp.query_start - 1
                        length = int(query.split('_')[-5])
                        prox_end = int(query.split('_')[-1]) + 1
                        prox_start = int(query.split('_')[-2])
                        dist_start = int(query.split('_')[-4])
                        dist_end = int(query.split('_')[-3]) + 1
                        sravn2 = length - best_hsp.query_end

                        # Count gaps in proximal and distal promoter regions
                        prox_gaps_left = 0
                        prox_gaps_right = 0
                        dist_gaps_left = 0
                        dist_gaps_right = 0

                        for i, letter in enumerate(best_hsp.query):
                            if letter == '-':
                                print(i, letter)
                                if i < prox_start:
                                    prox_gaps_left += 1
                                if i >= prox_end:
                                    prox_gaps_right += 1
                                if i < dist_start:
                                    dist_gaps_left += 1
                                if i >= dist_end:
                                    dist_gaps_right += 1

                        # Keep high-identity alignments with minimal unaligned query sequence
                        if (sravn1 < 2) and (sravn2 < 2) and identity > 0.95:
                            if best_hsp.frame[1] == 1:
                                d[query + '_1'] = [best_hsp.sbjct_start, best_hsp.sbjct_end, sravn1, sravn2, prox_gaps_left, prox_gaps_right, dist_gaps_left, dist_gaps_right]
                            if best_hsp.frame[1] == -1:
                                d[query + '_-1'] = [best_hsp.sbjct_start, best_hsp.sbjct_end, sravn1, sravn2, prox_gaps_left, prox_gaps_right, dist_gaps_left, dist_gaps_right]

                    # Calculate genomic coordinates and extract promoter regions with equal-length flanks
                    for key, el in d.items():
                        file_name = key
                        prox_end = int(key.split('_')[-2]) + 1
                        prox_start = int(key.split('_')[-3])
                        dist_start = int(key.split('_')[-5])
                        dist_end = int(key.split('_')[-4]) + 1
                        dlina = int(key.split('_')[-6])

                        reverse_start_prox = dlina - prox_end
                        reverse_start_dist = dlina - dist_end
                        true_leng_prox = prox_end - prox_start
                        true_leng_dist = dist_end - dist_start

                        with open(folder_path + '/promoters_alignments/' + 'prox_' + filename[0:-4] + '_' + file_name, 'w') as prox_site, open(folder_path + '/promoters_alignments/' + 'dist_' + filename[0:-4] + '_' + file_name, 'w') as dist_site:
                            if el:
                                start_query = int(el[2])
                                end_query = int(el[3])
                                gap_left_prox = int(el[4])
                                gap_right_prox = int(el[5])
                                gap_left_dist = int(el[6])
                                gap_right_dist = int(el[7])
                                left = int(el[0])
                                right = int(el[1])

                                # Convert BLAST coordinates to genomic coordinates
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

                                leng_prox = right_coordinate_prox - left_coordinate_prox
                                leng_dist = right_coordinate_dist - left_coordinate_dist

                                first_record_id = list(records.keys())[0]
                                gen_len = len(records[first_record_id])

                                # Extract proximal and distal promoter sequences
                                if (left_coordinate_prox - leng_prox >= 0) and (right_coordinate_prox + leng_prox <= gen_len):
                                    for id, seq in records.items():
                                        print('>' + id, file=prox_site)
                                        print(seq[left_coordinate_prox - leng_prox:right_coordinate_prox + leng_prox], file=prox_site)

                                if (left_coordinate_dist - leng_dist >= 0) and (right_coordinate_dist + leng_dist <= gen_len):
                                    for id, seq in records.items():
                                        print('>' + id, file=dist_site)
                                        print(seq[left_coordinate_dist - leng_dist:right_coordinate_dist + leng_dist], file=dist_site)
                except:
                    print('xml file empty')
