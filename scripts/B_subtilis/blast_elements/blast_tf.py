import subprocess
import os
from Bio.Blast import NCBIXML
from Bio import SeqIO


def filter_blast_results(blast_record):
    return any(alignment.hsps for alignment in blast_record.alignments)


# Input factor and alignment folders
tf_fna_path = 'data/B_subtilis/regulatory_elements/tf-factors.fna'
sin_minus_dir = 'intermediate_data/B_subtilis/intergenic_regions/ortologs_sinvergons-_align/'
sin_plus_dir = 'intermediate_data/B_subtilis/intergenic_regions/ortologs_sinvergons+_align/'
divergent_dir = 'intermediate_data/B_subtilis/intergenic_regions/ortologs_divergons_align/'
convergent_dir = 'intermediate_data/B_subtilis/intergenic_regions/ortologs_convergons_align/'
blast_path = 'path_to_blast'
folder_paths = [sin_minus_dir, sin_plus_dir, divergent_dir, convergent_dir]

for folder_path in folder_paths:
    os.makedirs(folder_path + '/tf_factors_alignments/', exist_ok=True)
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

            # Create a BLAST database and search for transcription factor binding sites
            cmd = blast_path + "makeblastdb -in " + folder_path + '/' + filename + " -dbtype nucl"
            subprocess.run(cmd, shell=True)
            cmd = blast_path + "blastn -query " + tf_fna_path + " -db "  + folder_path + '/' + filename + " -out out_blast -word_size 7 -evalue 1e-6 -outfmt 5 -num_threads 3 -max_target_seqs 500"
            subprocess.run(cmd, shell=True)

            # Parse BLAST results and select the best alignment
            with open("out_blast") as result_handle:
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
                        length = int(query.split('_')[-3])
                        tf_end = int(query.split('_')[-1])
                        tf_start = int(query.split('_')[-2])
                        sravn2 = length - best_hsp.query_end
                        gaps_left = 0
                        gaps_right = 0

                        # Count gaps 
                        for i, letter in enumerate(best_hsp.query):
                            if letter == '-':
                                if i < tf_start:
                                    gaps_left += 1
                                if i >= tf_end:
                                    gaps_right += 1

                        # Keep high-identity alignments with minimal unaligned query sequence
                        if (sravn1 < 2) and (sravn2 < 2) and identity > 0.95:
                            if best_hsp.frame[1] == 1:
                                d[query + '_1'] = [best_hsp.sbjct_start, best_hsp.sbjct_end, sravn1, sravn2, gaps_left, gaps_right]
                            if best_hsp.frame[1] == -1:
                                d[query + '_-1'] = [best_hsp.sbjct_start, best_hsp.sbjct_end, sravn1, sravn2, gaps_left, gaps_right]

                    # Calculate genomic coordinates and extract the binding site with equal-length flanks
                    for key, el in d.items():
                        file_name = key
                        dlina = int(key.split('_')[-4])
                        tf_start = int(key.split('_')[-3])
                        tf_end = int(key.split('_')[-2])
                        reverse_start = dlina - tf_end

                        with open(folder_path + '/tf_factors_alignments/' + filename[0:-4] + '_' + file_name, 'w') as site:

                            if el:
                                start_query = int(el[2])
                                end_query = int(el[3])
                                gap_left = int(el[4])
                                gap_right = int(el[5])
                                left = int(el[0])
                                right = int(el[1])

                                if left < right:
                                    left_coordinate = left - start_query + tf_start + gap_left - 1
                                    right_coordinate = right + end_query - reverse_start - gap_right
                                else:
                                    left_coordinate = right - end_query + reverse_start + gap_right - 1
                                    right_coordinate = left + start_query - tf_start - gap_left

                                leng = right_coordinate - left_coordinate
                                first_record_id = list(records.keys())[0]
                                gen_len = len(records[first_record_id])

                                # Extract the binding site with equal-length flanks
                                if (left_coordinate - leng >= 0) and (right_coordinate + leng <= gen_len):
                                    for id, seq in records.items():
                                        print('>' + id, file=site)
                                        print(seq[left_coordinate - leng:right_coordinate + leng], file=site)
                except:
                    print('xml file empty')
