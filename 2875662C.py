

import re
import sys

#creating a function to get the intron start and end position ny processing the CIGAR string
def get_cigar(cigar, pos):
    try:
        #Initialise the split read position as split_start and empty junctions list
        split_start= int(pos)
        junctions= []
        
        #Iterate over CIGAR string to find introns
        match= re.finditer(r'(\d+)([MDN])', cigar)
        
        for i in match:
            number= int(i.group(1))
            symbol= i.group(2)
            
            #getting the intron start and end positions by updating split_start for match and deletion, append to the junctions list
            if symbol in ['M', 'D']:
                split_start += number
            elif symbol == 'N':
                intron_start = split_start
                intron_end = split_start + number
                junctions.append((intron_start, intron_end))
        return junctions
    except ValueError:
        print(f'ValueError in get_cigar.\n')
        return []
    except Exception:
        print(f'Error in get_cigar.\n')
        return []
    

#Check if correct number of command line arguments is provided    
if len(sys.argv) != 3:
    print("Usage: python script.py sam_file gene_location_summary_file")
    sys.exit(1)

#Assign input file paths from command line arguments
sam_file_path = sys.argv[1]
gene_file_path = sys.argv[2]

#creating a dictionary for number of reads for each splice junction
intron_reads= {}

try:
    #reading the SAM file
    with open(sam_file_path) as sam_file:
        for line in sam_file:
            if line.startswith('@'):
                continue

            #Parsing the sam file
            line= line.rstrip()
            line= line.split('\t')
            rname = line[2]
            pos = int(line[3])
            cigar= line[5]
            read_times = line[-1]

            #check if there is an intron in the read and if the read is only mapped once
            if 'N' in cigar and re.search(r'1$', read_times):
                junctions= get_cigar(cigar, pos)
                        
                #Process each intron in all the reads 
                for intron_start, intron_end in junctions:
                    try: 
                        #reading the gene location file           
                        with open(gene_file_path) as gene_file:
                                    header= next(gene_file)
                                    try:

                                        #processing the lines in gene file to get the gene boundaries
                                        for gene_line in gene_file:
                                            gene_line= gene_line.rstrip()
                                            gene_line= gene_line.split()
                                            gene_id= gene_line[0]
                                            gene_loc_ch= gene_line[2]
                                            gene_loc_ch= gene_loc_ch.split(':')
                                            gene_loc= gene_loc_ch[1]
                                            gene_loc= gene_loc.split('..')
                                            gene_loc[1] = gene_loc[1].rstrip('(+)').rstrip('(-)')
                                            gene_start = int(gene_loc[0].replace(',', ''))
                                            gene_end = int(gene_loc[1].replace(',', ''))


                                            #creating a key for the junction_reads dictionary
                                            intron_key= f'{gene_id}\t{intron_start}\t{intron_end}'


                                            #check if the intron junctions are within the gene boundaries
                                            if gene_start <= intron_start <= gene_end and gene_start <= intron_end <= gene_end:

                                                #updating the junction_reads dictionary with the number of reads
                                                if intron_key in intron_reads:
                                                    intron_reads[intron_key] +=1
                                                else:
                                                    intron_reads[intron_key] = 1
                                    except Exception:
                                        print(f'Error in processing the junctions within the gene boundary: {gene_line}\n')
                    except FileNotFoundError:
                        print(f'File {gene_file_path} not found. Please check and try again.\n')
                        

except FileNotFoundError:
    print(f'File {sam_file_path} not found. Please check and try again.\n')

#Writing to the output file
with open('2875662C.txt', 'w') as sam_output:
    header = 'gene_ID\t' + 'Junction_start\t' + 'Junction_end\t' + 'Number_of_reads\n'
    sam_output.write(header)

    #initialise the current gene
    current_gene = None
    for key, num_of_reads in intron_reads.items():
       
       # Check if the current gene is different from the previous one
        if current_gene is not None and key.split('\t')[0] != current_gene:

            # Insert a newline after every block of entries with the same gene ID
            sam_output.write('\n')  

        sam_output.write(f'{key}\t{num_of_reads}\n')
        current_gene = key.split('\t')[0]
