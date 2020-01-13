from Bio import SeqIO
#newfile = open("Zea_mays.AGPv3.21.pep.longest.fa", "w")
records = list(SeqIO.parse("Zea_mays.AGPv3.21.pep.all.fa", "fasta"))


gene_names = set()
iso_len_dict = {}
i = 0
# This could probably be made much faster if I do the length comparisons here
# Within a dictionary, and then just write the dictionary to a fasta
# But this is not implemented yet. Takes about 15 minutes to run as-is...
while i <= len(records)-1:
#while i <= 10:
    name = records[i].name
    if name[0]=="G":
    	gene_names.add(name[:13])
    elif name[0]=="A" or name[0]=="E": 
    	gene_names.add(name) # I'm note shortening their names
    else:
    	print(name) # in case I didn't miss any genes
    iso_len_dict[name] = len(records[i].seq)
    i = i + 1

gene_dict = {i : -1 for i in gene_names} # creates an empty dictionary



for gene in gene_names:
	temp_list = []
	[temp_list.append([key, value]) for key, value in iso_len_dict.items() if gene in key]
	if len(temp_list)==1: # if there's only 1 isoform, add as longest
		gene_dict[gene]=temp_list[0][0]
		continue
	temp_list = list(map(list, zip(*temp_list)))
	max_value = max(temp_list[1])
	max_index = temp_list[1].index(max_value)
	gene_dict[gene]=temp_list[0][max_index]

iso_list = (list(gene_dict.values()))
gene_list = (list(gene_dict.keys()))

with open("Zea_mays.AGPv3.21.pep.longest.fa","w") as f:
	for seq_record in SeqIO.parse("Zea_mays.AGPv3.21.pep.all.fa", "fasta"):
		try:
			list_index = iso_list.index(seq_record.name)
			f.write(str(">"+gene_list[list_index]) + "\n")
			f.write(str(seq_record.seq) + "\n") 
		except:
			None

#print(gene_names)
#print(iso_names)
#print(seq_lens)