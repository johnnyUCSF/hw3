from Bio import SeqIO
#########
def import_blosum(filename):
    #######
    with open(filename) as f:
        content = f.readlines()
        # you may also want to remove whitespace characters like `\n` at the end of each line
        content = [x.strip() for x in content]
    #######
    blosum = []
    firstline = True
    for line in content:
        if line[0] != '#':
            if firstline == True:
                row = line.split()
                blosum.append(row)
                firstline = False
            else:
                row = line.split()
                row = [float(i) for i in row]
                blosum.append(row)
    return(blosum)

def import_pairs(filename):
    #####
    with open(filename) as f:
        content = f.readlines()
        content = [x.strip() for x in content]
    #####
    pairs = []
    for line in content:
        ###
        row = []
        ###
        fasta_files = line.split()
        for file in fasta_files:
            row.append(import_fasta(file))
        ###
        pairs.append(row)
    #####
    return(pairs)

def import_fasta(filename):
    for record in SeqIO.parse(filename, "fasta"):
        return(str(record.seq))

