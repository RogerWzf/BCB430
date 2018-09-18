# Purpose: Find the specific intron sequences for given protein name.
# Authors: Fan Shen, Zhifan Wu
# Date: 2018-03-25

def find_introns_for_protein(protein_name, result_file, write_file):
    """ (str, str, stre) --> nonetype
    
    Write in to write_file the name and sequences, where find protein_name
    in result_file.
    """
    file_handle = open(result_file)

    geneListFile = open(write_file, "w")
    str_len = 0

    flag = False
    for line in file_handle.readlines():
        if line.find(protein_name) != -1:
            geneListFile.write("%s\n" % line[:-1])
            flag = True
        else:
            if flag:
                geneListFile.write("%s\n" % line[:-1])
                str_len = str_len + len(line[:-1])
            flag = False

    geneListFile.write("The residues length is %s\n" % str_len)


if __name__ == '__main__':
    find_introns_for_protein("RRH", "result.fasta", "./RRH_putative_introns.fasta")
    find_introns_for_protein("RRH", "verify.fasta", "./RRH_verify_introns.fasta")
