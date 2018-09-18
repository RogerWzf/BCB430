# Purpose: Get the putative(intron-like) sequence from the BLAST results.
# Authors: Fan Shen, Zhifan Wu
# Date: 2018-03-25


# Import Biopython
from Bio import SearchIO


def getIntrons(blast_file, contigs_file, output_file):
    """ (str, str, str) -> NoneType
    
    Write the introns in output_file by analyze hits in blast_file and get the
    unmatched contigs/sequences from contigs_file.
    """
    # Parse the blast result.
    blast_qresults = SearchIO.parse(blast_file, 'blast-xml')

    # Create an empty dictionary
    query_dic = {}

    # Loop over the blast result to find all the hit query with hit's name and range
    for blast_qresult in blast_qresults:
        for hit in blast_qresult:
            for hsp in hit:
                if hit.query_id in query_dic:
                    query_dic[hit.query_id][1].append(hsp.query_range)
                else:
                    query_dic[hit.query_id] = (hit.id, [hsp.query_range])

    # Get the not match part from original file	    
    get_not_match(query_dic, contigs_file, output_file)


# Helper functions #
def get_not_match(query_dic, contigs_file, output_file):
    """ (dict{str: list}, str, str) -> NoneType
    
    Use the information in query_dic to get the not aligned sequences from the
    contig dictionary that is made by original file contigs_file. And then write the
    result to output_file.
    """
    # Open the file that contains the contigs
    file_handle = open(contigs_file)
    contigs_dic = read_contigs_data(file_handle)

    # Initialize the Intron list
    list_introns = []
    geneListFile = open(output_file, "w")

    # Initialize two variable for the analysis purpose.
    # Please note here please comment these two variable to firstly create
    # result.fasta fro further use, and then uncomment these two variables
    # to create result1.fasta for analysis purpose.
    # Number of hits.
    num_hit_n = 0
    # Number of nucleotides.
    num_all_n = 0

    for contig_name in query_dic:
        output = "|" + query_dic[contig_name][0] + "|intron"

        matches = query_dic[contig_name][1]
        orgnized_matches = organize_match(matches)
        num_intron = 1
        num_all_n = num_all_n + len(contigs_dic[contig_name])

        if orgnized_matches[0][0] != 0:
            name = ">" + contig_name + output + str(num_intron)
            sequence = contigs_dic[contig_name][:orgnized_matches[0][0] - 1]

            list_introns.append(name)

            geneListFile.write("%s\n" % name)
            geneListFile.write("%s\n" % sequence)

            num_hit_n = num_hit_n + len(sequence)
            num_intron = num_intron + 1

        for i in range(len(orgnized_matches)):
            start = orgnized_matches[i][1]
            if i == len(orgnized_matches) - 1:
                if start != len(contigs_dic[contig_name]):
                    name = ">" + contig_name + output + str(num_intron)
                    sequence = contigs_dic[contig_name][start:]

                    list_introns.append(name)

                    geneListFile.write("%s\n" % name)
                    geneListFile.write("%s\n" % sequence)

                    num_hit_n = num_hit_n + len(sequence)
                    num_intron = num_intron + 1
            else:
                end = orgnized_matches[i + 1][0] - 1
                if start != end:
                    name = ">" + contig_name + output + str(num_intron)
                    sequence = contigs_dic[contig_name][start:end]

                    list_introns.append(name)

                    geneListFile.write("%s\n" % name)
                    geneListFile.write("%s\n" % sequence)

                    num_hit_n = num_hit_n + len(sequence)
                    num_intron = num_intron + 1

    in_length = len(list_introns)
    geneListFile.write("The introns list is %d\n" % in_length)
    all_length = len(contigs_dic)
    geneListFile.write("All number of seq in contigs is %d\n" % all_length)
    hit_length = len(query_dic)
    geneListFile.write("All number of hit in contigs is %d\n" % hit_length)

    geneListFile.write("Num of hit nucleotides is %d\n" % num_hit_n)
    geneListFile.write("Num of all n is %d\n" % num_all_n)


# Helper Function of get_not_match function#
def read_contigs_data(file_handle):
    """ (file open for reading) -> dict
    
    Read from the file_handle and put the data into a dictionary that have name
    of contigs as key and sequence as value.
    """
    contigs_dic = {}
    name = ''

    for line in file_handle.readlines():
        if line.find('>') != -1:
            name = simple_name(line)
        else:
            # put the contigs seq to the dic
            if name in contigs_dic:
                contigs_dic[name] += line[:-1]
            else:
                contigs_dic[name] = line[:-1]
    return contigs_dic


# Helper function for read_contigs_data#
def simple_name(name):
    """ (str) -> str
    
    Simple the name of the contigs that match the name in hit
    """
    index = name.find(' ')
    return name[1:index]


def organize_match(matches):
    """ List of tuples -> List of tuples
    
    Organize the match tuples:
    1. Delete duplicates
    2. combine the matches to bigger one
    """
    # If there are only one matches
    if len(matches) == 1:
        return matches
    # If there are more than one matches
    else:
        sorted_matches = sorted(matches, key=lambda tup: tup[0])
        new_l = []
        compared = sorted_matches[0]

        for i in range(1, len(sorted_matches)):
            if (compared > sorted_matches[i]) - (compared < sorted_matches[i]):  # cmp(a, b) == (a > b) - (a < b)
                end0 = compared[1]
                end1 = sorted_matches[i][1]
                if sorted_matches[i][0] <= end0 + 1:
                    compared = (compared[0], max(end0, end1))
                else:
                    new_l.append(compared)
                    compared = sorted_matches[i]

        new_l.append(compared)

        return new_l


if __name__ == '__main__':
    getIntrons("contigs_against_exome.xml", "Anolis_RSO6ex3_5_Trinity_NR.fasta"
               , "./result1.fasta")
