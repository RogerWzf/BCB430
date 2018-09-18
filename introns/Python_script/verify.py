# Purpose: Verify the putative intron sequences.
# Authors: Fan Shen, Zhifan Wu
# Date: 2018-03-25


# Import Biopython and os.path
from Bio import SearchIO
import os.path


def verifyIntrons(blast_file, result_file, output_file):
    """ (str, str, str) -> NoneType
    
    Write the introns in output_file by analyze hits in blast_file and get the
    unmatched contigs/sequences from result_file.
    """
    # Get the Query Result
    blast_qresults = SearchIO.parse(blast_file, 'blast-xml')

    hit_query = {}

    for blast_qresult in blast_qresults:
        for hit in blast_qresult:
            for hsp in hit:
                if hit.query_id in hit_query:
                    hit_query[hit.query_id].append(hsp.query_range)
                else:
                    hit_query[hit.query_id] = [hsp.query_range]

    # Get the not match sequence.
    get_not_match(hit_query, result_file, output_file)


def get_not_match(contigs_dic, seq_file, output_file):
    """ (dic{str:list}, str, str) -> nonetype
    
    Use the information in contigs_dic to get the not aligned sequences from
    the list (no matter there is a hit or not) that is made by original file 
    seq_file.
    """
    # Open the file that contains the contigs
    file_handle = open(seq_file)

    # Get the unmatched seq list
    unmatched_list = []

    geneListFile = open(output_file, "w")
    if output_file.find("not") == -1:
        output = "|verify\n"
    else:
        output = "|not_introns\n"

    flag = False
    matched = False
    for line in file_handle.readlines():
        if line.find('>') != -1:
            name = line[1:-1]
            if not (name in contigs_dic):
                unmatched_list.append(name)
                geneListFile.write("%s\n" % line[:-1])
                flag = True
            else:
                unmatched_list.append(name)
                matches = contigs_dic[name]
                orgnized_matches = organize_match(matches)
                geneListFile.write("%s" % line[:-1] + output)
                matched = True
        else:
            if flag:
                geneListFile.write("%s\n" % line[:-1])
                flag = False
            elif matched:
                str_introns = ""
                if orgnized_matches[0][0] != 0:
                    str_introns = str_introns + line[:orgnized_matches[0][0] - 1]
                for i in range(len(orgnized_matches)):
                    start = orgnized_matches[i][1]
                    if i == len(orgnized_matches) - 1:
                        if start != len(line):
                            str_introns = str_introns + (line[start:])
                    else:
                        end = orgnized_matches[i + 1][0] - 1
                        if start != end:
                            str_introns = str_introns + (line[start:end])
                geneListFile.write("%s" % str_introns)
                matched = False


def organize_match(matches):
    """ List of tuples -> List of tuples
    
    Organize the match tuples:
    1. Delete duplicates
    2. combine the matches to bigger one
    """
    if len(matches) == 1:
        return matches
    else:
        sorted_matches = sorted(matches, key=lambda tup: tup[0])
        new_l = []
        compared = sorted_matches[0]

        for i in range(1, len(sorted_matches)):
            if (compared > sorted_matches[i]) - (compared < sorted_matches[i]):  # # cmp(a, b) == (a > b) - (a < b)

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
    verifyIntrons("./introns_aga_cds.xml", "./result.fasta", "./verify.fasta")

    if os.path.isfile("RRHintrons_aga_putative.xml"):
        verifyIntrons("RRHintrons_aga_putative.xml",
                      "RRH_putative_introns.fasta",
                      "./verify_RRH_not_putative_introns.fasta")
    if os.path.isfile("RRHintrons_aga_verify.xml"):
        verifyIntrons("RRHintrons_aga_verify.xml", "RRH_verify_introns.fasta",
                      "./verify_RRH_not_verify_introns.fasts")