#! /usr/bin/env python
# Author : Stefan Engelen - Genoscope - LBGB
# Contact : Stefan.ENGELEN@genoscope.cns.fr
# python/2.7.8
# pysam/0.8.4

import os
import argparse
import sys
import re
import pysam
from collections import Counter

# http://pysam.readthedocs.org/en/latest/api.html#pysam.AlignedRead.cigar
MATCH  = 0  # M
INS    = 1  # I
DEL    = 2  # D
SKIP   = 3  # N
SOFT   = 4  # S
HARD   = 5  # H
PAD    = 6  # P
EQUAL  = 7  # =
DIFF   = 9  # X


param="best"
outPrefix="test"
bamfile = ""
hitstat = False
readtype = "-"
orientation= False


parser = argparse.ArgumentParser(description='Requirement: module load python/3.7 pysam/0.15.4')

parser.add_argument( '-b', '--bam',          required=True,  type=str, metavar='BAM',    help="Bam file to process" )
parser.add_argument( '-o', '--outPrefix',   default='test', type=str, metavar='PREFIX', help='Prefix of the output file')
parser.add_argument( '-t', '--readtype',    default='-',    type=str,                   help='[nanopore_type | string]' )
parser.add_argument( '-s', '--statType',    default='best', type=str, choices=['default','best','allbest','all'] )
parser.add_argument( '-a', '--orientation', action='store_true',                        help='To print orientation of read mapping on reference' )
parser.add_argument( '-p', '--hitstat',     action='store_true',                        help='To calculate statistics per hit. Need \'allbest\' or \'all\' parameter to work' )
args = parser.parse_args()

def cigar_profile(cigar_tuples):
        """
        Return a dictionary that tabulates the total number
        of bases associated with each CIGAR operation.

        cigar_tuples is a list of (op, length) tuples.
        """
        cigar_prof = Counter()
        for cigar_tuple in cigar_tuples:
                cigar_prof[cigar_tuple[0]] += cigar_tuple[1]
        return cigar_prof

def get_total_unaligned(cigar_prof):
        """
        return the total number unaligned bases (hard or softclips.)
        """
        return cigar_prof[HARD] + cigar_prof[SOFT]

def get_query_alignment(read):
        """
        Return the query alignment start and end.
        """
        start = 0
        end = 0
        cigar_prof = cigar_profile(read.cigar)
        read_length = cigar_prof[INS] + cigar_prof[MATCH] + get_total_unaligned(cigar_prof)
        cigar = read.cigarstring
        regex = re.compile("([0-9]+)[H|S]")
        r = regex.findall(str(cigar))
        if r != []:
                if len(r)<2:
                        if re.search(re.compile("^([0-9]+)[H|S]"), cigar) is None:
                                r.insert(0,0)
                        elif re.search(re.compile("([0-9]+)[H|S]$"), cigar) is None:
                                r.insert(1,0)
        else :
                r.insert(0,0)
                r.insert(1,0)
        start = int(r[0])
        end = read_length - int(r[1])
        return start,end

# Pass 1:
# iterate through each BAM alignment
# and store the best alignment foR the read


outFileStat = "%s.stat" % (args.outPrefix)
outFileStatPerHit = "%s_perHit.stat" % (args.outPrefix)
outFileSam = "%s.bam" % (args.outPrefix)

inFileBamHandler = pysam.AlignmentFile(args.bam, "rb")
outFileStatHandler = open(outFileStat,'w')

if args.statType == "allbest" or args.statType == "best" or args.statType == "all":

        outFileSamHandler = pysam.AlignmentFile(outFileSam, "wb", template=inFileBamHandler)
        # best_align:
        #       key is query name
        #   val is tuple of (alignment length, cigar_prof)
        read_align = {}
        n=0

        for read in inFileBamHandler :
                n=n+1

                # flag 4 == read unmapped
                if not read.flag & 4:
                        cigar_prof = cigar_profile(read.cigar)

                        if read.qname not in read_align:
                                read_align[read.qname] = [read]
                        elif read.qname in read_align:
                                read_align[read.qname].append(read)

        # Pass 2:
        # Report the alignment and error profile for each read's best alignment

        if args.orientation == True and args.statType != "allbest" :
                outFileStatHandler.write('\t'.join(['query', 'chr_name', 'read_type', 'reference_start', 'reference_end', 'read_start', 'read_end', 'read_len', 'align_len', 'unalign_len', 'matches', 'mismatches', 'insertions', 'deleti
ons_', 'tot_errors','hit_nb','orientation']))
                outFileStatHandler.write('\n')
        else:
                outFileStatHandler.write('\t'.join(['query', 'chr_name', 'read_type', 'reference_start', 'reference_end', 'read_start', 'read_end', 'read_len', 'align_len', 'unalign_len', 'matches', 'mismatches', 'insertions', 'deleti
ons_', 'tot_errors','hit_nb']))
                outFileStatHandler.write('\n')

        if (args.hitstat == True):
                outFileStatPerHitHandler = open(outFileStatPerHit,'w')
                if args.orientation == False :
                        outFileStatPerHitHandler.write('\t'.join(['query', 'chr_name', 'read_type', 'reference_start', 'reference_end', 'read_start', 'read_end', 'read_len', 'align_len', 'unalign_len', 'matches', 'mismatches', 'insert
ions', 'deletions_', 'tot_errors','score','hit_nb']))
                        outFileStatPerHitHandler.write('\n')
                else :
                        outFileStatPerHitHandler.write('\t'.join(['query', 'chr_name', 'read_type', 'reference_start', 'reference_end', 'read_start', 'read_end', 'read_len', 'align_len', 'unalign_len', 'matches', 'mismatches', 'insert
ions', 'deletions_', 'tot_errors','score','hit_nb','orientation']))
                        outFileStatPerHitHandler.write('\n')

        for read_name in read_align: # on recupere les cles (read_name)
                list_read_filter = []
                carac_list_sorted = []
                pos_init=0

                for read in read_align[read_name]: # chaque hit du read read_name
#                       cigar_prof = cigar_profile(read.cigar)
#                       read_length = read.query_length + get_total_unaligned(cigar_profile(read.cigar))
                        score = read.opt('AS')
#                       score = cigar_prof[MATCH] - read.opt('NM') - read.opt('NM') + cigar_prof[INS] + cigar_prof[DEL]
                        start,end = get_query_alignment(read)
                        pos = 0
                        for carac in carac_list_sorted:
                                if score > carac[4]:
                                    break
                                pos=pos+1
                        carac = (read.reference_start,read.reference_end,start,end,score,pos_init,read.reference_name)
                        pos_init = pos_init + 1
                        carac_list_sorted.insert(pos,carac)

                carac_list_filter = [carac_list_sorted[0]]
                list_read_filter.append(read_align[read_name][carac_list_sorted[0][5]])

                if args.statType == "all":
                        for i in range(1,len(carac_list_sorted)):
                                list_read_filter.append(read_align[read_name][carac_list_sorted[i][5]])
                                carac_list_filter.append(carac_list_sorted[i])

                if args.statType == "allbest":
                        for i in range(1,len(carac_list_sorted)):
                                j=0
                                while j <len(carac_list_filter):
                                        if ((carac_list_filter[j][6] == carac_list_sorted[i][6]) and ((carac_list_filter[j][0] <= carac_list_sorted[i][0] <=  carac_list_filter[j][1]) or (carac_list_filter[j][0] <= carac_list_sorted[i]
[1] <=  carac_list_filter[j][1]))) \
                                        or ((carac_list_filter[j][2] <= carac_list_sorted[i][2] <=  carac_list_filter[j][3]) or (carac_list_filter[j][2] <= carac_list_sorted[i][3] <=  carac_list_filter[j][3])) \
                                        or ((carac_list_filter[j][6] == carac_list_sorted[i][6]) and ((carac_list_sorted[i][0] <= carac_list_filter[j][0] <=  carac_list_sorted[i][1]) or (carac_list_sorted[i][0] <= carac_list_filter[j]
[1] <=  carac_list_sorted[i][1]))) \
                                        or ((carac_list_sorted[i][2] <= carac_list_filter[j][2] <=  carac_list_sorted[i][3]) or (carac_list_sorted[i][2] <= carac_list_filter[j][3] <=  carac_list_sorted[i][3])) :
                                                break
                                        j=j+1
                                if j == len(carac_list_filter):
                                        list_read_filter.append(read_align[read_name][carac_list_sorted[i][5]])
                                        carac_list_filter.append(carac_list_sorted[i])

                name = list_read_filter[0].qname
                if args.readtype == "nanopore_type" :
                        regex = re.compile(".*_(pass|fail)_.*")
                        r = regex.search(name)
                        if r : read_type = r.group(1)
                        else : read_type = "-"

                else :
                    read_type = "-"

#               read_length = list_read_filter[0].query_length + get_total_unaligned(cigar_prof)
                read = list_read_filter[0]
                cigar_prof = cigar_profile(read.cigar)
                read_length = cigar_prof[INS] + cigar_prof[MATCH] + get_total_unaligned(cigar_prof)
                ref_name_tot = read.reference_name
                ref_name_hash = {}

                match_tot = 0
                hit_nb = 0
                mismatch_tot = 0
                total_errors_tot = 0
                read_align_len_tot = 0
                insertion_tot = 0
                deletion_tot = 0

                sens = "NA"

                for read in list_read_filter:
                        outFileSamHandler.write(read)
                        cigar_prof = cigar_profile(read.cigar)
                        insertion = cigar_prof[INS]
                        deletion = cigar_prof[DEL]
                        mismatch = read.opt('NM') - insertion - deletion
                        match = cigar_prof[MATCH] - mismatch
                        align_length = insertion + cigar_prof[MATCH]
                        unaligned_length = get_total_unaligned(cigar_prof)
                        total_errors_tot = total_errors_tot + read.opt('NM')
                        hit_nb = hit_nb + 1
                        match_tot = match_tot + match
                        mismatch_tot = mismatch_tot + mismatch
                        read_align_len_tot = read_align_len_tot + align_length
                        insertion_tot = insertion_tot + insertion
                        deletion_tot = deletion_tot + deletion
                        ref_name = inFileBamHandler.getrname(read.tid)
                        if read.is_reverse :
                                sens = "-"
                        else :
                                sens = "+"
                        if ref_name in ref_name_hash :
                                ref_name_hash[ref_name] = ref_name_hash[ref_name] + 1
                        else :
                                ref_name_hash[ref_name] = 1
                        start,end = get_query_alignment(read)

                        if args.statType != "best" and args.hitstat == True:
                                if args.orientation == False :
                                        outFileStatPerHitHandler.write('\t'.join(str(s) for s in [name, ref_name, read_type, read.reference_start, read.reference_end, start, end,read_length,align_length,unaligned_length,match,mismatch
,insertion,deletion,read.opt('NM'),read.opt('AS'),hit_nb]))
                                        outFileStatPerHitHandler.write('\n')
                                else :
                                        outFileStatPerHitHandler.write('\t'.join(str(s) for s in [name, ref_name, read_type, read.reference_start, read.reference_end, start, end,read_length,align_length,unaligned_length,match,mismatch
,insertion,deletion,read.opt('NM'),read.opt('AS'),hit_nb,sens]))
                                        outFileStatPerHitHandler.write('\n')
                        elif args.hitstat == True and args.statType == "best":
                                print( "\n*** \'histat\' option works only with \'all\' and \'allbest\' parameters. *** \n" )
                                parser.print_help()
                                sys.exit()

                l = len(ref_name_hash)
                if  l > 1:
                        ref_name_tot = ref_name_tot + ";+" + str(l-1)
                unaligned_len_tot = read_length - read_align_len_tot
                if args.orientation == True and args.statType != "allbest" :
                        outFileStatHandler.write('\t'.join(str(s) for s in [name, ref_name_tot, read_type, read.reference_start,read.reference_end, start, end, read_length,read_align_len_tot,unaligned_len_tot,match_tot,mismatch_tot,in
sertion_tot,deletion_tot,total_errors_tot,hit_nb,sens]))
                        outFileStatHandler.write('\n')
                else :
                        outFileStatHandler.write('\t'.join(str(s) for s in [name, ref_name_tot, read_type, read.reference_start,read.reference_end, start, end, read_length,read_align_len_tot,unaligned_len_tot,match_tot,mismatch_tot,in
sertion_tot,deletion_tot,total_errors_tot,hit_nb]))
                        outFileStatHandler.write('\n')

        outFileSamHandler.close()
outFileStatHandler.close()
if (args.hitstat == True): outFileStatPerHitHandler.close()
inFileBamHandler.close()
