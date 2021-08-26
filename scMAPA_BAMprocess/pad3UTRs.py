#1. read in a bam file
#2. read in the ref. file
#3. For each 3'UTR, mask the bam file with the the region 
#4. For each masked region, retrieve a read and pad them 
import pybedtools
import os
import argparse
parser = argparse.ArgumentParser(description='Description of this program')
parser.add_argument('-r','--refbed', help='bed file including only 3UTR region of each gene', required=True)
parser.add_argument('-i','--inbam', help='input bam file of the reads to be padded', required=True)
args = vars(parser.parse_args())
refbed = args['refbed']
inbam=args['inbam']

#threeUTRs=pybedtools.example_bedtool('/ihome/hpark/hyp15/storage/scMAPA/mm10_m18_refseq_extracted_3UTR.bed')
#threeUTRs=pybedtools.example_bedtool('/ihome/hpark/hyp15/storage/scMAPA/mm10_test1.bed')
bams4threeUTRs=pybedtools.BedTool(inbam).intersect(refbed,s=True,bed=True,wb=True)
'''outbedOrig=open('orig'+refbed+'.bed', "w")
for bam4threeUTR in bams4threeUTRs:
    outbedOrig.write("\t".join([bam4threeUTR['chrom'], str(bam4threeUTR['start']), str(bam4threeUTR['end']), bam4threeUTR['name']+"_"+bam4threeUTR['fields'][15], bam4threeUTR['score'], bam4threeUTR['strand']])+"\n")
outbedOrig.close()'''

outbedTrans=open(inbam+'.trans.bed', "w")
for bam4threeUTR in bams4threeUTRs:
    if bam4threeUTR['strand']=="+":
        startPos=bam4threeUTR['start']
        if startPos > bam4threeUTR['fields'][13]:
            startPos=bam4threeUTR['fields'][13]#is the start coordinate in the reference
        outbedTrans.write("\t".join([bam4threeUTR['chrom'], str(startPos), str(bam4threeUTR['end']), bam4threeUTR['name']+"_"+bam4threeUTR['fields'][15], bam4threeUTR['score'], bam4threeUTR['strand']])+"\n")
    else:#minus strand
        endPos=bam4threeUTR['end']
        if endPos < bam4threeUTR['fields'][14]:
            endPos=bam4threeUTR['fields'][14]#is the end coordinate in the reference
        outbedTrans.write("\t".join([bam4threeUTR['chrom'], str(bam4threeUTR['start']), str(endPos), bam4threeUTR['name']+"_"+bam4threeUTR['fields'][15], bam4threeUTR['score'], bam4threeUTR['strand']])+"\n")
outbedTrans.close()
