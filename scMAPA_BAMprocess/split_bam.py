import pysam
import csv
import sys

cluster_file = sys.argv[1];
bam_file = sys.argv[2];

print("Cluster file: "+cluster_file);
print("Bam file: "+bam_file)

cluster_dict = {}
with open(cluster_file) as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    #skip header
    header = next(csv_reader)
    for row in csv_reader:
        #print(row[0])
        #print(row[1])
        cluster_dict[row[0]] = row[1]

clusters = set(x for x in cluster_dict.values())

fin = pysam.AlignmentFile(bam_file, "rb")

fouts_dict = {}
for cluster in clusters:
    fout = pysam.AlignmentFile("cluster" + cluster + ".bam", "wb", template = fin)
    fouts_dict[cluster] = fout

for read in fin:
    tags = read.tags
    CB_list = [ x for x in tags if x[0] == "CB"]
    if CB_list:
        cell_barcode = CB_list[0][1]
    # the bam files may contain reads not in the final clustered barcodes
    # will be None if the barcode is not in the clusters.csv file
    else: 
        continue
    cluster_id = cluster_dict.get(cell_barcode)
    if cluster_id:
        fouts_dict[cluster_id].write(read)
## do not forget to close the files
fin.close()
for fout in fouts_dict.values():
    fout.close()