import pysam, argparse, sys, gzip, os



parser=argparse.ArgumentParser()

parser.add_argument("-i", "--inBam", help="Input bam file", action = "store", required = True)
parser.add_argument("-o", "--outBam", help="Output bam file", action = "store", required = True)
parser.add_argument("-t", "--target_alleles", help="Target alleles tsv file (gzip not supported)", action = "store", required=True)

args = parser.parse_args()

#input bam reader
inBam = pysam.AlignmentFile(args.inBam, "rb")

#output bam writer
outBam = pysam.AlignmentFile(args.outBam + "_unsorted", "wb", template = inBam)

targetsFile = open(args.target_alleles, "r")

selectedReadNames=set()

for targetsLine in targetsFile:
    
    targetContig, targetPos, targetBase = targetsLine.split()
    targetPos = int(targetPos) - 1
    entries = inBam.fetch(contig=targetContig, start=targetPos, stop = targetPos+1)
    
    for entry in entries:
        try: queryPairedPositions,refPairedPositions = zip(*entry.get_aligned_pairs())
        except: continue
        if targetPos in refPairedPositions:
            readTargetPos = queryPairedPositions[refPairedPositions.index(targetPos)]
            if readTargetPos != None:
                readTargetBase = entry.query_sequence[readTargetPos]
                if readTargetBase.upper() == targetBase:
                    selectedReadNames.add(entry.query_name)


sys.stderr.write("\nFound {} entries carrying a target base\n".format(len(selectedReadNames)))

#index file for pulling out reads
sys.stderr.write("\nBuilding read name index for extracting selected pairs...\n")
readIndex = pysam.IndexedReads(inBam)
readIndex.build()

sys.stderr.write("\nWriting {} selected entries...\n".format(len(selectedReadNames)))
for readName in selectedReadNames:
    entries=readIndex.find(readName)
    for entry in entries: outBam.write(entry)

inBam.close()
outBam.close()
targetsFile.close()

sys.stderr.write("\nSorting output bam file...\n")

pysam.sort("-o", args.outBam, args.outBam + "_unsorted")

os.remove(args.outBam + "_unsorted")

sys.stderr.write("\nIndexing sorted bam file...\n")

pysam.index(args.outBam)

sys.stderr.write("\nDone.\n")
