import argparse
from collections import defaultdict,Counter
import csv
import json

def processKmers(kmerLength, cdr, kmerDict, key, cdrDict):
    totalLength = len(cdr)
    getValid = cdr[3:(totalLength-3)] ## discard first3 and last3 aa
    for i in range(0, len(getValid), 1):
        kmer = getValid[i:i+kmerLength]
        if len(kmer) == kmerLength:
            kmerDict[kmer][key] += 1
            cdrDict[kmer].add(cdr)
        ## make discontinuous 
        if len(kmer) == kmerLength:
            kmerDisc = kmer[0]+'.'+kmer[2:]
            kmerDict[kmerDisc][key] += 1
            cdrDict[kmerDisc].add(cdr)

def processSseq(sseqFile):
    #sseqFile = 'data/pnasRaw.csv'
    with open(sseqFile) as f:
        reader = csv.reader(f)
        kmerDictBeta = defaultdict(lambda:defaultdict(int))
        kmerDictAlpha = defaultdict(lambda:defaultdict(int))
        cdrDictBeta = defaultdict(set)
        cdrDictAlpha = defaultdict(set)
        for n, row in enumerate(reader):
            if n > 0:
                if row[0] and row[1] and row[2]:
                    key = ','.join(row[:3])
                    if row[3]:
                        cdrb = row[3]
                        for kmerL in range(3, 6):
                            processKmers(kmerLength=kmerL, cdr=cdrb, kmerDict=kmerDictBeta, key=key, cdrDict=cdrDictBeta)
                    if row[4]:
                        cdra = row[4]
                        for kmerL in range(3, 6):
                            processKmers(kmerLength=kmerL, cdr=cdra, kmerDict=kmerDictAlpha, key=key, cdrDict=cdrDictAlpha)
                    if row[5]:
                        cdra = row[5]
                        for kmerL in range(3, 6):
                            processKmers(kmerLength=kmerL, cdr=cdra, kmerDict=kmerDictAlpha, key=key, cdrDict=cdrDictAlpha)
    print("PROCESSED {}".format(sseqFile))
    return [kmerDictBeta, kmerDictAlpha, cdrDictAlpha, cdrDictBeta]


def writeDicts(kmerDict, chain):
    outFile = open('outputs/KMERS_{}.csv'.format(chain), 'w')
    writer = csv.writer(outFile)
    header= ['KMER','MER','KEY','COUNT']
    writer.writerow(header)
    for k, v in kmerDict.items():
        for k2, v2 in v.items():
            row = [k, len(k), k2, v2]
            writer.writerow(row)

def jsonDump(dict, jsonFile):
	with open(jsonFile, 'w') as F:
		json.dump(dict, F)

def main():
    parser = argparse.ArgumentParser(description='Script to make KMERS from single cell TCR sequencing')
    parser.add_argument('-S', help='csv file with identifiers and cdr3b, cdr3a and alt cdr3a, see example in data folder', required=True)
    args=parser.parse_args()
    sseqFile=args.S
    processedDicts=processSseq(sseqFile = sseqFile)
    writeDicts(kmerDict=processedDicts[0], chain='BETA') 
    writeDicts(kmerDict=processedDicts[1], chain='ALPHA')
    alphaDump = {i:','.join(list(j)) for i, j in processedDicts[2].items()}
    betaDump = {i:','.join(list(j)) for i, j in processedDicts[3].items()}
    jsonDump(dict=alphaDump, jsonFile='outputs/ALPHA_DUMP.json')
    jsonDump(dict=betaDump, jsonFile='outputs/BETA_DUMP.json')

if __name__ == "__main__":main()