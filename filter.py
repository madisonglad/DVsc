import subprocess
import optparse
import os
import math

class lf_opt:
    def __init__(self):

        parser = optparse.OptionParser()
        parser.add_option("-s", "--sample", dest="sampleid", help="sample id")
        parser.add_option("-m", "--method", dest="method", help="RNA-Seq method")
        parser.add_option("-o", "--outpath", dest="outpath", help="path of output")

        self.options, self.args = parser.parse_args()
        print(self.options, self.args)

def totalreadscount(id):
    n = 0
    num = 0
    with open(id +  '/virus_count.txt') as f:
        for line in f:
            num = num + int(line.split('\t')[-1])
            n = n + 1
    return n, num


def qc(sampleid, n, totalnum, method, outpath):
    outf = open(outpath + '/' + sampleid + '_virus_count_filter.txt', 'w')
    with open(sampleid + '/virus_count.txt') as f:
        for line in f:
            filter = 'N'
            readnums = int(line.split('\t')[2])
            
            content = readnums / totalnum
            
            if method == 'scRNAseq':
                if content < 1 and readnums < 50:
                    filter = 'N'
                else:
                    filter = 'Y'
                    
            if method == 'bulkRNAseq':
                if totalnum > 100000:
                    if content > 0.04:
                        filter = 'Y'
                if totalnum > 3000 and totalnum <= 100000:
                    if content > 0.20:
                        filter = 'Y'
                if totalnum <= 3000:
                    if content > 0.70:
                        filter = 'Y'
                        
            if filter == 'N':
                    continue
            outf.write(line.strip() + '\t' + str(content) + str(n) + '\n')                   

def main():
    opt = lf_opt()
    sampleid = opt.options.sampleid
    method = opt.options.method
    outpath = opt.options.outpath

    n, totalnum = totalreadscount(sampleid)
    qc(sampleid, n, totalnum, method, outpath)    
    
if __name__ == "__main__":
    main()

    
