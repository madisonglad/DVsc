import subprocess
import optparse
import os
import math

class lf_opt:
    def __init__(self):

        parser = optparse.OptionParser()
        parser.add_option("-m", "--method", dest="method", help="RNA-Seq method")
        parser.add_option("-i", "--input", dest="input", help="input of viruscountfile")
        parser.add_option("-o", "--output", dest="output", help="output of viruscountfilterfile")

        self.options, self.args = parser.parse_args()
        print(self.options, self.args)

def totalreadscount(viruscountfile):
    n = 0
    num = 0
    with open(viruscountfile) as f:
        for line in f:
            num = num + int(line.split('\t')[-1])
            n = n + 1
    return n, num


def qc(viruscountfile, n, totalnum, method, viruscountfilterfile):
    outf = open(viruscountfilterfile, 'w')
    with open(viruscountfile) as f:
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
                if totalnum >= 2981 or content >= 0.6:
                    filter = 'Y'
                        
            if filter == 'N':
                    continue
            outf.write(line.strip() + '\t' + str(content) + str(n) + '\n')                   

def main():
    opt = lf_opt()
    method = opt.options.method
    viruscountfile = opt.options.input
    viruscountfilterfile = opt.options.output

    n, totalnum = totalreadscount(viruscountfile)
    qc(viruscountfile, n, totalnum, method, viruscountfilterfile)    
    
if __name__ == "__main__":
    main()

    
