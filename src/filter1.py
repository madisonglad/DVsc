import subprocess
import optparse
import os
import math

class lf_opt:
    def __init__(self):

        parser = optparse.OptionParser()
        parser.add_option("-s", "--sample", dest="sampleid", help="sample id")
        parser.add_option("-f", "--annotatefile", dest="annotatefile", help="path of annotatefile")
        parser.add_option("-r", "--readlength", dest="readlength", help="reads length")
        parser.add_option("-n", "--nclistoutfile", dest="nclistoutfile", help="outfile of nclist")
        parser.add_option("-c", "--nclistfilteroutfile", dest="nclistfilteroutfile", help="outfile of nclistfilter")
        parser.add_option("-a", "--annotatedfilteroutfile", dest="annotatedfilteroutfile", help="outfile of annotatedfilter")


        self.options, self.args = parser.parse_args()
        print(self.options, self.args)

def anno_dict(file):
    dict_length = {}
    dict_name = {}
    with open(file) as f:
        for line in f:
            ncid = line.split('\t')[0]
            length = line.split('\t')[1]
            name = line.split('\t')[2]            
            dict_length[ncid] = length
            dict_name[ncid] = name
    return dict_length, dict_name

def filter(nclistoutfile, nclistfilteroutfile, readlength):
    outf = open(nclistfilteroutfile, 'w')
    with open(nclistoutfile) as f:
        for line in f:
            filter = 'Y'
            qualifyreads=line.split('\t')[1]
            cov=line.split('\t')[4]
            maxstr=line.split('\t')[5]
            
            if float(qualifyreads) < 10:
                filter = 'N'
               
            if float(cov) < 0.011554:
                filter = 'N'
            
            if float(readlength) <75:
                if float(maxstr) < 50:
                    filter = 'N'
                    
            if float(readlength) > 75:
                if float(maxstr) < 165.925:
                    filter = 'N'
                        
            if filter == 'N':
                    continue
            outf.write(line)                   

def anno(sampleid, annotatedfilteroutfile, dict):
    outf = open(annotatedfilteroutfile, 'w')
    for root, dirs, files in os.walk(sampleid):
        for file in files:
            if file.endswith('nclist.filter.txt'):
                filepath = root + '/' + file
                with open(filepath) as f:
                    for line in f:
                        ncid = 'NC_' + line.split('\t')[0].split('_')[1]
                        outf.write(dict[ncid] + '\t' + line)


def main():
    opt = lf_opt()
    sampleid = opt.options.sampleid
    annotatefile = opt.options.annotatefile
    readlength = opt.options.readlength
    
    nclistoutfile = opt.options.nclistoutfile
    nclistfilteroutfile = opt.options.nclistfilteroutfile
    annotatedfilteroutfile = opt.options.annotatedfilteroutfile

    dict_length, dict_name = anno_dict(annotatefile)
    filter(nclistoutfile, nclistfilteroutfile,readlength)    
    anno(sampleid, annotatedfilteroutfile, dict_name)
    
if __name__ == "__main__":
    main()

    
