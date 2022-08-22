import subprocess
import optparse
import os
import pandas as pd
from collections import defaultdict


class lf_opt:
    def __init__(self):
        parser = optparse.OptionParser()
        parser.add_option("-s", "--sample", dest="sampleid", help="sample id")
        parser.add_option("-a", "--annotatedfilterfile", dest="annotatedfilterfile", help="annotatedfilterfile")
        parser.add_option("-v", "--viruscountfilterfile", dest="viruscountfilterfile", help="viruscountfilterfile")
        parser.add_option("-o", "--output", dest="output", help="output of assemble file")

        self.options, self.args = parser.parse_args()
        print(self.options, self.args)
        
def runcmd(command):
    ret = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=10000)
    if ret.returncode == 0:
        pass
    else:
        print(ret)
        
def virus_dict(annotatedfilterfile):
    dict = defaultdict(list)
    with open(annotatedfilterfile) as f:
        for line in f:
            virusname = line.split('\t')[0].replace(' ', '_').replace('(', '_').replace(')', '_').replace('/','_')
            ncid = line.split('\t')[1]
            dict[virusname].append(ncid)
    return dict
   
def assemble(sample, viruscountfilterfile, dict_virus, assemblefile):
    with open(viruscountfilterfile) as f:
        for line in f:
            #virusname = line.split('\t')[1].replace('_', ' ')
            virusname = line.split('\t')[1]
            ncids = dict_virus[virusname]
            for i in ncids:
                cmd1 = 'stringtie -o ' + sample + '_' + i + '.single.gtf ' + sample + '/viral_bam/' + i + '.bam'
                #print(cmd1)
                runcmd(cmd1)
                
    cmd2 = 'cat ' + sample + "*.single.gtf | sed -r '/^#/d' > " + assemblefile
    #print(cmd2)
    runcmd(cmd2)
    
    cmd3= 'rm *.single.gtf'            
    runcmd(cmd3)
            
def main():
    opt = lf_opt()
    sampleid = opt.options.sampleid
    annotatedfilterfile = opt.options.annotatedfilterfile
    viruscountfilterfile = opt.options.viruscountfilterfile
    assemblefile = opt.options.output

    dict_virus = virus_dict(annotatedfilterfile)

    assemble(sampleid, viruscountfilterfile, dict_virus, assemblefile)
    

    
if __name__ == "__main__":
    main()
