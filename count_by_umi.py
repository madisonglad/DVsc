import subprocess
import optparse
import os
import pandas as pd
from collections import defaultdict


class lf_opt:
    def __init__(self):
        parser = optparse.OptionParser()
        parser.add_option("-s", "--sample", dest="sampleid", help="sample id")
        parser.add_option("-o", "--outpath", dest="outpath", help="path of output")

        self.options, self.args = parser.parse_args()
        print(self.options, self.args)

def runcmd(command):
    ret = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=10000)
    if ret.returncode == 0:
        pass
    else:
        print(ret)

def virus(id):
    with open(id + 'filter.anno.txt') as f:
        dict_virus = defaultdict(list)
        nc_list = []
        for line in f:
            virus_name = line.split('\t')[0].replace(' ', '_').replace('/', '_').replace('(', '_').replace(')', '_')
            virus_nc = line.split('\t')[1]
            dict_virus[virus_name].append(virus_nc)
            nc_list.append(virus_nc)
    return dict_virus, nc_list

def virus_filter(id,outpath):
    with open(outpath + '/' + id + '_virus_count_filter.txt') as f:
        dict_virus_filter = defaultdict(list)
        for line in f:
            virus_name = line.split('\t')[1]
            dict_virus_filter[virus_name] = ''
    return dict_virus_filter
    
def create_umi(id, nc_list):
    path = id + '/viral_bam'
    for root, dirs, files in os.walk(path):
        for file in files:
            if file.split('.')[0] in nc_list:
                outf = open(id + '/tmp/' + file.split('.')[0] + '.ncid.1', 'w')
                cmd = "samtools view " + id + '/viral_bam/' + file + " | awk '{print $1}' > " + id + '/tmp/' + file.split('.')[0] + '.ncid'

                # sub barcode and umi
                with open(id + '/tmp/' + file.split('.')[0] + '.ncid') as f:
                    for line in f:
                        line_new = line.split('_')[-2] + '_' + line.split('_')[-1]
                        outf.write(line_new)

def create_uniq_umi(id, dict_virus):
    for virus in dict_virus:
        virus = virus.replace(' ', '_').replace('/', '_').replace('(', '_').replace(')', '_')
        print(dict_virus[virus])
        for nc in dict_virus[virus]:
            cmd = 'cp ' + id + '/tmp/' + nc + '.ncid.1 ' + id + '/tmp2/' + virus + '/'
            runcmd(cmd)
        cmd = 'cat ' + id + '/tmp2/' + virus +'/*.ncid.1 > ' + id + '/tmp2/' + virus + '.merge.ncid'
        runcmd(cmd)
        cmd = 'sort ' + id + '/tmp2/'  + virus + '.merge.ncid | uniq > ' + id + '/single_count/' + virus + '.uniq.ncid'
        runcmd(cmd)

def count(id, dict_virus):
    for virus in dict_virus:
        outf = open(id + '/single_count/' + virus + '.tmp.count', 'w')
        dict  = defaultdict(list)
        with open(id + '/single_count/' + virus + '.uniq.ncid') as f:
            for line in f:
                barcode = line.split('_')[0]
                dict[barcode].append(line.split()[0])
        for barcode in dict:
            line_new = virus + '\t' + barcode + '\t' + str(len(dict[barcode])) + '\n'
            outf.write(line_new)

def filter_countbyumi(id, dict_virus_filter, outpath):
    outf = open(outpath + '/' + id + '_countbyumi_filter.txt','w')
    with open(id + '/countbyumi.txt') as f:
        for line in f:
            virus_name = line.split('\t')[0]
            if virus_name in dict_virus_filter:
                outf.write(line)

def main():
    opt = lf_opt()
    sampleid = opt.options.sampleid
    outpath = opt.options.outpath

    dict_virus, nc_list = virus(sampleid)

    create_umi(sampleid, nc_list)
    
    create_uniq_umi(sampleid, dict_virus)

    count(sampleid, dict_virus)

    cmd = 'cat ' + sampleid + '/single_count/*.tmp.count > ' + sampleid + '/countbyumi.txt'
    runcmd(cmd)
    
    dict_virus_filter = virus_filter(sampleid,outpath)
    
    filter_countbyumi(sampleid, dict_virus_filter, outpath)
    

    
if __name__ == "__main__":
    main()
