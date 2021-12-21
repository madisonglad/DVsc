import subprocess
import optparse
import os
import pandas as pd
from collections import defaultdict


class lf_opt:
    def __init__(self):
        parser = optparse.OptionParser()
        parser.add_option("-s", "--sample", dest="sampleid", help="sample id")

        self.options, self.args = parser.parse_args()
        print(self.options, self.args)

def runcmd(command):
    ret = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=10000)
    if ret.returncode == 0:
        pass
    else:
        print(ret)

def filter_list(id):
    with open(id + 'filter.anno.txt') as f:
        dict_virus = defaultdict(list)
        nc_list = []
        for line in f:
            virus_name = line.split('\t')[0].replace(' ', '_').replace('/', '_').replace('(', '_').replace(')', '_')
            virus_nc = line.split('\t')[1]
            dict_virus[virus_name].append(virus_nc)
            nc_list.append(virus_nc)
    return dict_virus, nc_list


def viral_count(id, nc_list, dict_virus):
    path = id + '/viral_bam'
    for root, dirs, files in os.walk(path):
        for file in files:
            if file.split('.')[0] in nc_list:
                cmd = "samtools view " + id + '/viral_bam/' + file + " | awk '{print $1}' > " + id + '/tmp/' + file.split('.')[0] + '.ncid'
                runcmd(cmd)

    dict = {}
    for virus in dict_virus:
        virus = virus.replace(' ', '_').replace('/', '_').replace('(', '_').replace(')', '_')
        cmd = 'mkdir -p ' + id + '/tmp2/' + virus
        runcmd(cmd)

        for nc in dict_virus[virus]:
            cmd = 'cp ' + id + '/tmp/' + nc + '.ncid ' + id + '/tmp2/' + virus + '/'
            runcmd(cmd)
        cmd = 'cat ' + id + '/tmp2/' + virus +'/*.ncid > ' + id + '/tmp2/' + virus + '.merge.ncid'
        runcmd(cmd)
        cmd = 'sort ' + id + '/tmp2/'  + virus + '.merge.ncid | uniq > ' + id + '/single_count/' + virus + '.uniq.ncid'
        runcmd(cmd)
        cmd = 'wc -l ' + id + '/single_count/' + virus + '.uniq.ncid > ' + id + '/single_count/' + virus + '.uniq.ncid.count'
        runcmd(cmd)
        
        with open(id + '/single_count/' + virus + '.uniq.ncid.count') as f:
            for line in f:
                count = line.split()[0]
                dict[virus] = count
                
        cmd = 'rm ' + id + '/single_count/' + virus + '.uniq.ncid.count' 
        runcmd(cmd)                

    outf = open(id + '/virus_count.txt', 'w')
    for virus in dict_virus:
        outf.write(id + '\t' + virus + '\t' + dict[virus] + '\n')



def main():
    opt = lf_opt()
    sampleid = opt.options.sampleid

    dict_virus, nc_list = filter_list(sampleid)

    viral_count(sampleid, nc_list, dict_virus)
    

    
if __name__ == "__main__":
    main()
