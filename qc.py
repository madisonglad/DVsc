import subprocess
import optparse
import os
import pandas as pd

class lf_opt:
    def __init__(self):
        parser = optparse.OptionParser()
        parser.add_option("-s", "--sample", dest="sampleid", help="sample id")
        parser.add_option("-l", "--lengthcut", dest="lengthcut", help="continue length cutoff")
        parser.add_option("-f", "--annotatefile", dest="annotatefile", help="path of annotatefile")

        self.options, self.args = parser.parse_args()
        print(self.options, self.args)

def runcmd(command):
    ret = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=10000)
    if ret.returncode == 0:
        pass
    else:
        print(ret)

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
    
def create_bam(id):
    sorted_bam = id + '/' + id + '.sort.bam'
    dict_ratio = {}
    dict_mappingreads = {}
    with open(id + '/count_ref.txt', 'r') as f:
        for line in f:
            if not line.startswith('NC'):
                continue
            if line.startswith('NC_012920'):
                continue
            virus_name = line.split('\t')[0]
            mapping_reads = line.split('\t')[2]
            if int(mapping_reads) < 5:
                continue
            outfile = id + '/viral_bam/' + virus_name + '.bam'
            
            cmd_viral = 'samtools view -b ' + sorted_bam + ' \'' + virus_name + '\' > ' + outfile
            runcmd(cmd_viral)
            
            tmptxt1 = id + '/viral_bam/' + virus_name +'.tmp.1.txt'
            cmd_sta = 'samtools view ' + outfile + ' | wc -l > ' + tmptxt1
            runcmd(cmd_sta)
            
            with open(tmptxt1) as f:
                for line in f:
                    mapping_num = line.strip()
            cmd_rm1 = 'rm ' + tmptxt1
            runcmd(cmd_rm1)
            
            tmptxt2 = id + '/viral_bam/' + virus_name +'.tmp.2.txt'
            cmd_filter = "samtools view " + outfile + "|awk -F'\t' '{if($5==60) print $0}' | wc -l > " + tmptxt2
            runcmd(cmd_filter)
            with open(tmptxt2) as f:
                for line in f:
                    qualify_num = line.strip()
            cmd_rm2 = 'rm ' + tmptxt2
            runcmd(cmd_rm2)
            #print(virus_name, mapping_num, qualify_num)            

            qualify_ratio = float(qualify_num) / float(mapping_num)
            dict_ratio[virus_name] = qualify_ratio
            
            dict_mappingreads[virus_name] = mapping_reads
            
    return dict_ratio, dict_mappingreads
            
def get_continue_str(orilist):
    length = len(orilist)
    index_count_dict = {}
    for i in range(length):
        key = index_count_dict.keys()
        if key:
            dif = orilist[i] - orilist[i-1]
            if dif == 1:
                index_count_dict[max(key)].append(orilist[i])
            else:
                index_count_dict[i] = [orilist[i]]
        else:
            index_count_dict[i] = [orilist[i]]
            
    maxstr = index_count_dict[max(index_count_dict, key=lambda x: len(index_count_dict[x]))]
    
    return index_count_dict, maxstr

def qc(id, dict_ratio, dict_mappingreads,lengthcut,dict_length):
    outf = open(id + '/filterlist.txt','w')

    for root,dirs,files in os.walk(id + '/viral_bam'):
        for filename in files:
            if filename.endswith('.bam'):
                name = filename.split('.')[0]
                filepath = root + '/' + filename
                covfile = id + '/bed/' + name + '.cov'
                
                cmd = 'samtools depth ' + filepath + ' > ' + covfile

                runcmd(cmd)

                #length = int(name.split('_')[-1])
                length = int(dict_length[name.split('_')[0] + '_' + name.split('_')[1]])

                df = pd.read_csv(covfile, sep='\t', names=['name','pos','depth'])

                cov = df['depth'][ df['depth'] > 0].count() / length

                if length < 15000:
                    cutoff = 200/length
                else:
                    cutoff = 200/15000
                 
                if cov < cutoff:
                    continue
                
                substr = df['pos'][ df['depth'] > 0].tolist()
                
                index, maxstr = get_continue_str(substr)
                
                if len(maxstr) < int(lengthcut):
                    continue                                              
                
                if float(dict_ratio[name]) < 0.8 and int(dict_mappingreads[name]) < 1000:
                    continue
                    
                if float(dict_ratio[name]) < 0.6 and int(dict_mappingreads[name]) < 10000 and int(dict_mappingreads[name]) >= 1000:
                    continue
                
                if float(dict_ratio[name]) < 0.1 and int(dict_mappingreads[name]) >= 10000:
                    continue
                    
                outf.write(name + '\t' + str(dict_mappingreads[name]) + '\t' + str(dict_ratio[name]) + '\t' + str(cov) + '\t' + str(len(maxstr)) + '\n')
                   
def anno(sampleid, dict):
    outf = open(sampleid + 'filter.anno.txt', 'w')
    for root, dirs, files in os.walk(sampleid):
        for file in files:
            if file.startswith('filterlist'):
                filepath = root + '/' + file
                with open(filepath) as f:
                    for line in f:
                        ncid = 'NC_' + line.split('\t')[0].split('_')[1]
                        outf.write(dict[ncid] + '\t' + line)
                                              
def main():
    opt = lf_opt()
    sampleid = opt.options.sampleid
    lengthcut = opt.options.lengthcut
    annotatefile = opt.options.annotatefile
    
    dict_length, dict_name = anno_dict(annotatefile)
    
    print('create bam files start...')
    ratio_dict, mapping_dict = create_bam(sampleid)
    print('create bam files success!')
    
    print('qc start...')
    qc(sampleid, ratio_dict, mapping_dict, lengthcut, dict_length)   
    print('qc success!')
    
    anno(sampleid, dict_name)
    
if __name__ == "__main__":
    main()
