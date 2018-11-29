import sys,os
import re
import argparse

def safe_open(file_name,mode):
    try:
        if not file_name.endswith('.gz'):
            return open(file_name,mode)
        else:
            import gzip
            return gzip.open(file_name,mode)
    except IOError:
        print file_name + ' do not exist!'

def getHet(str):
    if re.search(r'\d',str):
        sample_format = re.search('(\d+)[\/|\|](\d+)',str)
        if sample_format.group(1) == sample_format.group(2):
            return 'Hom'
        else:
            return 'Het'
    else:
        return '.'
                                  
def getRegion(infile,Snum):
    ReLi = []
    ofile = safe_open(infile,'r')
    for line in ofile:
        uline = line.strip().split('\t')
        if 'CHROM' in uline:
            print uline
        #if line.startswith('Priority') or line.startswith('CHROM'):
            chrpos = uline.index('CHROM')
            ppos = uline.index('POS')
            refpos = uline.index('REF')
            altpos = uline.index('ALT')
            formatpos = '.'
            if Snum == '1':
                if 'FORMAT' in uline:
                    formatpos = uline.index('FORMAT')
            continue
        #if type == 'snp':
        else:
            ref = uline[refpos]
            alt = uline[altpos]
            region = str(uline[chrpos])+':'+str(uline[ppos])+'-'+str(uline[ppos])
            sformat = '.'
            if formatpos != '.':
                sformat = getHet(uline[formatpos+1])
            uinfo = '|'.join([region,ref,alt,sformat])
            if (len(ref) != len(alt)) or (not re.search(r'[ATCGatcg]',ref)) or (not re.search(r'[ATCGatcg]',alt)):
                if re.search(r'[ATCGatcg]',ref) or re.search(r'[ATCGatcg]',alt):  ##both '.'
                    mtype = 'indel'
                else:
                    mtype = 'snp'
            else:
                mtype = 'snp'
            if uinfo not in ReLi:
                ReLi.append(uinfo+'|'+mtype)
    return ReLi

def GetMpileup(ulist,infile,Sam,Snum):
    cmdLi = []
    #count = 0
    for i in ulist:
        info = i.strip().split('|')
        regionInfo = info[0]
        refInfo = info[1]
        altInfo = info[2]
        genotypeInfo = info[3]
        mtype = info[-1]
        Mcmd = ''
        ##4,5,6 for the first samples
        ##add type recognize and mpileup for ~5bp region, do not out put ref,alt and information can't get
        if mtype == 'indel':
            thisSite = regionInfo.split(':')[1].split('-')[0]
            if int(thisSite)-5 > 0:
                tmpregions = regionInfo.split(':')[0]+':'+str(int(thisSite)-5)+'-'+str(int(thisSite)-1)
            else:
                tmpregions = regionInfo.split(':')[0]+':'+'0-'+str(int(thisSite)-1)
            tmpregione =  regionInfo.split(':')[0]+':'+str(int(thisSite)+1)+'-'+str(int(thisSite)+5)
            Mcmd += 'samtools mpileup -r %s  -q 1 -C 50 -t DP,DV -m 2 -F 0.002 -f %s %s |  awk -F "\\t" \'{print "%s\\t"$1"\\t"$2"\\t%s\\t%s\\t%s\\t"$4"\\t"$5"\\t"$6' %(tmpregions, reffasta, Sam,'-5 tmp indel','/','/','/')
            if Snum > 1:
                for i in range(2,Snum+1):
                    Mcmd += '"\\t"$%s' %str(6+(i-1)*3-2)
                    Mcmd += '"\\t"$%s' %str(6+(i-1)*3-1)
                    Mcmd += '"\\t"$%s' %str(6+(i-1)*3)
            Mcmd += '}\' >> %s\n' %infile
        Mcmd += 'samtools mpileup -r %s  -q 1 -C 50 -t DP,DV -m 2 -F 0.002 -f %s %s |  awk -F "\\t" \'{print "%s\\t"$1"\\t"$2"\\t%s\\t%s\\t%s\\t"$4"\\t"$5"\\t"$6' %(regionInfo, reffasta, Sam,mtype,refInfo,altInfo,genotypeInfo)
        if Snum > 1:
            for i in range(2,Snum+1):
                Mcmd += '"\\t"$%s' %str(6+(i-1)*3-2)
                Mcmd += '"\\t"$%s' %str(6+(i-1)*3-1)
                Mcmd += '"\\t"$%s' %str(6+(i-1)*3)
        Mcmd += '}\' >> %s' %infile
        if mtype == 'indel':
            Mcmd += '\nsamtools mpileup -r %s  -q 1 -C 50 -t DP,DV -m 2 -F 0.002 -f %s %s |  awk -F "\\t" \'{print "%s\\t"$1"\\t"$2"\\t%s\\t%s\\t%s\\t"$4"\\t"$5"\\t"$6' %(tmpregione, reffasta, Sam,'+5 tmp indel','/','/','/')
            if Snum > 1:
                for i in range(2,Snum+1):
                    Mcmd += '"\\t"$%s' %str(6+(i-1)*3-2)
                    Mcmd += '"\\t"$%s' %str(6+(i-1)*3-1)
                    Mcmd += '"\\t"$%s' %str(6+(i-1)*3)
            Mcmd += '}\' >> %s' %infile
        cmdLi.append(Mcmd)
        #cmdLi.append('samtools mpileup -r %s  -q 1 -C 50 -t DP,DV -m 2 -F 0.002 -f /PUBLIC/database/HUMAN/genome/Human/human_g1k_v37_decoy.fasta  %s |  awk -F "\\t" \'{print "%s\\t"$1"\\t"$2"\\t"$3"\\t%s\\t%s\\t"$4"\\t"$5"\\t"$6"\\t"$7"\\t"$8"\\t"$9"\\t"$10"\\t"$11"\\t"$12}\' >> %s' %(regionInfo,Sam,type,altInfo,genotypeInfo,file))
    return cmdLi
    
from argparse import RawTextHelpFormatter
parser = argparse.ArgumentParser(description = 'Get mpileup format result.')
#parser.add_argument('-P', '--Propositus', metavar = 'String', help="The Propositus's sample name.", required=True)
#parser.add_argument('-s', '--samples',metavar = 'String', help="The samples name.", required=True)
parser.add_argument('-B', '--Bams', metavar = 'String/File', help="All samples bam file need to generate mpileup file.Comma separated or per line in a file. eg: test1.bam,test2.bam,test3.bam", required=True)
#parser.add_argument(
parser.add_argument('-input', '--input', metavar = 'File', help="Annotated file.", required=True)
parser.add_argument('-r', '--reffasta', metavar = 'File', help="The reference fasta file.", required=True)
parser.add_argument('-output','--outfile' ,  metavar = 'File', help="Out put result of mpileup", required=False)
parser.add_argument('-O', '--outshell', metavar = 'File', help="Output run shell.", required=False)
argv = vars(parser.parse_args())

inputf = argv['input'].strip()
outf = argv['outfile'].strip()
reffasta = argv['reffasta']

##bams
sam = argv['Bams']
if os.path.isfile(sam) and (not sam.endswith('.bam')):
    tempSam = []
    for line in safe_open(sam,'r'):
        if line.strip() not in tempSam:
            tempSam.append(line.strip())
    Sam = ' '.join(tempSam)
else:
    tempSam = []
    for i in sam.strip().split(','):
        tempSam.append(i.strip())
    print tempSam
    Sam = ' '.join(tempSam)

Snum = len(Sam.strip().split(' '))

print Snum,' samples: ', Sam
##output shell
if argv['outshell']:
    output = os.path.abspath(argv['outshell'])
else:
    output = 'runMpileup.sh'

outhead = 'MutationType\tChr\tPosition\tRef\tAlt\tGenotype'
for samp in Sam.strip().split(' '):
    sname = os.path.basename(samp).split('.')[0]
    outhead += '\t%s_ReadsNumber\t%s_Base\t%s_BaseQuality' %(sname,sname,sname)
with safe_open(outf,'w') as file:
    file.write(outhead+'\n')

OutLi = []
outputshell = open(output,'w')
region = getRegion(inputf,Snum)
OutLi = GetMpileup(region,outf,Sam,Snum)
outputshell.write(' && \\\n'.join(OutLi)+'\n')
outputshell.close()
os.system('sh %s' %output)
