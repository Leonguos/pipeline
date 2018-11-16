#!/usr/bin/python
# vim: set fileencoding=utf-8 :

#modify: 20170606 for some description
import sys
import os
import re
import codecs
import argparse
from argparse import RawTextHelpFormatter


BASE_DIR = os.path.dirname(os.path.abspath(__file__))

def safe_open(file,mode):
    file = os.path.abspath(file)
    if not os.path.exists(file):
        exit('%s is not exists.' %file)
    elif file.endswith('.gz'):
        import gzip
        return gzip.open(file,mode)
    else:
        return open(file,mode)

def getpos(name,title,name2=''):
    if name in title:
        pos = title.index(name)
    elif name2 != '' and name2 in title:
        pos = title.index(name2)
    else:
        if name2 != '':
           exit('Warning: %s and %s not in title.' %(name,name2))
        else:
           exit('Warning: %s not in title.' %name)
    return pos

def getpos2(name,title,name2=''):
    ntitle = [i.lower() for i in title]
    if name.lower() in ntitle:
        pos = ntitle.index(name.lower())
    elif name2 != '' and name2.lower() in ntitle:
        pos = ntitle.index(name2.lower())
    else:
        if name2 != '':
           exit('Warning: %s and %s not in title.' %(name,name2))
        else:
           exit('Warning: %s not in title.' %name)
    return pos

def geth(sample,formatinfo,sampleformat,pos,mpileupDic):
    if not re.match(r'\d+', sampleformat) or re.match(r'0/0',sampleformat):
        if mpileupDic == {}:
            return ['']
        else:
            if sample in mpileupDic.keys():
                if mpileupDic[sample].get(pos):
                    if mpileupDic[sample][pos][1] > 0:
                        return ['',mpileupDic[sample][pos][0],mpileupDic[sample][pos][1]]
                    else:
                        return ['/',mpileupDic[sample][pos][0],mpileupDic[sample][pos][1]]
                else:
                    return ['/','/','/']
        ##根据mpileup文件中的点判断是没覆盖还是没突变，没覆盖返回'/'，没突变返回''和AF,深度
    else:
        dpdic = dict(zip(formatinfo.split(':'), sampleformat.split(':')))
        if ('DV' in dpdic) and re.match(r'\d+',dpdic['DV']) and re.match(r'\d+',dpdic['DP']) and int(dpdic['DP']) > 0:
            dvinfo = float(dpdic['DV'])
            dpinfo = int(dpdic['DP'])
            AF = dvinfo/dpinfo
        else:
            if ('AD' in dpdic) and  re.match(r'\d+',dpdic['AD']) and re.match(r'\d+',dpdic['DP']) and int(dpdic['DP']) > 0:
                Adp = 0
                for adp in range(1,len(dpdic['AD'].split(','))):
                    if dpdic['AD'].split(',')[adp]!=".":
                        Adp += float(dpdic['AD'].split(',')[adp])
            else:
                exit('Wrong FORMAT infomation %s' %dpdic)
            dpinfo = int(dpdic['DP'])
            dvinfo = float(Adp) 
            AF = dvinfo/dpinfo
        genotype = dpdic['GT']
        sample_format = re.search('(\d+)[\/|\|](\d+)', genotype)
        if sample_format.group(1) == sample_format.group(2) and sample_format.group(1) != '0':
            return ['hom', AF, dpinfo]
        elif sample_format.group(1) == sample_format.group(2) and sample_format.group(1) == '0':
            return ['', 0, dpinfo]
        else:
            if AF > 0.8:
                return ['hom', AF, dpinfo]
            else:
                return ['het', AF, dpinfo]

##def get gene list based on provided
def getgli(genef):
    ##geneli is provided gene list.And convert all of provided gene to upper case.
    geneli = []
    ugene = []
    for line in safe_open(genef,'r'):
        if line.startswith('#'):continue
        genetmp = line.strip().split(',')
        for ge in genetmp:
            if ge not in geneli:
                geneli.append(ge.upper())  ##convert to upper case 
                
    ##genedb is a file contain standard gene name and id and aliases.            
    genedb = safe_open(os.path.join(BASE_DIR, 'GeneDic.xls'), 'r')
    ##geneali is the dic of {same standard name: ID0\tstandard1 | ID1\tstandard1,aliases1,aliases2 | ID2\tstandard2,aliases1,aliases2} .And the gene of genedb is upper case.
    global geneali
    geneali = {}
    for line in genedb:
        if line.startswith('#'):continue
        line = line.strip().split(':')
        geneali[line[0].strip()] = line[1].strip()
     
    ##tfugl: tmp finally used gene list , get finally use gene list from provided gene list and standard-aliases geneali. 
    tfugl = []
    for pgene in geneli:                
        if pgene in geneali.keys():
            tfugl.append(geneali[pgene])
        else:
            tmp = []
            for k in geneali.keys():
                ##geneali[k].split(',') are only aliases
                if pgene in geneali[k].split(','):
                    if geneali[k] not in tmp:
                        tmp.append(geneali[k])
            if len(tmp) > 1:
                print '%s is aliases more than one gene. %s' %(pgene,' | '.join(tmp))
            if len(tmp) == 1:
                print '%s is not standard but aliases of: %s' %(pgene,' | '.join(tmp))
            if len(tmp) == 0:
                print '%s is not standard gene name nor aliases.' %pgene
                ugene.append(pgene)
                tmp.append('Unknown\t'+pgene)
            tfugl.append(' | '.join(tmp))
    
    ##sfugl for standard finally used gene list.
    ##afugl for aliases finally used gene list.
    ##we default consider the gene list provided is standard gene name.
    sfugl = []
    afugl = []
    for g in tfugl:
        tmp = g.strip().split(' | ')
        for n in range(0,len(tmp)):
            ttmp = tmp[n].strip().split('\t')[1]
            gene = ttmp.strip().split(',')
            if len(gene) == 1:
                sfugl.append(gene[0])
            else:
                sfugl.append(gene[0])
                afugl.extend(gene[1:])
    sfugl = list(set(sfugl))
    afugl = list(set(afugl))
    ugene = list(set(ugene))
    return geneli,sfugl,afugl,ugene
    
##def annotated information
##samin p.name.xian;p.name.relation;n.name.relation
def getinfo(suffix,fid,normalLi,patientLi,mpileupDic,checkConfidence,other = []):
    otherm = other
    input = safe_open(suffix,'r')
    allinfo = []
    Nnum = len(normalLi)
    Npat = len(patientLi)
    AllLi = []
    AllLi.extend(normalLi)
    AllLi.extend(patientLi)
    for line in input:
        cflag = 0
        mflag = 0
        wflag = 0
        #if line.startswith('CHROM') or line.startswith('Priority'):
        if line.startswith('GeneName') or line.startswith('GeneCandidate'):
            SamIndex = []
            #NormIndex = []
            #PaIndex = []
            head = line.strip().split('\t')
            title = line.strip().split('\t')
            gpos = getpos('GeneName',head)
            pripos = getpos('Priority',head)#Priority
            cpos = getpos('CHROM',head)
            ppos = getpos('POS',head)
            rpos = getpos('REF',head)
            apos = getpos('ALT',head)
            fpos = getpos('FORMAT',head)
            for sam in AllLi:
                SamIndex.append((getpos(sam,head),sam))
            #title[1:1] = ['Confidence\tCheckLog']
            title[pripos:pripos] = ['Confidence\tCheckLog']
            allinfo.append(title)
            continue
        else:
            SamInfo = []
            uline = line.strip().split('\t')
            refinfo = uline[rpos]
            altinfo = uline[apos]
            if len(refinfo) != len(altinfo):
                mtype = 'indel'
            else:
                mtype = 'snp'
            posinfo = '_'.join([uline[cpos],uline[ppos],uline[rpos],uline[apos]])
            chrpos = 'chr'+uline[cpos]+':'+uline[ppos]
            nucpos = uline[rpos]+'>'+uline[apos]
            chet = ''
            ##得到正常人和患者的突变信息，AF,DP信息，未覆盖的显示为'/'
            for dex in SamIndex:
                gtinfo = uline[dex[0]]
                gttype = geth(head[dex[0]],uline[fpos],gtinfo,posinfo,mpileupDic)  #根据mpileup文件判断纯杂合还是未覆盖, [type,AF,depth]
                gttype.append(gtinfo)
                gttype.append(dex[1])
                SamInfo.append(gttype)  ## [type,AF,depth,sample_format,sample_format_index]
            log = []
            if mpileupDic != {}:
                if checkConfidence == 'Y':
                    for samtype in SamInfo: ##samtype [type,AF,depth,sam]
                        if not re.match(r'\d+',samtype[3].split(':')[0]) or re.match(r'0/0',samtype[3].split(':')[0]): 
                            #样本原始没有检出有突变,判断是不是可信的没有检出突变，如果有遗传模式，根据遗传模式判断是不是可信。
                            #if float(samtype[1]) <0.25 and float(samtype[2]) >= udepth:  ##当前位点没有突变，继续判断上下游有没有突变。
                            #if posinfo == '22_20779975_GC_G':
                            #    print 'test', samtype[-1],mpileupDic[samtype[-1]][posinfo]
                            if (samtype[1] != '/' and float(samtype[1]) < 0.25) or samtype[1] == '/' :
                                if mtype == 'indel':
                                    tmpflag = 0
                                    for num in range(int(uline[ppos])-bp,int(uline[ppos])):
                                        npos = '_'.join([uline[cpos],str(num),'/','/'])
                                        #if posinfo == '22_20779975_GC_G':
                                        #    print 'test2', samtype[-1],mpileupDic[samtype[-1]][npos]
                                        if mpileupDic[samtype[-1]][npos][2] == 'indel' and mpileupDic[samtype[-1]][npos][0] != '/' \
                                            and float(mpileupDic[samtype[-1]][npos][0]) >= af \
                                            and int(mpileupDic[samtype[-1]][npos][1]) >= udepth:  
                                            ##没有检出indel突变，如果上下游检出可信的突变，则认为这个indel可能不可信
                                            if argv['hmode']:
                                                if hmode == 'R':
                                                    if (samtype[-1] in normalLi) and float(mpileupDic[samtype[-1]][npos][0]) < 0.75:
                                                        tmpflag += 0
                                                        mflag  += 1
                                                        log.append('Notice: %s of %s may be a het mutation (nearby indel %s)' \
                                                        %(posinfo,samtype[-1],mpileupDic[samtype[-1]][npos]))
                                                    elif (samtype[-1] in patientLi) and float(mpileupDic[samtype[-1]][npos][0]) >= 0.75:
                                                        tmpflag += 0
                                                        mflag += 1
                                                        log.append('Notice: %s of %s may be a hom mutation (nearby indel %s)' \
                                                        %(posinfo,samtype[-1],mpileupDic[samtype[-1]][npos]))
                                                    else:
                                                        tmpflag += 1
                                                        wflag += 1
                                                        log.append('Notice: %s of %s may not be a confidence mutation type for recessive mode (nearby indel %s)' \
                                                        %(posinfo,samtype[-1],mpileupDic[samtype[-1]][npos]))
                                                else:##显性遗传模式
                                                    if samtype[-1] in normalLi:
                                                        tmpflag += 1
                                                        wflag += 1
                                                        ##没有检出indel突变，如果上下游检出可信的突变，则认为这个indel可能不可信
                                                        log.append('%s of %s may not be a confidence ref (nearby indel　%s)' \
                                                        %(posinfo,samtype[-1],mpileupDic[samtype[-1]][npos])) 
                                                    elif (samtype[-1] in patientLi) and float(mpileupDic[samtype[-1]][npos][0]) < 0.75:
                                                        tmpflag += 0
                                                        mflag += 1
                                                        log.append('Notice: %s of %s may be a het mutation (nearby indel %s)' \
                                                        %(posinfo,samtype[-1],mpileupDic[samtype[-1]][npos]))
                                                    else:
                                                        tmpflag += 1
                                                        wflag += 1
                                                        log.append('Notice: %s of %s may not be a confidence mutation type for dominant mode (nearby indel %s)' \
                                                        %(posinfo,samtype[-1]),mpileupDic[samtype[-1]][npos])
                                            else:
                                                tmpflag += 1
                                                wflag += 1
                                    for num in range(int(uline[ppos])+1, int(uline[ppos])+bp):
                                        npos = '_'.join([uline[cpos],str(num),'/','/'])
                                        if mpileupDic[samtype[-1]][npos][2] == 'indel' and mpileupDic[samtype[-1]][npos][0] != '/' \
                                            and float(mpileupDic[samtype[-1]][npos][0]) >= af and mpileupDic[samtype[-1]][npos][1] >= udepth:
                                            if argv['hmode']:
                                                if hmode == 'R':
                                                    if (samtype[-1] in normalLi) and float(mpileupDic[samtype[-1]][npos][0]) < 0.75:
                                                        tmpflag += 0
                                                        mflag += 1
                                                        log.append('Notice: %s of %s may be a het mutation (nearby indel %s)' \
                                                        %(posinfo,samtype[-1],mpileupDic[samtype[-1]][npos]))
                                                    elif (samtype[-1] in patientLi) and float(mpileupDic[samtype[-1]][npos][0]) >= 0.75:
                                                        tmpflag += 0
                                                        mflag += 1
                                                        log.append('Notice: %s of %s may be a hom mutation (nearby indel %s)' \
                                                        %(posinfo,samtype[-1],mpileupDic[samtype[-1]][npos]))
                                                    else:
                                                        tmpflag += 1
                                                        wflag += 1
                                                        log.append('Notice: %s of %s may not be a confidence mutation type for recessive mode (nearby indel %s)' \
                                                        %(posinfo,samtype[-1],mpileupDic[samtype[-1]][npos]))
                                                else:##显性遗传模式
                                                    if samtype[-1] in normalLi:
                                                        tmpflag += 1
                                                        wflag += 1
                                                        log.append('%s of %s may not be a confidence ref (nearby indel　%s)' \
                                                         %(posinfo,samtype[-1],mpileupDic[samtype[-1]][npos]))
                                                    elif samtype[-1] in patientLi and float(mpileupDic[samtype[-1]][npos][0]) < 0.75:
                                                        tmpflag += 0
                                                        mflag += 1
                                                        log.append('Notice: %s of %s may be a het mutation (nearby indel %s)' \
                                                        %(posinfo,samtype[-1],mpileupDic[samtype[-1]][npos]))
                                                    else:
                                                        tmpflag += 1
                                                        wflag += 1
                                                        log.append('Notice: %s of %s may not be a confidence mutation type for dominant mode (nearby indel %s)' \
                                                        %(posinfo,samtype[-1],mpileupDic[samtype[-1]][npos]))
                                            else:
                                                tmpflag += 1
                                                wflag += 2
                                    if tmpflag == 0:
                                        cflag += 0  ##mpileupDic[pat][altinfo] = [mpileupaf,mpileupdp,mtype]
                                    else:
                                        cflag += 1
                                        if not argv['hmode']:
                                            log.append('%s of %s may not be a confidence ref %s (nearby indel)' \
                                            %(posinfo,samtype[-1],mpileupDic[samtype[-1]][posinfo]))
                                            print posinfo, 'in',  samtype[-1], 'may not be a confidence ref (nearby indel)'
                            #elif samtype[2] < udepth:
                            #    cflag += 1
                            #    wflag += 1
                            #    log.append('%s of %s may not be a confidence ref %s' %(posinfo, samtype[-1],mpileupDic[samtype[-1]][posinfo]))                                
                            else: ##虽然变异检测没有，但当前位点其实有突变,或者小于udepth的有突变
                                if argv['hmode']:
                                    if hmode == 'R' and mpileupDic[samtype[-1]][posinfo][0] != '/':
                                        if (samtype[-1] in normalLi) and  mpileupDic[samtype[-1]][posinfo][0] != '/' and \
                                            float(mpileupDic[samtype[-1]][posinfo][0]) >= af and \
                                            float(mpileupDic[samtype[-1]][posinfo][0]) < 0.75 \
                                            and int(mpileupDic[samtype[-1]][posinfo][1]) >= udepth:
                                            cflag += 0
                                            mflag += 1
                                            log.append('Notice: %s of %s may be a het mutation %s' \
                                            %(posinfo,samtype[-1],mpileupDic[samtype[-1]][posinfo]))
                                        elif (samtype[-1] in patientLi) and float(mpileupDic[samtype[-1]][posinfo][0]) >= 0.75 and \
                                            int(mpileupDic[samtype[-1]][posinfo][1]) >= udepth:
                                            cflag += 0
                                            mflag += 1
                                            log.append('Notice: %s of %s may be a hom mutation %s' \
                                            %(posinfo,samtype[-1],mpileupDic[samtype[-1]][posinfo]))
                                        else:
                                            cflag += 1
                                            wflag += 1
                                            log.append('%s of %s may not be a confidence mutation type for recessive mode %s' \
                                            %(posinfo, samtype[-1],mpileupDic[samtype[-1]][posinfo]))
                                            print posinfo, 'in',  samtype[-1], 'may not be a confidence ref',mpileupDic[samtype[-1]][posinfo]
                                    elif mpileupDic[samtype[-1]][posinfo][0] != '/':
                                        if samtype[-1] in normalLi:
                                            cflag += 1
                                            wflag += 1
                                            log.append('%s of %s may not be a confidence ref %s' \
                                            %(posinfo,samtype[-1]),mpileupDic[samtype[-1]][posinfo])
                                        elif (samtype[-1] in patientLi) and float(mpileupDic[samtype[-1]][posinfo][0]) >= af and \
                                            float(mpileupDic[samtype[-1]][posinfo][0]) < 0.75 and \
                                            int(mpileupDic[samtype[-1]][posinfo][1]) >= udepth:
                                            cflag += 0
                                            mflag += 1
                                            log.append('Notice: %s of %s may be a het mutation %s' \
                                            %(posinfo,samtype[-1],mpileupDic[samtype[-1]][posinfo]))
                                        else:
                                            cflag += 1
                                            wflag += 1
                                            log.append('%s of %s may not be a confidence mutation type for dominant mode %s' \
                                            %(posinfo, samtype[-1],mpileupDic[samtype[-1]][posinfo]))
                                            print posinfo, 'in',  samtype[-1], 'may not be a confidence mutation type for dominant mode ',mpileupDic[samtype[-1]][posinfo]
                                else:
                                    cflag += 1
                                    wflag += 1
                                    log.append('%s of %s may not be a confidence ref %s (allele frequency)' \
                                    %(posinfo, samtype[-1],mpileupDic[samtype[-1]][posinfo]))
                                    print posinfo, 'in',  samtype[-1], 'may not be a confidence ref',mpileupDic[samtype[-1]][posinfo]
                        else:##检出突变，看是不是可信突变
                            #if posinfo == '6_1612017_A_ACGG':
                            #    print 'test3', samtype[-1], mpileupDic[samtype[-1]][posinfo]
                            ##检出突变，并且mpileup有才继续判断，否则认为pileup的时候没有，是以为reads质量问题
                            if mpileupDic[samtype[-1]].get(posinfo):
                                if mpileupDic[samtype[-1]][posinfo][0] != '/' and float(mpileupDic[samtype[-1]][posinfo][0]) >= af \
                                    and int(mpileupDic[samtype[-1]][posinfo][1]) >= udepth:  ##检出可信突变
                                    if argv['hmode']:
                                        if (samtype[0] == 'hom' and float(mpileupDic[samtype[-1]][posinfo][0]) >= 0.75) or \
                                            (samtype[0] == 'het' and float(mpileupDic[samtype[-1]][posinfo][0]) >= af and \
                                            float(mpileupDic[samtype[-1]][posinfo][0]) < 0.75): ##有突变且符合hom & het
                                            cflag += 0
                                        else:
                                            cflag += 1
                                            wflag += 1
                                            log.append('%s of %s may not be a confidence mutation type corresponding to the het/hom genotype %s' \
                                            %(posinfo, samtype[-1], mpileupDic[samtype[-1]][posinfo]))
                                            print posinfo, 'in',  samtype[-1], 'may not be a confidence mutation', mpileupDic[samtype[-1]][posinfo]
                                    else:
                                        cflag += 0
                                        if (samtype[0] == 'hom' and float(mpileupDic[samtype[-1]][posinfo][0]) >= 0.75) or \
                                            (samtype[0] == 'het' and float(mpileupDic[samtype[-1]][posinfo][0]) >= af and \
                                            float(mpileupDic[samtype[-1]][posinfo][0]) < 0.75): ##有突变且符合hom & het
                                            mflag += 0
                                        else:
                                            mflag += 1
                                            log.append('Notice: %s of %s may not be a confidence mutation type corresponding to the het/hom genotype %s' \
                                            %(posinfo, samtype[-1], mpileupDic[samtype[-1]][posinfo]))
                                            print posinfo, 'in',  samtype[-1], 'may not be a confidence mutation', mpileupDic[samtype[-1]][posinfo]
                                else:
                                    if mtype == 'indel':
                                        tmpflag2 = 0
                                        tmpflag3 = 0
                                        for num in range(int(uline[ppos])-bp,int(uline[ppos])):
                                            npos = '_'.join([uline[cpos],str(num),'/','/'])
                                            if mpileupDic[samtype[-1]][npos][2] == 'indel' and mpileupDic[samtype[-1]][npos][0] != '/' and \
                                                float(mpileupDic[samtype[-1]][npos][0]) >= af and int(mpileupDic[samtype[-1]][npos][1]) >= udepth:##上下游有可信的indel突变 
                                                if argv['hmode']:
                                                    if hmode == 'R' and (samtype[-1] in patientLi) and float(mpileupDic[samtype[-1]][npos][0]) >= 0.75:
                                                        tmpflag2 += 0
                                                        tmpflag3 += 1
                                                        print 'Nearby confidence hom mutation at %s of %s: %s' \
                                                        %(npos, samtype[-1],mpileupDic[samtype[-1]][npos])
                                                        log.append('Nearby confidence hom mutation at %s of %s: %s' \
                                                        %(npos, samtype[-1],mpileupDic[samtype[-1]][npos]))
                                                    elif hmode == 'R' and (samtype[-1] in normalLi) and float(mpileupDic[samtype[-1]][npos][0]) < 0.75:
                                                        tmpflag2 += 0
                                                        tmpflag3 += 1
                                                        print 'Nearby confidence het mutation at %s of %s: %s' \
                                                        %(npos, samtype[-1],mpileupDic[samtype[-1]][npos])
                                                        log.append('Nearby confidence het mutation at %s of %s: %s' \
                                                        %(npos, samtype[-1],mpileupDic[samtype[-1]][npos]))
                                                    elif hmode == 'D' and (samtype[-1] in patientLi) and float(mpileupDic[samtype[-1]][npos][0]) < 0.75:
                                                        tmpflag2 += 0
                                                        tmpflag3 += 1
                                                        print 'Nearby confidence het mutation at %s of %s: %s' \
                                                        %(npos, samtype[-1],mpileupDic[samtype[-1]][npos])
                                                        log.append('Nearby confidence het mutation at %s of %s: %s' \
                                                        %(npos, samtype[-1],mpileupDic[samtype[-1]][npos]))
                                                    else:
                                                        tmpflag2 += 1
                                                else:
                                                    tmpflag2 += 0
                                                    tmpflag3 += 1
                                                    log.append('Nearby confidence mutation at %s of %s: %s' \
                                                    %(npos, samtype[-1],mpileupDic[samtype[-1]][npos])) 
                                            else:
                                                tmpflag2 += 1
                                        for num in range(int(uline[ppos])+1,int(uline[ppos])+3):
                                            npos = '_'.join([uline[cpos],str(num),'/','/'])
                                            if mpileupDic[samtype[-1]][npos][2] == 'indel' and mpileupDic[samtype[-1]][npos][0] != '/'  \
                                                and float(mpileupDic[samtype[-1]][npos][0]) >= af and \
                                                int(mpileupDic[samtype[-1]][npos][1]) >= udepth:
                                                if argv['hmode']:
                                                    if hmode == 'R' and (samtype[-1] in patientLi) and float(mpileupDic[samtype[-1]][npos][0]) >= 0.75: 
                                                        tmpflag2 += 0
                                                        tmpflag3 += 1 
                                                        print 'Nearby confidence hom mutation at %s of %s: %s' \
                                                        %(npos, samtype[-1],mpileupDic[samtype[-1]][npos])
                                                        log.append('Nearby confidence hom mutation at %s of %s: %s' \
                                                        %(npos, samtype[-1],mpileupDic[samtype[-1]][npos]))
                                                    elif hmode == 'R' and (samtype[-1] in normalLi) and float(mpileupDic[samtype[-1]][npos][0]) < 0.75:
                                                        tmpflag2 += 0
                                                        tmpflag3 += 1
                                                        print 'Nearby confidence het mutation at %s of %s: %s' \
                                                        %(npos, samtype[-1],mpileupDic[samtype[-1]][npos])
                                                        log.append('Nearby confidence het mutation at %s of %s: %s' \
                                                        %(npos, samtype[-1],mpileupDic[samtype[-1]][npos]))
                                                    elif hmode == 'D' and (samtype[-1] in patientLi) and float(mpileupDic[samtype[-1]][npos][0]) < 0.75:
                                                        tmpflag2 += 0
                                                        tmpflag3 += 1
                                                        print 'Nearby confidence het mutation at %s of %s: %s' \
                                                        %(npos, samtype[-1],mpileupDic[samtype[-1]][npos])
                                                        log.append('Nearby confidence het mutation at %s of %s: %s' \
                                                        %(npos, samtype[-1],mpileupDic[samtype[-1]][npos]))
                                                    else:
                                                        tmpflag2 += 1
                                                else:
                                                    tmpflag2 += 0
                                                    tmpflag3 += 1
                                                    log.append('Nearby confidence mutation at %s of %s: %s' \
                                                    %(npos, samtype[-1],mpileupDic[samtype[-1]][npos]))
                                            else:
                                                tmpflag2 += 1
                                        if tmpflag2 == 0:  ##虽然这个位置mpileup的结果显示没突变，但是上下游有可信突变，也认为是可信突变。
                                            cflag += 0
                                            mflag += 1
                                            print 'Nearby confidence mutation:', 'at ',posinfo, 'of ',  samtype[-1]
                                        else:
                                            if tmpflag3 == 0:  ##一个样本，虽然本身突变不可靠，但是indel上下游有可靠的时候，不再输出自己不可靠
                                                cflag += 1
                                                wflag += 1
                                                log.append('%s of %s may not be a confidence mutation %s' \
                                                %(posinfo, samtype[-1], mpileupDic[samtype[-1]][posinfo]))
                                                print posinfo, 'in',  samtype[-1], 'may not be a confidence mutation', mpileupDic[samtype[-1]][posinfo]
                                    else:
                                        cflag += 1
                                        wflag += 1
                                        log.append('%s of %s may not be a confidence mutation %s' %(posinfo, samtype[-1], mpileupDic[samtype[-1]][posinfo]))
                                        print posinfo, 'in',  samtype[-1], 'may not be a confidence mutation', mpileupDic[samtype[-1]][posinfo]
                            else:##mpileup文件没有这个点，可能是由于-C50 参数被限制掉了。
                                cflag += 1
                                wflag += 1
                                log.append('%s of %s may not be a confidence mutation because of mapping quality.' %(posinfo, samtype[-1]))
                    ##cflag == 0的时候是最理想的时候
                    if cflag == 0 and mflag == 0 and wflag == 0:
                        confidenceLevel = "H"
                    ##nearby indel或者遗传模式不符合的时候 
                    elif cflag == 0 and mflag > 0 and wflag == 0:
                        confidenceLevel = "M"
                    else:
                        confidenceLevel = "L"
                else:
                    confidenceLevel = ''
            else:
                confidenceLevel = ''
        #uline[1:1] = [' | '.join(log)]
        #uline[1:1] = [confidenceLevel]
        uline[pripos:pripos] = [' | '.join(log)]
        uline[pripos:pripos] = [confidenceLevel]
        allinfo.append(uline)
    return allinfo

parser = argparse.ArgumentParser(description = ' .\nContact: yuhuan@novogene.com', formatter_class = RawTextHelpFormatter)
parser.add_argument('-R','--relation', metavar = 'File', required = True, help = 'The sample info file include #FamilyID       SampleID        SEX     Normal/Patient')
parser.add_argument('-F','--familyID', metavar = 'String', required = True, help = 'The familyID want to do')
parser.add_argument('-I','--input', metavar = 'String',required = True,help = 'Input file.And the samples in the relation file must in the first line of the input file.')
parser.add_argument('-O', '--out',  metavar = 'String',  required = False, help = 'Outputfile prefix.' )
parser.add_argument('-Pileup','--mpileup',metavar = 'File',required = False, help = 'The mpileup file help to show confidence.Alternative with -mvcf')
parser.add_argument('-mvcf','--mvcf',metavar = 'File',required = False, help = 'The mpileup vcf file help to show confidence.Alternative with -Pileup')
parser.add_argument('-Check','--CheckConfidende',metavar = 'String',required = False, help = 'Show the confidence information or not when -Pileup used.',choices = ['Y','N'],default = 'N')
#parser.add_argument('-g','--genelist',metavar = 'File',required = False, help = 'The provided special gene list.Per gene for line or \',\' separated genes per line.')
#parser.add_argument('-PO','--PatientOnly',metavar = 'String', required = False, help = 'Show site in propositus only or no(show both sites in all members,eg:show site in normal sample but not in patient also).',choices = ['Y','N'],default = 'Y')
parser.add_argument('-M', '--hmode', metavar = 'String', required = False, help = 'If need check heredity mode,give the heretance mode,choices = [R,D], R for recessive, D for dominant',choices = ['R','D'])
parser.add_argument('-AF','--AF', metavar = 'Float', required = False, help = 'Min AF value to judge a confidence mutaion,default=0.25', default='0.25')
parser.add_argument('-soft','--soft',metavar = 'String',required = False, help = 'The software be used for snp/indel detecting.',choices = ['samtools','gatk'],default = 'samtools')
parser.add_argument('-D','--depth',metavar = 'Int', required = False, help = 'Min depth to judge a confidence mutaion,default=10',default='10')
parser.add_argument('-bp','--bp',metavar = 'Int', required = False, help = 'Range from indel to judge confidence indel, must smaller than 5. default=3',default='3')
argv = vars(parser.parse_args())

##get parameter
samin = argv['relation']
famid = argv['familyID']
suffix = argv['input']
outpath = os.path.split(os.path.abspath(argv['out']))[0]
prefix = os.path.split(os.path.abspath(argv['out']))[1]
#pocheck = argv['PatientOnly']
if argv['hmode']:
    hmode = argv['hmode']
af = float(argv['AF'])
udepth = int(argv['depth'])
bp = int(argv['bp'])
##check output path
#if not os.path.exists(out):
#    os.makedirs(out)
#out = os.path.abspath(out)

##get gene list and tableth
#if argv['genelist']:  ##genelist for 2 tables: [in gene list] and [not in gene list]
#    genef = argv['genelist']
#    geneli,sfugl,afugl,ugene = getgli(genef)


checkConfidence = argv['CheckConfidende']
if checkConfidence == 'Y' and (not argv['mpileup'] and not argv['mvcf']):
    exit('-Check parameter Y must use based on -mpileup or -mvcf')
if argv['mpileup']:
    mipleupf = os.path.abspath(argv['mpileup'])
if argv['mvcf']:
    mvcff = os.path.abspath(argv['mvcf'])

soft = argv['soft'].strip()

D = {}
count = 0
##read samp_info
##D {familyID:{sample1:relation|sex|type|chinesename,sample2:relation|sex|type|chinesename}}
#titleLi = ('#familyid','sampleid','name','sex','relation','normal/patient')
titleLi = ('#familyid','sampleid','sex','normal/patient')
titleLi2 = ('#familyid','sampleid','sex','n/p')
patientLi = []
normalLi = []
allsamLi = []
for n in (safe_open(samin,'r')):
    firstline = [i.lower() for i in n.strip().split('\t')]
    if n.startswith('#') and '#familyid' in firstline:
        #firstline = [i.lower() for i in n.strip().split('\t')]
        #if not set(titleLi).issubset(set(firstline)) and not set(titleLi2).issubset(set(firstline)):
        #    exit('Wrong sample info format: #FamilyID, SampleID, SEX, Relation and Normal/Patient(N/P) must in the firt row.Ignore case')
        fpos = getpos2('#FamilyID',firstline)
        spos = getpos2('SampleID',firstline)
        dpos = getpos2('Normal/Patient',firstline,'N/P')
        #rpos = getpos2('Relation',firstline)
        xpos = getpos2('SEX',firstline)
        #npos = getpos2('Name',firstline)
    elif n.startswith('#') and '#familyid' not in n:continue    
    else:
        nLine = n.strip().split('\t')
        if nLine[fpos] == famid:
            if nLine[spos] not in allsamLi:
                allsamLi.append(nLine[spos])
            if nLine[dpos].lower() != 'p':
                if nLine[spos] not in normalLi:
                    normalLi.append(nLine[spos])
            if nLine[dpos].lower() == 'p':
                if nLine[spos] not in patientLi:
                    patientLi.append(nLine[spos])
            if nLine[fpos] == famid and nLine[fpos] not in D:
                D[nLine[fpos]] = {}
            
        #distype[nLine[spos]] = nLine[dpos]
        #D[nLine[fpos]][nLine[spos]] = nLine[rpos]+'|'+nLine[xpos]+'|'+nLine[dpos]
            D[nLine[fpos]][nLine[spos]] = nLine[dpos]
print D

if checkConfidence == 'Y':
    mpileupDic = {}
    mposDic = {}        
    if argv['mpileup']:
        mpileupinfo = safe_open(mipleupf,'r')
        print 'Not patient:',normalLi
        print 'Patient:',patientLi
        mtype = '/'
        for line in mpileupinfo:
            if line.startswith('Mutation'):
                mpileuptitle = line.strip().split('\t')
                typos = getpos2('MutationType',mpileuptitle)
                mchrpos = getpos2('Chr',mpileuptitle)
                mpospos = getpos2('Position',mpileuptitle)
                mrefpos = getpos2('Ref',mpileuptitle)
                maltpos = getpos2('Alt',mpileuptitle)
                for nor in normalLi:
                    mposDic[nor] = [getpos(nor+'_ReadsNumber',mpileuptitle),getpos(nor+'_Base',mpileuptitle)]
                for pat in patientLi:
                    mposDic[pat] = [getpos(pat+'_ReadsNumber',mpileuptitle),getpos(pat+'_Base',mpileuptitle)]
                continue
            else:
                uline = line.strip().split('\t')
                altinfo = '_'.join([uline[mchrpos],uline[mpospos],uline[mrefpos],uline[maltpos]])
                ttype = uline[typos]
                ##mpileupDic for {sample:{alt:af,dp}},alt is chr_pos_re_alt
                for nor in normalLi:
                    if nor not in mpileupDic:
                        mpileupDic[nor] = {}
                    nmpileupdp = int(uline[mposDic[nor][0]])  ##total coverage
                    nmpileupaf = '/'   ## when no coverage or
                    if nmpileupdp != 0: ###count af when there are reads coverd
                        if ttype == 'snp':  ##对于snp,数有多少个突变
                            nvdp = len(re.findall(uline[maltpos].upper(),uline[mposDic[nor][1]])) + \
                                   len(re.findall(uline[maltpos].lower(),uline[mposDic[nor][1]]))  ##count mutaion number   
                            needrmLi = (re.findall(r'[\+\-](\d+)(\w+)',uline[mposDic[nor][1]]))  ##一行read同时有snp与indel时，snp数量可能多数，需要减去一行里面的Indel中的对应的snp突变
                            needrmcount = 0
                            for rmitem in needrmLi:
                                needrmcount += len(re.findall(uline[maltpos].upper(),rmitem[1][:int(rmitem[0])])) + \
                                len(re.findall(uline[maltpos].lower(),rmitem[1][:int(rmitem[0])]))
                            nvdp -= needrmcount   ##减去indel中的碱基突变
                            mtype = 'snp'
                        else: ##对于indel， 只数+-后面的突变有多少个，也就是有多少个加减，indel不一样的情况不考虑，有时候容易个别比对错误，只考虑个数
                            nvdp = len(re.findall(r'[\+\-](\d+)',uline[mposDic[nor][1]]))
                            mtype = 'indel'
                        #else:

                            #nvdp = len(re.findall(r'[ATGCatgc]',uline[mposDic[nor][1]]))
                            #notmut = re.findall(r'[\+\-](\d+)',uline[mposDic[nor][1]])  ##indel 突变不能计算多次dp，一个indel对dp的贡献只能是1，先减去indel的，再加回indel的个数
                            #for ct in notmut:
                            #    nvdp -= int(ct)
                            #if len(re.findall(r'[\+\-](\d+)',uline[mposDic[nor][1]])) >= 3: ##大于3个indel，认为是indel,如果碰见复杂的，又有snp又有indel的，这个地方可能不准确，但是对后面的突变的准确性的判断影响不会太大
                            #    nvdp += len(notmut) ##indel 要加回发生indel的个数
                            #    mtype = 'indel'
                        nmpileupaf = float(nvdp)/nmpileupdp
                    mpileupDic[nor][altinfo] = [nmpileupaf,nmpileupdp,mtype]
                for pat in patientLi:
                    if pat not in mpileupDic:
                        mpileupDic[pat] = {}
                    mpileupdp = int(uline[mposDic[pat][0]])
                    mpileupaf = '/'
                    if mpileupdp != 0:
                        if ttype == 'snp':  ##对于snp,数有多少个突变
                            vdp = len(re.findall(uline[maltpos],uline[mposDic[pat][1]])) + \
                            len(re.findall(uline[maltpos].lower(),uline[mposDic[pat][1]]))  ##count mutaion number   
                            pneedrmLi = (re.findall(r'[\+\-](\d+)(\w+)',uline[mposDic[pat][1]]))  ##一行read同时有snp与indel时，snp数量可能多数，需要减去一行里面的Indel中的对应的snp突变
                            pneedrmcount = 0
                            for prmitem in pneedrmLi:
                                pneedrmcount += len(re.findall(uline[maltpos].upper(),prmitem[1][:int(prmitem[0])])) + \
                                len(re.findall(uline[maltpos].lower(),prmitem[1][:int(prmitem[0])]))

                            vdp -= pneedrmcount   ##减去indel中的碱基突变
                            mtype = 'snp'
                        else: ##对于indel， 只数+-后面的突变有多少个，也就是有多少个加减
                            vdp = len(re.findall(r'[\+\-](\d+)',uline[mposDic[pat][1]]))
                            mtype = 'indel'
                        #else:
                        #    vdp = len(re.findall(r'[ATGCatgc]',uline[mposDic[pat][1]]))
                        #    notmut = re.findall(r'[\+\-](\d+)',uline[mposDic[pat][1]])  ##indel 突变不能计算多次dp，一个indel对dp的贡献只能是1，先减去indel的，再加回indel的个数
                        #    for ct in notmut:
                        #        vdp -= int(ct)
                        #    if len(re.findall(r'[\+\-](\d+)',uline[mposDic[pat][1]])) >= 3:
                        #        vdp += len(notmut) ##indel 要加回发生indel的个数
                        #        mtype = 'indel'
                        mpileupaf = float(vdp)/mpileupdp
                    mpileupDic[pat][altinfo] = [mpileupaf,mpileupdp,mtype]

        #print 'Mpileup count information:',mpileupDic        
    elif argv['mvcf']:
        mvcfinfo = safe_open(mvcff,'r')
        for line in mvcfinfo:
            if line.startswith('#CHROM'):
                title = line.strip().split('\t')
                chrpos = title.index('#CHROM')
                pospos = title.index('POS')
                refpos = title.index('REF')
                altpos = title.index('ALT')
                formatpos = title.index('FORMAT')
            else:
                if line.startswith('#'):continue
                uline = line.strip().split('\t')
                formatinfo = uline[formatpos].split(':')
                for sam in allsamLi:
                    if sam not in mpileupDic:
                        mpileupDic[sam] = {}
                    saminfo = uline[title.index(sam)].strip().split(':')
                    #if re.search(r'\d+',saminfo[0]):
                    dpdic = dict(zip(formatinfo,saminfo))
                    af = '/'
                    for alt in uline[altpos].strip().split(','):
                        altinfo = '_'.join([uline[chrpos],uline[pospos],uline[refpos],alt])
                        samdepth = int(dpdic['DP'])
                        if re.match(r'\d+',saminfo[0]):
                            if soft == 'samtools':
                                altdp = float(dpdic['DV'])
                            else:
                                for adnum in range(1,len(dpdic['AD'].strip().split(','))):
                                    altdp += float(dpdic['AD'].strip().split(',')[adnum])
                            af = altdp/samdepth
                        mpileupDic[sam][altinfo] = [af,samdepth]
    #print 'mpileupDic:', mpileupDic

for fid in D.keys():
    member = []
    othmember = []
    confiinfo = getinfo(suffix,fid,normalLi,patientLi,mpileupDic,checkConfidence,othmember)
    member.extend(othmember)
    with open(os.path.join(outpath,prefix+'.confidence.xls'),'w') as file:
        for n in confiinfo:
            file.write('\t'.join(n)+'\n')


