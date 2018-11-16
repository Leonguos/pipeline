#!/usr/bin/python
# -*- coding=utf-8 -*-
#==============================================================================================================================
# Name               : IntegrateFile_pep.v1.py
# Created On         : 2018-04-20
# Author             : lili
# Comment            : Integrate the results after pipeline.
# Last Modified By   : lili
# Last Modified On   : 
# Modified           :2018.7.2: modify the column of the file. Info,Format after the expreesion info using 4.6 annovar pipeline. See /NJPROJ2/DISEASE/WORK/lili/test/integrate/test1/Integrate.xls 
# Test Usage         :python IntegrateFile --pwd --infile --saminfo --analy_array --outpath
#==============================================================================================================================
import os
import sys
import argparse
import linecache
import pprint
from operator import itemgetter, attrgetter

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.dirname(os.path.dirname(BASE_DIR))
sys.path.append(ROOT_DIR)

from config import config
import utils


global ConfigDic
ConfigDic = {}

class File:
    def __init__(self):
        pass
    def safeOpen(self,fileName,mode):
        try:
            if not fileName.endswith('.gz'):
                return open(fileName,mode)
            else:
                import gzip
                return gzip.open(fileName,mode)
        except IOError:
            print fileName + ' do not exist!'
    def getpos2(self,name,title,name2=''):
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
    def getMpos(self,title):
        chrp = title.index('CHROM')
        posp = title.index('POS')
        refp = title.index('REF')
        altp = title.index('ALT')
        genep = title.index('GeneName')
        formatp = title.index('FORMAT')
        orirefp = title.index('Ori_REF')
        #4.5 no tissure info
        #pathp = title.index('REACTOME_PATHWAY')
        pathp = title.index('Subcellular_Location')
        return chrp,posp,refp,altp,genep,formatp,orirefp,pathp    
    def dealSaminfo(self,ProCheckInfo,ProFile,sameinfo,ANALY_DICT):
    #########################################################################################################################
    ##fidDic={fid1:{
    #		    'sidList':[sid1,sid2...],
    #               sid1:{'sex':sex1,'pn':pn1},
    #		    sid2:{'sex':sex2,'pn':pn2},
    #		    ...
    #		    'allfiles':[filesList1],
    #		    'tagfiles':[file1,file2,...],
    #               },
    #         fid2:{
    #		    'sidList':[sid1,sid2...],
    #               sid1:{'sex':sex1,'pn':pn1},
    #		    sid2:{'sex':sex2,'pn':pn2},
    #		    ...
    #		    'allfiles':[filesList],
    #               'tagfiles':[file1,file2,...],
    #               },
    #         ......
    #           }
    #tagDic = {file1:tag1, file2:tag2, file3:tag2, ...}
    ##fidList=[fid1,fid2,...]
    ###########################################################################################################################
        fidInfo = {}
        fidList = []
        tagDic = {}
        for line in self.safeOpen(sameinfo,'r'):
            tmpline = [i.lower() for i in line.strip().split('\t')]
            if '#familyid' in tmpline:
                head = line.strip().split('\t')
                pfid = self.getpos2('#FamilyID',head)
                psid = self.getpos2('SampleID',head)
                psex = self.getpos2('sex',head)
                pnp = self.getpos2('Normal/Patient',head)
                if ANALY_DICT['denovo_samtools']:
                    ppa = self.getpos2('Pa',head)
                    pma = self.getpos2('Ma',head)
            elif tmpline[0].startswith('#disease'):
                diseases = tmpline[0].split(':')[1]
                diseases = '"{}"'.format(diseases.strip().strip('"').strip("'"))
                fidInfo['diseases'] = diseases
            elif tmpline[0].startswith('#b'):
                pass
            else:
                uline = line.strip().split('\t')
                fid = uline[pfid]
                sid = uline[psid]
                sexinfo = uline[psex]
                pninfo = uline[pnp]
                if fid not in fidInfo:
                    fidInfo[fid] = {}
                    fidList.append(fid)
                    fidInfo[fid]['allfiles'] = []
                if not fidInfo[fid].get(sid,False):
                    fidInfo[fid].setdefault('sidList',[]).append(sid)#fidInfo[fid]['sidList'] contain all samples in fid
                    fidInfo[fid][sid] = {'sex':sexinfo,'pn':pninfo}
                    # print fidInfo[fid][sid]
                    # print uline, ppa
                    if ANALY_DICT['denovo_samtools'] and fidInfo[fid][sid]['pn'] == 'P':
                        if (len(uline) > ppa) and (len(uline) > pma):
                            fidInfo[fid][sid]['pa'] = uline[ppa]
                            fidInfo[fid][sid]['ma'] = uline[pma]
                        else:
                            fidInfo[fid][sid]['pa'] = ''
                            fidInfo[fid][sid]['ma'] = ''
                        #print fidInfo[fid]
                else:
                    print('Warnning:\tThe SampleID is dup.Please check:\tfid,\tsid')
                if ANALY_DICT['filter_db']:
                    filterDBsnpfile = os.path.join(ConfigDic['filterDBsnpDir'],'snp.merged.freq.func.syn.deleterious.xls')
                    tagDic[filterDBsnpfile] = 'FilterDB_SNP'
                    if filterDBsnpfile not in fidInfo[fid]['allfiles']:
                        fidInfo[fid]['allfiles'].append(filterDBsnpfile)
                    filterDBindelfile = os.path.join(ConfigDic['filterDBindelDir'],'indel.merged.freq.func.syn.deleterious.xls')
                    tagDic[filterDBindelfile] = 'FilterDB_InDel'
                    if filterDBindelfile not in fidInfo[fid]['allfiles']:
                        fidInfo[fid]['allfiles'].append(filterDBindelfile)
                if ANALY_DICT['model_dominant']:
                    modelDsnpfile = os.path.join(ConfigDic['modelFDir'],fid,fid+'.snp.dominant.xls')
                    tagDic[modelDsnpfile] = 'Dominant_SNP'
                    if modelDsnpfile not in fidInfo[fid]['allfiles']:
                        fidInfo[fid]['allfiles'].append(modelDsnpfile)
                    modelDindelfile = os.path.join(ConfigDic['modelFDir'],fid,fid+'.indel.dominant.xls')
                    tagDic[modelDindelfile] = 'Dominant_InDel'
                    if modelDindelfile not in fidInfo[fid]['allfiles']:
                        fidInfo[fid]['allfiles'].append(modelDindelfile)
                if ANALY_DICT['model_recessive']:
                    modelRsnpfile = os.path.join(ConfigDic['modelFDir'],fid,fid+'.snp.recessive.xls')
                    tagDic[modelRsnpfile] = 'Recessive_SNP'
                    if modelRsnpfile not in fidInfo[fid]['allfiles']:
                        fidInfo[fid]['allfiles'].append(modelRsnpfile)
                    modelRindelfile = os.path.join(ConfigDic['modelFDir'],fid,fid+'.indel.recessive.xls')
                    tagDic[modelRindelfile] = 'Recessive_InDel'
                    if modelRindelfile not in fidInfo[fid]['allfiles']:
                        fidInfo[fid]['allfiles'].append(modelRindelfile)
                    modelCfile = os.path.join(ConfigDic['modelFDir'],fid,fid+'.snp.indel.compound_heterozygous.xls')
                    tagDic[modelCfile] = 'Compound_Heterozygous_SNPInDel'
                    if modelCfile not in fidInfo[fid]['allfiles']:
                        fidInfo[fid]['allfiles'].append(modelCfile)

                if ANALY_DICT['denovo_samtools']:
                    #print "Yes"
                    # print "Denovo Intersect ..."
                    if fidInfo[fid][sid]['pn'] == 'P':
                        # print '>>>', sid, fidInfo[fid][sid]['pa'], fidInfo[fid][sid]['ma']
                        if (fidInfo[fid][sid]['pa'] != '') and (fidInfo[fid][sid]['ma'] != ''):
                            DenovoSNPfile = os.path.join(ConfigDic['denovoDir'], fid+'_'+sid, 'SNP', fid+'_'+sid+'.denovo.snp.intersect.xls')
                            tagDic[DenovoSNPfile] = 'Denovo_SNP_'+sid
                            if DenovoSNPfile not in fidInfo[fid]['allfiles']:
                                fidInfo[fid]['allfiles'].append(DenovoSNPfile)
                            DenovoINDELfile = os.path.join(ConfigDic['denovoDir'], fid+'_'+sid, 'INDEL', fid+'_'+sid+'.denovo.indel.intersect.xls')
                            tagDic[DenovoINDELfile] = 'Denovo_INDEL_'+sid
                            if DenovoINDELfile not in fidInfo[fid]['allfiles']:
                                fidInfo[fid]['allfiles'].append(DenovoINDELfile)
                                #print fidInfo[fid]['allfiles']
                if ANALY_DICT['filter_acmg']:
#                    acmgPfile = os.path.join(ConfigDic['acmgDir'],fid,fid+'.snp.indel.Pathogenic.xls')
                    acmgPfile = os.path.join(ConfigDic['acmgDir'],fid,fid,fid+'.snp.indel.Pathogenic.xls')
                    tagDic[acmgPfile] = 'ACMG_Pathogenic'
                    if acmgPfile not in fidInfo[fid]['allfiles']:
                        fidInfo[fid]['allfiles'].append(acmgPfile)
#                    acmgLPfile = os.path.join(ConfigDic['acmgDir'],fid,fid+'.snp.indel.LikelyPathogenic.xls')
                    acmgLPfile = os.path.join(ConfigDic['acmgDir'],fid,fid,fid+'.snp.indel.LikelyPathogenic.xls')
                    tagDic[acmgLPfile] = 'ACMG_LikelyPathogenic'
                    if acmgLPfile not in fidInfo[fid]['allfiles']:
                        fidInfo[fid]['allfiles'].append(acmgLPfile)

        for fid in fidList:
             for ifile in fidInfo[fid]['allfiles']:
                 done = ProCheckInfo.checkFile(ifile)
                 if (done != 0):
                     fidInfo[fid].setdefault('tagfiles',[]).append(ifile)
                     fidInfo[fid][ifile] = tagDic[ifile]
             if ANALY_DICT['phenolyzer']:
                 phenolyzerfile = os.path.join(ConfigDic['phenolyzerDir'],'AllGene_list.xls')
                 if not os.path.exists(phenolyzerfile):
                     phenolyzerfile = os.path.join(ConfigDic['integrateDir'],'AllGene_list.xls')
                 tagDic[phenolyzerfile] = 'Phenolyzer'
                 if phenolyzerfile not in fidInfo[fid]['allfiles']:
                     fidInfo[fid]['allfiles'].append(phenolyzerfile)
                     fidInfo[fid].setdefault('tagfiles',[]).append(phenolyzerfile)
                     fidInfo[fid][phenolyzerfile] = tagDic[phenolyzerfile]
                 

        print fidInfo
        print tagDic
        return fidInfo,fidList
########check info: analyarray,directory
class CheckInfo:

    def __int__(self):
        pass
        
    def checkDir(self,Dir,war):
        if not os.path.exists(Dir):
            if war == 'Y':
                print('Warning:\tThe directory is not exit.\t'+Dir+'\tPlease check it!')
            else:
                print('Note:\tThe directory is not exit.\t'+Dir+'\tAnd will make new.')
                os.makedirs(Dir)

    def checkFile(self,file1,done=1):
        if not os.path.exists(file1):
            exit('Error:\tThe file is not exit.\t'+file1)
            done = 0
        else:    
            done = 1
            print ('Check\t'+file1)
            aline = linecache.getline(file1,2)
            print aline
            if aline == '':
                done = 0
                print('Warnning :\tThere is only one title in the file.\t'+file1)

class GetInfo(object):

    def __int__(self):
        pass

    def getDir(self, ProCheckInfo, analydir, ANALY_DICT, newjob, analy_array):

        if max(analy_array) > 5 :
            mappingDir = os.path.join(analydir,'Mapping')
            ConfigDic['mappingDir'] = mappingDir
            advanceDir = os.path.join(analydir,'Advance',newjob)
            ProCheckInfo.checkDir(advanceDir,'Y') 
            ConfigDic['advanceDir'] = advanceDir
            integrateDir = os.path.join(advanceDir,'IntegrateResult')
            ConfigDic['integrateDir'] = integrateDir
            confidenceDir = os.path.join(advanceDir,'IntegrateResult/mkconfidence')
            ConfigDic['confidenceDir'] = confidenceDir
            if ANALY_DICT['filter_db']:
                filterDBsnpDir = os.path.join(advanceDir,'Merged_vcf/Filter/SNP')
                ProCheckInfo.checkDir(filterDBsnpDir,'Y')
                ConfigDic['filterDBsnpDir'] = filterDBsnpDir
                filterDBindelDir = os.path.join(advanceDir,'Merged_vcf/Filter/INDEL')
                ProCheckInfo.checkDir(filterDBindelDir,'Y')
                ConfigDic['filterDBindelDir'] = filterDBindelDir
            if ANALY_DICT['filter_model']:
                modelFDir = os.path.join(advanceDir,'ModelF')#+familyID
                ProCheckInfo.checkDir(modelFDir,'Y')
                ConfigDic['modelFDir'] = modelFDir
            if ANALY_DICT['filter_acmg']:
                ##FilterDB
                filterDBsnpDir = os.path.join(advanceDir,'Merged_vcf/Filter/SNP')
                ProCheckInfo.checkDir(filterDBsnpDir,'Y')
                ConfigDic['filterDBsnpDir'] = filterDBsnpDir
                filterDBindelDir = os.path.join(advanceDir,'Merged_vcf/Filter/INDEL')
                ProCheckInfo.checkDir(filterDBindelDir,'Y')
                ConfigDic['filterDBindelDir'] = filterDBindelDir
                ##ACMG
                acmgDir = os.path.join(advanceDir,'ACMG')#+familyID
                ProCheckInfo.checkDir(acmgDir,'Y')
                ConfigDic['acmgDir'] = acmgDir

            if ANALY_DICT['denovo']:
                denovoDir = os.path.join(advanceDir,'Denovo/Intersect')#DenovoSam,DenovoF,Intersect,  DenovoSam/fid_sid/fid_sid.*.xls
                ProCheckInfo.checkDir(denovoDir,'Y')
                ConfigDic['denovoDir'] = denovoDir
            if ANALY_DICT['phenolyzer']:
                phenolyzerDir = os.path.join(advanceDir,'Network')#AllGene_list.xls
                ProCheckInfo.checkDir(phenolyzerDir,'Y')
                ConfigDic['phenolyzerDir'] = phenolyzerDir
        return ConfigDic
    #return advanceDir,integrateDir,filterDBsnpDir,filterDBindelDir,modelFDir,acmgDir,denovoDir,phenolyzerDir

    def getsex(self,string):
        if string.lower() == 'm' or string.lower() == 'male':
            return '1'
        elif string.lower() == 'f' or string.lower() == 'female':
            return '2'
        else:
            return '0'

    def getnp(self,string):
        if string.lower() == 'n':
            return '1'
        elif string.lower() == 'p':
            return '2'
        else:
            return '0'
    
    def getSam(File,SamLi):
        for line in self.safeOpen(File,'r'):
            if line.startswith('Priority') or line.startswith('CHROM'):
                title = line.strip().split('\t')
                chrp,posp,refp,altp,formatp,orirefp,pathp = ProFile.getMpos(title)
                for i in range(formatp+1,orirefp):
                    sam = title[i]
                    if sam not in SamLi:
                        SamLi.append(sam)
        return SamLi
##handle file for fiuline,SamInfo(GT),filterpos
def DealFile(ProFile,File,Filterpos,SamLi,SamInfo,tagDic={}):#File:fidInfo   ,SamLi=fidDic[fid]['sidList'],  SamInfo=fidDic[fid][sid],  tagDic=fidDic[fid][file1]
    for line in ProFile.safeOpen(File,'r'):
        if line.startswith('Priority') or line.startswith('CHROM'):
            title = line.strip().split('\t')
            chrp,posp,refp,altp,genep,formatp,orirefp,pathp = ProFile.getMpos(title)
            fititle = title[:formatp-1]
            fititle.extend(title[orirefp:pathp+1])
            fititle.extend(title[formatp-1:formatp+1])
            fititle.extend(SamLi)
            fititle.extend(title[pathp+1:])
            #fititle.extend(title[orirefp:])
            continue
        else:
            uline = line.strip().split('\t')
            ##final uline
            fiuline = uline[:formatp-1]
            fiuline.extend(uline[orirefp:pathp+1])
            fiuline.extend(uline[formatp-1:formatp+1])
            for sam in SamLi:##所有的样本，在此文件中的信息，没有则标记成'/'，按照SamLi的顺序来。后面的文件的title，最好有前面的文件的样本，不然前面的样本的信息会被标记成'/'
                if sam not in SamInfo:
                    SamInfo[sam] = {}
            for altinfo in uline[altp].strip().split(','):
                if 'X' in uline[chrp]:
                    chr = 23
                elif 'Y' in uline[chrp]:
                    chr = 24
                else:
                    if 'chr' in uline[chrp]:
                        chr = int(uline[chrp][3:])
                    else:
                        chr = int(uline[chrp])
                pos = (chr,int(uline[posp]),uline[refp],altinfo)
                for sam in SamLi: ##获得每个文件中，对应的样本的突变信息，最开始初始化为'/'，读取到对应的文件，有此样本的信息的时候，样本的突变信息将会覆盖'/'
                    if pos not in SamInfo[sam]:
                        SamInfo[sam][pos] = '/'
                        if sam in title:
                            SamInfo[sam][pos] = uline[title.index(sam)]
                    else:
                        if sam in title:
                            SamInfo[sam][pos] = uline[title.index(sam)]
                    sinfo = SamInfo[sam][pos]
                    #print 'final sinfo: ',sam,pos,sinfo
                    fiuline.append(sinfo)
                #fiuline.extend(uline[orirefp:])
                fiuline.extend(uline[pathp+1:])
                if pos not in Filterpos:  ## for the first time apperence
                    if tagDic != {}:
                        givetag = tagDic[File]
                        #if 'novo' not in givetag.lower():
                        if chr == 23:
                            givetag = 'X link'
                        if chr == 24:
                            givetag = 'Y link'
                        Filterpos[pos] = [fiuline,1,givetag]
                    else:
                        Filterpos[pos] = [fiuline,1]
                else:
                    if tagDic != {}:
                        givetag = tagDic[File]
                        Filterpos[pos][0] = fiuline ##位点就算已经在FilterPos的字典中，这一行的信息fiuline还是要重新写入，因为样本的突变信息可能有变化（从'/'变成后来的有突变信息
                        if 'novo' not in givetag.lower():
                            if chr == 23:
                                givetag = 'X link'
                            if chr == 24:
                                givetag = 'Y link'
                        if givetag not in Filterpos[pos][2]:
                            Filterpos[pos][2] += ','+givetag
                        #if not FilterPos[pos][2].get(givetag,False):
                        #    if givefam not in FilterPos[pos][2][givetag]:
                        #        FilterPos[pos][2][givetag].append(givefam)
                        #    else:
                        #        print 'Please check the input file has been repeated: '+File+' of '+givetag+' from '+givefam
                    else:
                        Filterpos[pos][0] = fiuline
                    Filterpos[pos][1] += 1  ##位点已经存在，再次写入的时候，记录+1
            Filterpos[pos].append(uline[genep])#Filterpos[pos]:[fiuline,count,genename]
    return fititle,Filterpos,SamInfo

def addTag(ProFile,fidInfo,fidList,samp_info,configDic,confidence='N'): 
    ##files=fidDic[fid1]['files']
    if fidList == []:
        exit('-sam_info must be gived, and famliyid must be not empty.')
    ##get positon and line for all files of one family
    genelist = []
    for sfid in fidList:
        FilterPos = {}
        SamLi = []
        SamInfo = {}
        SamLi = fidInfo[sfid]['sidList']
        files = fidInfo[sfid].get('tagfiles')
        tagDic = fidInfo[sfid]
        outpath = configDic['integrateDir']
        confidencePath = configDic['confidenceDir']
        with open(outpath+'/'+sfid+'.integrate.xls','w') as outfid:
            outtitle = []
            if files:
                if 'Phenolyzer' in tagDic.itervalues():
                #if 'AllGene_list.xls' in files:
                    outtitle.extend(['GeneCandidate','FamilyID','Tag','PhenolyzerScore','PhenolyzerGene'])
                else:
                    outtitle.extend(['GeneCandidate','FamilyID','Tag'])
            ##title : 
            dicGene = {}
            title = []
            if not files:
                continue
            for file in files:
                if 'AllGene_list.xls' in file:##get phenolyzer dict for gene
                    for iline in ProFile.safeOpen(file,'r'):
                        if iline.startswith('Rank'):
                            iiline = iline.strip().split('\t')
                            genep = iiline.index('Gene')
                            scorep = iiline.index('Score')
                            statusp = iiline.index('Status')
                        else:
                            iiline = iline.strip().split('\t')
                            gene = iiline[genep]
                            score = iiline[scorep]
                            status = iiline[statusp]
                            if not dicGene.get(gene,False):
                                dicGene.setdefault(gene,[]).append(score)
                                dicGene.setdefault(gene,[]).append(status)
                else:##other files for pos-line-tag
                    titletmp,FilterPos,SamInfo = DealFile(ProFile,file,FilterPos,SamLi,SamInfo,tagDic)
                    #title from the longest title of the files:
                    if set(title).issubset(set(titletmp)):
                        title = titletmp
                   
            outtitle.extend(title)
            outfid.write('\t'.join(outtitle)+'\n')
            #FilterPos[pos] = [fiuline,1,givetag,gene]
            ##write union,sort by chr,then by position
            for key in sorted(FilterPos.keys(),key=itemgetter(0,1)):
                gene = FilterPos[key][3]
                if dicGene.get(gene,False):
                    phenolyzer = dicGene[gene]##[score,status]
                else:
                    phenolyzer = ['--','--']
                tag = FilterPos[key][2]
                line = FilterPos[key][0]
                finalline = []
                finalline.append(gene)
                finalline.append(sfid)
                finalline.append(tag)
                if 'Phenolyzer' in tagDic.itervalues():
                    finalline.extend(phenolyzer)
                finalline.extend(FilterPos[key][0])
                outfid.write('\t'.join(finalline)+'\n')
                ##candidategene:
                for igene in gene.split(','):
                    if igene not in genelist:
                        genelist.append(gene)
        #make confidence:    
        #if confidence == 'Y':
        bamlist = [configDic['mappingDir']+'/'+i+'.'+i+'/'+i+'.final.bam' for i in SamLi]
        step1 = 'cd {path} \n##step1: generate pileup file\npython {BASE_DIR}/mpileup.py -B {bamfiles} -input {infile} -output {outfile} -O {outshell}'.format(
            path=outpath, bamfiles=','.join(bamlist), infile = outpath+'/'+sfid+'.integrate.xls', outfile=confidencePath+'/'+sfid+'.Mpileup.xls',
            outshell = confidencePath+'/'+sfid+'.mpileup.sh', BASE_DIR=BASE_DIR)
        step2 = '##step2: run\nsh {outshell}'.format(outshell=confidencePath+'/'+sfid+'.mpileup.sh')
        step3 = '##step3: mark confidence\npython {BASE_DIR}/MarkConfidence.py -I {infile} -O {outfile}  -Pileup {pileupfile} -R {sampleinfo} -F {fid} -Check Y'.format(
            infile=outpath+'/'+sfid+'.integrate.xls', outfile=confidencePath+'/'+sfid+'.integrate',
            pileupfile=confidencePath+'/'+sfid+'.Mpileup.xls', sampleinfo=samp_info, fid=sfid, BASE_DIR=BASE_DIR)
        with open(outpath+'/'+sfid+'_mkconfidence.sh','w') as mkconfidenceShell:
            #mkconfidenceShell.write(' && \\\n'.join([step1,step2,step3])+'\n')
            mkconfidenceShell.write(' && \\\n'.join([step1,step3])+'\n')
        if confidence == 'Y' and len(SamLi)!=1:
            print("make confidence:...")
            cmd = 'cd {path} \nsh {workshell}'.format(path=outpath,workshell=outpath+'/'+sfid+'_mkconfidence.sh')
            os.system(cmd)
        else:
            print("Warnning: You did not make confidence or the family just has only one sample.\nIt is suggested that confidence markers should be according to the number of variants.")
    with open(outpath+'/total.candidate.gene.xls','w') as candidategene:
        candidategene.write('\n'.join(genelist))

def catALL(fidList,analypath):
    cmd = 'rm -f {outfile}'.format(outfile=analypath+'/'+'Integrate.xls')
    os.system(cmd)
    #fidfile = analypath+'/'+sfid+'.integrate.xls'
    outfile = analypath+'/'+'Integrate.xls'
    #fidconfidencefile = analypath+'/'+'mkconfidence/'+sfid+'.integrate.confidence.xls'
    confidencefile = analypath+'/'+'Integrate.confidence.xls'
    if len(fidList) !=1:
        for sfid in fidList:
            catfidfilecmd = 'cat {fidfile} >> {outfile} &&\\\nrm -f {fidfile} '.format(fidfile=analypath+'/'+sfid+'.integrate.xls',outfile=outfile)
            os.system(catfidfilecmd)
            if os.path.exists(analypath+'/'+'mkconfidence/'+sfid+'.integrate.confidence.xls'):
                catconfidencecmd = '&&\\\ncat {fidconfidence} >> {confidencefile} && \\\nrm -f {fidconfidence}'.format(fidconfidence=analypath+'/'+'mkconfidence/'+sfid+'.integrate.confidence.xls', confidencefile=confidencefile)
                os.system(catconfidencecmd)
    else:
        cmd = 'mv -f {fidfile} {outfile}'.format(fidfile=analypath+'/'+fidList[0]+'.integrate.xls',outfile=outfile)
        os.system(cmd)
        if os.path.exists(analypath+'/'+'mkconfidence/'+fidList[0]+'.integrate.confidence.xls'):
            lnconfidencecmd = 'mv -f {fidconfidence} {confidencefile}'.format(fidconfidence=analypath+'/'+'mkconfidence/'+fidList[0]+'.integrate.confidence.xls', confidencefile=confidencefile)
            os.system(lnconfidencecmd)

#result just for family filter, excluding ACMG
def ecludeACMG(analypath):
    allcatfile = analypath+'/'+'Integrate.xls'
    extractfile = analypath+'/'+'Integrate.candidate.xls'
    allconfidencefile = analypath+'/'+'Integrate.confidence.xls'
    extractconfidencefile = analypath+'/'+'Integrate.confidence.candidate.xls'
    cmd = 'grep -v ACMG_ {catfile} > {outfile} '.format(catfile=allcatfile, outfile=extractfile)
    os.system(cmd)
    if os.path.exists(allconfidencefile):
        cmd1 = 'grep -v ACMG_ {allconfidence} > {Hconfidence}'.format(allconfidence=allconfidencefile, Hconfidence=extractconfidencefile)
        os.system(cmd1)

#def txt2excel():


def main():

    confidence = args['confidence']
    moduledir = args['moduledir']
    newjob = args['newjob']
    analydir = args['analydir']
    analy_array = args['analy_array'].split(',')

    ANALY_DICT = utils.get_analysis_dict(analy_array, config.ANALYSIS_CODE)
    pprint.pprint(ANALY_DICT)

    ConfigDic = {}

    info = GetInfo()

    ProCheckInfo = CheckInfo()
    ProFile = File()

    ConfigDic = info.getDir(ProCheckInfo, analydir, ANALY_DICT, newjob, analy_array)
    if args['out']:
        integrateDir = os.path.abspath(args['out'].strip())
        confidenceDir = os.path.join(integrateDir,'mkconfidence')
        ConfigDic['integrateDir'] = integrateDir
        ConfigDic['confidenceDir'] = confidenceDir
    print ConfigDic
    ProCheckInfo.checkDir(ConfigDic['confidenceDir'],'N')
    print ConfigDic
    ##########################################depend on ANALY_DICT, get sample info.#########################################
    #fidInfo:
    ########
    samp_info = os.path.abspath(args['samp_info'].strip())
    if os.path.isfile(samp_info):
        fidList = []
        fidInfo = {}
        fidInfo,fidList = ProFile.dealSaminfo(ProCheckInfo,ProFile,samp_info,ANALY_DICT)
    else:
        sys.exit('The '+samp_info+' is not a file.')

    #phenolyzer for all.gene.list.xls 
    if ANALY_DICT['phenolyzer'] and not os.path.exists(os.path.join(ConfigDic['phenolyzerDir'],'AllGene_list.xls')):
        cmd = 'python '+moduledir+'/Phenolyzer/phenolyzer-0.1.5/phenolyzer_pipe4.6.py --disease '+fidInfo['diseases']+' --dir '+ConfigDic['integrateDir']
        print cmd
        os.system(cmd)

    #######################################################addfile for every family###############################################
    addTag(ProFile,fidInfo,fidList,samp_info,ConfigDic,confidence)
    ##cat all family result:
    catALL(fidList,ConfigDic['integrateDir'])
    ##extract no ACMG for candidate:
    ecludeACMG(ConfigDic['integrateDir'])
    if ANALY_DICT['phenolyzer'] and not os.path.exists(os.path.join(ConfigDic['phenolyzerDir'],'AllGene_list.xls')):
        os.system('cd '+ConfigDic['integrateDir']+'; rm -rf AllGene_list.xls AllGene_score.xls phenolyzer.tar html network.js out*')


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description="Integrate all result based on the sampleinfo and analyarry",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog='Contact: lili@novogene.com')

    parser.add_argument('--analydir', help="the analysis directory")
    parser.add_argument('--out', help='the output directory', default='.')
    parser.add_argument('--samp-info', help="the sample_info file")
    parser.add_argument('--newjob', help="the job name")

    parser.add_argument(
        '--confidence',
        help="whether do confidence analysis or not",
        choices=['Y', 'N'],
        default='N')

    parser.add_argument(
        '--moduledir',
        help="The module directory[default=%(default)s]",
        default=config.CONFIG.get('software', 'moduledir'))
        
    parser.add_argument('--analy-array', help="the analysis code")

    args = vars(parser.parse_args())

    main()




