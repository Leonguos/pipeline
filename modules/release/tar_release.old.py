#!/usr/bin/env python
# -*- coding: utf-8 -*-
#change tar_release_v4.6_NJ_20180325.py to tar_release_v4.6_NJ_20180423.py
# Author: shexuan
# Date: 2017.9.17
# Version: 0.3
# Function : yourself look
# Add function: release several jobs together
# Add function: release files under */Result/job/ named by yourself

from __future__ import print_function, division
import os
import commands
import time
from datetime import datetime
# import gzip
import argparse
from argparse import RawTextHelpFormatter
# import smtplib
# from email.mime.text import MIMEText
# from email.header import Header
# from email.mime.multipart import MIMEMultipart
# from email.mime.application import MIMEApplication
#add separate to tar package separately by 20180514
import sys
# Solve the encoding problem
reload(sys)
sys.setdefaultencoding('utf-8')


# just suit for /PROJ/HUMAN/share/Disease_pipeline/Human_reseq/Human_reseq_pipeline.py pipeline.
# from django.template import Context, Template, loader
# from django.conf import settings
# add data_release by tar file less than 55G __date__20170707


class Release(object):
    '''Relsease whatever you want! :D '''

    def __init__(self):
        self.pn = list(os.popen("awk '{{print $1}}' {0}".format(args["pn"])))[
            0].strip()
        self.proj_name = list(os.popen("awk '{{print $2}}' {0}".format(args["pn"])))[
            0].strip()
        self.proj_dir = args["projdir"].strip()
        # modify by yincan 20171120
        # job_dir means */Result/*job_dir where exists data that we want to release
        self.job_dir = self.proj_dir + '/Result/' + args["date"].strip()
        # outdir means release_dir where we release results
        self.outdir = args["odir"].strip()
        self.now = datetime.now().strftime("%Y%m%d")
        self.analy_array = (args["analy_array"].strip()).split(",")
        self.specified_release = (args["release"].strip()).split(",")
        self.qc_list = self.proj_dir + "/" + "qc_list_" + args["pre"].strip()
        self.multi_checklist = self.proj_dir + '/Result/checkSize.xls'
        # self.multi_checklist = self.outdir + '/checkSize.xls'
        self.job_name = args["date"].strip().split(',')
        self.all_data = ["RawData", "QC", "Mapping",
                         "Variation", "Advance", "Readme"]
        self.mail = args["mail"].strip()
        self.separate = args["separate"].strip()
        self.yymail = args["yymail"].strip()
        self.checkpy = "/NJPROJ1/DISEASE/share/Disease/QC/CheckResult_v1.0.mod.20171127.py"
        self.force = args["force"].strip()
        self.checkpy_res = ''
        ############################################################
        # data_size < 55G, just tar all data
        # if release rawdata, the tar file named data_PN_date.tar
        self.raw_tardir = self.outdir + "/data_" + self.pn + "_" + args["pre"].strip()+"_"+self.now
        ##add 20180512
        self.rawtarpackage = "data_" + self.pn + "_" + args["pre"].strip()+"_"+self.now + ".tar"
        self.varadvpackage="analysis_" + self.pn + "_" + args["pre"].strip()+"_"+self.now + ".tar"
        self.cleanpackage="cleandata_" + self.pn + "_" + args["pre"].strip()+"_"+self.now + ".tar"
        self.bampackage="bamdata_" + self.pn + "_" + args["pre"].strip()+"_"+self.now + ".tar"
        ##add 20180512
        ##modify 20180412 add fenqi
        self.raw_release = self.raw_tardir + "/data_" + self.pn + "_" + args["pre"].strip()+"_"+self.now + ".tar"
        self.raw_tar_checklist = self.raw_tardir + '/checkSize.xls'
        # if not release rawdata, the tar file named analysis_PN_date.tar
        self.tardir = self.outdir + "/analysis_" + self.pn + "_" +args["pre"].strip()+"_"+ self.now
        self.release = self.tardir + "/analysis_" + self.pn + "_" +args["pre"].strip()+"_"+ self.now + ".tar"
        self.tar_checklist = self.tardir + '/checkSize.xls'
        ############################################################
        # data_size >= 55G, make a directory and link the file to the directory
        # if release rawdata, the dir named data_pn_date
        self.raw_release_dir = self.outdir + "/data_" + self.pn +"_"+args["pre"].strip()+ "_" +self.now
        self.not_tar_raw_checklist = self.raw_release_dir + '/checkSize.xls'
        # if not release rawdata, the dir named analysis_pn_date
        self.release_dir = self.outdir + "/" + "analysis_" + self.pn + "_" +args["pre"].strip()+ "_"+self.now
        self.not_tar_checklist = self.release_dir + '/checkSize.xls'

    def releaseType(self):
        '''if --release argument are not assigned, defaultRelease.'''
        analy_module = [int(float(i)) for i in self.analy_array]
        # specifify release data
        if self.specified_release != [""]:
            data_release = self.specified_release
            return data_release
        elif self.specified_release == [""]:
            # just QC
            if len(analy_module) == 1:
                data_release = ["RawData"]
                return data_release
            # primary analysis
            if 1 < len(analy_module) < 6:
                data_release = ["RawData", "Mapping", "Variation", "Readme"]
                return data_release
            # advance analysis
            else:
                data_release = ["RawData", "Mapping", "Variation",
                                "Advance", "Readme"]
                return data_release
    def sampleSum(self):
        
        data_release=self.releaseType()
        if "RawData" in data_release:
            sum= commands.getoutput(
                    ''' ls -l {0} | grep "^d" | wc -l '''.format(Release.outdir +
                        "/ReleaseResult/Data/RawData" ))
           
            self.rawtarpackage = "data_" + self.pn + "_" + args["pre"].strip()+"_"+"S"+sum+"_"+self.now + ".tar"##
            print(self.rawtarpackage)
           
        if "QC" in releasetype:
        
            sum= commands.getoutput(
                    ''' ls -l {0} | grep "^d" | wc -l '''.format(Release.outdir +
                        "/ReleaseResult/Data/CleanData" ))
           
            
            self.cleanpackage="cleandata_" + self.pn + "_" + args["pre"].strip()+"_"+"S"+sum+"_"+self.now + ".tar"
            print(self.cleanpackage)
        if "Mapping" in releasetype:

            sum= commands.getoutput(
                    ''' ls -l {0} | grep "^d" | wc -l '''.format(Release.outdir +
                        "/ReleaseResult/Data/BamData" ))
           
            self.bampackage="bamdata_" + self.pn + "_" + args["pre"].strip()+"_"+"S"+sum+"_"+self.now + ".tar"
            print(self.bampackage)
        else:
            pass
        

    def releaseSum(self):
        
        data_release=self.releaseType()
        if "RawData" in data_release:
            sum= commands.getoutput(
                    ''' ls -l {0} | grep "^d" | wc -l '''.format(Release.outdir +
                        "/ReleaseResult/Data/RawData" ))
            self.raw_tardir = self.outdir + "/data_" + self.pn + "_" + args["pre"].strip()+"_"+"S"+sum+"_"+self.now
            print(self.raw_tardir)
            self.raw_release_dir=self.outdir + "/data_" + self.pn + "_" + args["pre"].strip()+"_"+"S"+sum+"_"+self.now
            self.raw_release=self.raw_tardir+"/data_" + self.pn + "_" + args["pre"].strip()+"_"+"S"+sum+"_"+self.now+".tar"
           
            self.not_tar_raw_checklist = self.raw_release_dir + '/checkSize.xls'
            self.raw_tar_checklist = self.raw_tardir + '/checkSize.xls'
        elif "RawData"  not in releasetype and "QC" in releasetype:
        
            sum= commands.getoutput(
                    ''' ls -l {0} | grep "^d" | wc -l '''.format(Release.outdir +
                        "/ReleaseResult/Data/CleanData" ))
            self.tardir=self.outdir + "/analysis_" + self.pn + "_" + args["pre"].strip()+"_"+"S"+sum+"_"+self.now
            self.tar_checklist = self.tardir + '/checkSize.xls'
            print(self.tardir)
            self.release=self.tardir + "/analysis_" + self.pn + "_" + args["pre"].strip()+"_"+"S"+sum+"_"+self.now+".tar"
            
            self.release_dir=self.outdir + "/analysis_" + self.pn + "_" + args["pre"].strip()+"_"+"S"+sum+"_"+self.now
            self.not_tar_checklist = self.release_dir + '/checkSize.xls'
        elif "RawData"  not in releasetype and "Mapping" in releasetype:

            sum= commands.getoutput(
                    ''' ls -l {0} | grep "^d" | wc -l '''.format(Release.outdir +
                        "/ReleaseResult/Data/BamData" ))
            self.tardir=self.outdir + "/analysis_" + self.pn + "_" + args["pre"].strip()+"_"+"S"+sum+"_"+self.now
            self.tar_checklist = self.tardir + '/checkSize.xls'
            print(self.tardir)
            self.release=self.tardir + "/analysis_" + self.pn + "_" + args["pre"].strip()+"_"+"S"+sum+"_"+self.now+".tar"
            self.release_dir=self.outdir + "/analysis_" + self.pn + "_" + args["pre"].strip()+"_"+"S"+sum+"_"+self.now
            self.not_tar_checklist = self.release_dir + '/checkSize.xls'
        else:
            pass
            
    def argCheck(self):
        '''check if Invalid Arguments are input.'''
        # check args for --release
        data_release = self.releaseType()
        data_release_dict = self.data_release_tran(data_release)
        if len(self.job_name) == 1:
            data_dirs = os.listdir(self.job_dir)
            for data in data_release_dict:
                if self.specified_release != [""]:
                    if data_release_dict[data] not in data_dirs:
                        sys.exit('Invalid Arguments: ' + data)
                else:
                    if data_release_dict[data] not in data_dirs:
                        sys.exit(
                            'Default release failed because the Defualt release file {} dose Not exist, please use "--release" to release data.'.format(data))
        # check args for --date
            job_dirs = os.listdir(self.proj_dir + '/Result')
            for job in self.job_name:
                if job not in job_dirs:
                    sys.exit('Invalid job name: ' + job)
            return None

    def dataSizeCheck(self):
        '''checkout data size whether less than 55G to decide tar or copy on disk.'''
        ######### method 1, use os.path.getsize() ##########################################
        # data_release = self.releaseType()
        # all_files = []
        # for data in data_release:
        #     for root, dirs, files in os.walk(data):
        #         all_files.extend([os.path.join(root, file) for file in files])
        # size = sum([os.path.getsize(file) for file in all_files])
        # return size / (1024**3)
        ######### method 2, call linux command du -shL #####################################
        data_release = self.releaseType()
        data_release_dict = self.data_release_tran(data_release)
        data_dirs = [(self.job_dir + "/" + dir_) for dir_ in data_release_dict.values()]
        size = 0
        for data in data_dirs:
            res = os.popen("du -smL --apparent-size {0}".format(data))
            size += float((list(res)[0]).split()[0])
        return size / 1024, data_release

    def checkRelease(self):
        '''execute the check scripts.'''
        for job in self.job_name:
            res = os.popen(
                "python {0} -P {1} -J {2} -L {3}".format(self.checkpy, self.proj_dir, job, self.qc_list))
            self.checkpy_res = "".join(res)
        return self.checkpy_res


    def data_release_tran(self, data_release):
        dict_template = {"RawData": "ReleaseResult/Data/RawData",
                        "QC": "ReleaseResult/Data/CleanData",
                        "Mapping": "ReleaseResult/Data/BamData",
                        "Variation": "ReleaseResult/PrimaryAnalysis",
                        "Advance": "ReleaseResult/FinalResult",
                        "Readme": "ReleaseResult/Readme"}
        data_release_dict = {i:dict_template[i] for i in data_release}
        return data_release_dict
    def dataRelease(self):
        '''if data_release < 55G, tar all file, else cp on disk.'''
        self.sampleSum()
        self.releaseSum()
        data_size, data_release = self.dataSizeCheck()
        self.releaseSum()
        data_release_dict = self.data_release_tran(data_release)
        self.releaseSum()
        if os.getcwd() != self.job_dir:
            os.chdir(self.job_dir)
        else:
            pass
        # whether force tar/not_tar
        if self.force:
            if self.force == "tar":
                if "RawData" in data_release:
                    os.system("rm -rf {0}".format(self.raw_release_dir))
                    os.mkdir(self.raw_release_dir)
                    os.system("tar cphf {0} {1}".format(
                        self.raw_release, " ".join(data_release_dict.values())))
                    # os.system('''less {0} | awk '{{print $3"\t"$6}}' > {1} && unix2dos {1}'''.format(
                    #     self.raw_release, self.raw_tar_checklist))
                    os.system("perl /NJPROJ1/DISEASE/share/Software/bin/dirCheckSize.pl {0} -s {1}".format(
                        self.raw_tardir, self.raw_tar_checklist))
                    return self.raw_tardir, self.raw_tar_checklist
                else:
                    os.system("rm -rf {0}".format(self.release_dir))
                    os.mkdir(self.release_dir)
                    os.system("tar cphf {0} {1}".format(
                        self.release, " ".join(data_release_dict.values())))
                    # os.system('''less {0} | awk '{{print $3"\t"$6}}' > {1} && unix2dos {1}'''.format(
                    #     self.release, self.tar_checklist))
                    os.system("perl /NJPROJ1/DISEASE/share/Software/bin/dirCheckSize.pl {0} -s {1}".format(
                        self.tardir, self.tar_checklist))
                    return self.tardir, self.tar_checklist
            else:  # force not_tar
                if "RawData" in data_release:
                    # because copy just add file to the dir instead of overwrite,
                    # here remove the existed dir and recreated, and the following is the same
                    # tar small files into a tar
                    os.system("rm -rf {}".format(self.raw_release_dir))
                    os.mkdir(self.raw_release_dir)
                    for i in ['RawData', 'QC', 'Mapping']:
                        if i in data_release_dict:
                            os.system(
                                "cp -a {0} {1}".format(data_release_dict[i], self.raw_release_dir))
                            data_release_dict.pop(i)
                    os.system("tar cphf {0} {1}".format(
                        self.raw_release, " ".join(data_release_dict.values())))
                    if len(self.job_name) == 1:
                        # dirCheckSize.pl
                        os.system("perl /NJPROJ1/DISEASE/share/Software/bin/dirCheckSize.pl {0} -s {1}".format(
                            self.raw_release_dir, self.not_tar_raw_checklist))
                        return self.raw_release_dir, self.not_tar_raw_checklist
                    else:
                        return self.raw_release_dir, None
                else:  # not release RawData
                    os.system("rm -rf {0}".format(self.release_dir))
                    os.mkdir(self.release_dir)
                    for i in ['RawData', 'QC', 'Mapping']:
                        if i in data_release_dict:
                            os.system(
                                "cp -a {0} {1}".format(data_release_dict[i], self.release_dir))
                            data_release_dict.pop(i)
                    os.system("tar cphf {0} {1}".format(
                        self.release, " ".join(data_release_dict.values())))
                    if len(self.job_name) == 1:
                        os.system("perl /NJPROJ1/DISEASE/share/Software/bin/dirCheckSize.pl {0} -s {1}".format(
                            self.release_dir, self.not_tar_checklist))
                        return self.release_dir, self.not_tar_checklist
                    else:
                        return self.release_dir, None
        # tar/not_tar rely on data size (55G)
        #modify 20180512
        #print(data_size)
        if data_size < 50 and self.separate=="N": # tar
            if "RawData" in data_release:
                os.system("rm -rf {0}".format(self.raw_release_dir))
                os.mkdir(self.raw_release_dir)
                os.system("tar cphf {0} {1}".format(
                    self.raw_release, " ".join(data_release_dict.values())))
                # os.system('''less {0} | awk '{{print $3"\t"$6}}' > {1} && unix2dos {1}'''.format(
                #     self.raw_release, self.raw_tar_checklist))
                os.system("perl /NJPROJ1/DISEASE/share/Software/bin/dirCheckSize.pl {0} -s {1}".format(
                    self.raw_tardir, self.raw_tar_checklist))
                return self.raw_tardir, self.raw_tar_checklist
            else:
                os.system("rm -rf {0}".format(self.release_dir))
                os.mkdir(self.release_dir)
                os.system("tar cphf {0} {1}".format(
                    self.release, " ".join(data_release_dict.values())))
                # os.system('''less {0} | awk '{{print $3"\t"$6}}' > {1} && unix2dos {1}'''.format(
                #     self.release, self.tar_checklist))
                os.system("perl /NJPROJ1/DISEASE/share/Software/bin/dirCheckSize.pl {0} -s {1}".format(
                    self.tardir, self.tar_checklist))
                return self.tardir, self.tar_checklist
        elif data_size < 50 and self.separate=="Y":# tar separately modify 20180512
            
            if "RawData" in data_release:
                print(self.raw_release_dir)
                os.system("rm -rf {0}".format(self.raw_release_dir))
                os.mkdir(self.raw_release_dir)
                os.system("tar cphf {0} {1}".format(self.raw_release_dir+"/"+self.rawtarpackage,args["odir"]+"/ReleaseResult/Data/RawData"))
                if "QC" in data_release:
                    os.system("tar cphf {0} {1}".format(self.raw_release_dir+"/"+self.cleanpackage,args["odir"]+"/ReleaseResult/Data/CleanData"))
                if "Mapping" in data_release:
                    os.system("tar cphf {0} {1}".format(self.raw_release_dir+"/"+self.bampackage,args["odir"]+"/ReleaseResult/Data/BamData"))
                for i in ['RawData', 'QC', 'Mapping']:
                    if i in data_release:
                        data_release_dict.pop(i)
                        os.system("tar cphf {0} {1}".format(
                            self.raw_release_dir+"/"+self.varadvpackage," ".join(data_release_dict.values())))

                # os.system('''less {0} | awk '{{print $3"\t"$6}}' > {1} && unix2dos {1}'''.format(
                #     self.raw_release, self.raw_tar_checklist))
                print(self.raw_release_dir)
                os.system("perl /NJPROJ1/DISEASE/share/Software/bin/dirCheckSize.pl {0} -s {1}".format(
                    self.raw_release_dir, self.raw_tar_checklist))
                return self.raw_tardir, self.raw_tar_checklist
            else:
                os.system("rm -rf {0}".format(self.release_dir))
                os.mkdir(self.release_dir)
                if "QC" in data_release:
                    os.system("tar cphf {0} {1}".format(self.release_dir+"/"+self.cleanpackage,args["odir"]+"/ReleaseResult/Data/CleanData"))
                if "Mapping" in data_release:
                    os.system("tar cphf {0} {1}".format(self.release_dir+"/"+self.bampackage,args["odir"]+"/ReleaseResult/Data/BamData"))
                for i in ['RawData', 'QC', 'Mapping']:
                    if i in data_release:
                        data_release_dict.pop(i)
                        os.system("tar cphf {0} {1}".format(
                            self.release_dir+"/"+self.varadvpackage," ".join(data_release_dict.values())))
                

                
                # os.system('''less {0} | awk '{{print $3"\t"$6}}' > {1} && unix2dos {1}'''.format(
                #     self.release, self.tar_checklist))
                os.system("perl /NJPROJ1/DISEASE/share/Software/bin/dirCheckSize.pl {0} -s {1}".format(
                    self.tardir, self.tar_checklist))
                return self.tardir, self.tar_checklist
        if data_size >= 50:
                if "RawData" in data_release:
                    # because copy just add file to the dir instead of overwrite,
                    # here remove the existed dir and recreated, and the following is the same
                    # tar small files into a tar
                    os.system("rm -rf {}".format(self.raw_release_dir))
                    os.mkdir(self.raw_release_dir)
                    for i in ['RawData', 'QC', 'Mapping']:
                        if i in data_release:
                            os.system(
                                "cp -a {0} {1}".format(data_release_dict[i], self.raw_release_dir))
                            data_release_dict.pop(i)
                    os.system("tar cphf {0} {1}".format(
                        self.raw_release, " ".join(data_release_dict.values())))
                    if len(self.job_name) == 1:
                        # dirCheckSize.pl
                        os.system("perl /NJPROJ1/DISEASE/share/Software/bin/dirCheckSize.pl {0} -s {1}".format(
                            self.raw_release_dir, self.not_tar_raw_checklist))
                        return self.raw_release_dir, self.not_tar_raw_checklist
                    else:
                        return self.raw_release_dir, None
                else:  # not release RawData
                    os.system("rm -rf {0}".format(self.release_dir))
                    os.mkdir(self.release_dir)
                    for i in ['RawData', 'QC', 'Mapping']:
                        if i in data_release:
                            os.system(
                                "cp -a {0} {1}".format(data_release_dict[i], self.release_dir))
                            data_release_dict.pop(i)
                    os.system("tar cphf {0} {1}".format(
                        self.release, " ".join(data_release_dict.values())))
                    if len(self.job_name) == 1:
                        os.system("perl /NJPROJ1/DISEASE/share/Software/bin/dirCheckSize.pl {0} -s {1}".format(
                            self.release_dir, self.not_tar_checklist))
                        return self.release_dir, self.not_tar_checklist
                    else:
                        return self.release_dir, None

    def rawRelease(self, multi_release):
        '''For multiRelease, if any job release RawData, named data_pn_date, else analysis_pn_date.'''
        for res in multi_release:
            if (res.split('/')[-1]).startswith('data'):
                return True
        else:
            return False

    def multiRelease(self):
        '''release several jobs together.'''
        if len(self.job_name) > 1:
            os.system("rm -rf {0}".format(self.multi_checklist))
            multi_release = []
            for job_dir in self.job_name:
                self.job_dir = self.proj_dir + '/Result/' + job_dir
                self.raw_release = self.job_dir + "/data_" + \
                    job_dir + "_" + self.pn + "_" + self.now + ".tar"
                self.release = self.job_dir + "/analysis_" + \
                    job_dir + "_" + self.pn + "_" + self.now + ".tar"
                self.raw_release_dir = self.job_dir + "/data_" + \
                    job_dir + "_" + self.pn + "_" + self.now
                self.release_dir = self.job_dir + "/analysis_" + \
                    job_dir + "_" + self.pn + "_" + self.now
                release_res, checklist = self.dataRelease()
                multi_release.append(release_res)
                if self.force == 'tar':
                    # merge each job's tar checklist
                    os.system('''awk 'BEGIN{{print "check for {0}"}}{{print $0}}END{{ print "  "}}' {1} >> {2} && rm {1}'''.format(
                        job_dir, checklist, self.multi_checklist))
            # merge all release data together
            if self.force == 'tar':
                if self.rawRelease(multi_release):
                    tar_name = self.outdir + "/data_" + self.pn + "_" + self.now + ".tar"

                    os.system("tar cphf {0} {1}".format(
                        tar_name, " ".join(multi_release)))
                    #os.system('''less {0} | awk '{{print $3"\t"$6}}' > {1} && unix2dos {1}'''.format(tar_name, self.multi_checklist))
                else:
                    tar_name = self.outdir + "/analysis_" + self.pn + "_" + self.now + ".tar"

                    os.system("tar cphf {0} {1}".format(
                        tar_name, " ".join(multi_release)))
                    #os.system('''less {0} | awk '{{print $3"\t"$6}}' > {1} && unix2dos {1}'''.format(tar_name, self.multi_checklist))
                for i in multi_release:
                    os.remove(i)
                return tar_name, self.multi_checklist
            if self.force == 'not_tar':
                if self.rawRelease(multi_release):
                    release_d = self.outdir + "/data_" + self.pn + "_" + self.now
                    os.system("rm -rf {0}".format(release_d))
                    os.mkdir(release_d)
                    for dir_ in multi_release:
                        os.system("cp -a {0} {1}".format(dir_, release_d))
                    os.system(
                        "perl /NJPROJ1/DISEASE/share/Software/bin/dirCheckSize.pl {0}".format(release_d))
                    self.multi_checklist = release_d + '/checkSize.xls'
                else:
                    release_d = self.outdir + "/analysis_" + self.pn + "_" + self.now
                    os.system('rm -rf {0}'.format(release_d))
                    os.mkdir(release_d)
                    for dir_ in multi_release:
                        os.system("cp -a {0} {1}".format(dir_, release_d))
                    os.system(
                        "perl /NJPROJ1/DISEASE/share/Software/bin/dirCheckSize.pl {0}".format(release_d))
                    self.multi_checklist = release_d + '/checkSize.xls'
                for i in multi_release:
                    os.system("rm -rf {0}".format(i))
                return release_d, self.multi_checklist
        else:
            return self.dataRelease()

    def sendEmail(self):
        '''send email to project manager after checkout and tar'''
        # self.argCheck()
        release_dir, checklist = self.multiRelease()
        release_dir_list = release_dir.split("/")
        # mod to jundge release position
        release_position = "南京集群"

        # release_dir =  self.release_dir
        # check_result = self.checkRelease()
        # if check_result(from CheckResult.py) does not contain warning, show passed
        # if check_result.find("Warning") == -1:
        check_result = "Automatic check unavailable for 4.6 pepline now."
        # mail_host = "163.177.72.143"
        # mail_user = "humaninfo@novogene.com"
        # mail_pass = "DhumanB0206"
        # sender = "humaninfo@novogene.com"
        receivers = self.mail
        title = "【项目释放】" + self.pn + " " + self.proj_name
        title = title.encode('gb2312')
        # write the checkSize.xls in email with text.
        ###modyfied by yincan 20171121
        # check_list_str = ""
        # if len(self.job_name) == 1:
        #     for ck in [self.tar_checklist, self.not_tar_checklist, self.not_tar_raw_checklist]:
        #         if os.path.exists(ck):
        #             check_list_str = commands.getoutput(
        #                 '''awk '{{if($1~/^[0-9]/){{$1=$1/(1024*1024);printf "%.2f %s \\n",$1,$2;}}else{{print $0;}}}}' {0} | awk '{{if($1~/^[0-9]/){{printf "%-20s %s\\n",$1" M",$2}}else{{print $0}}}}' '''.format(ck))
        # else:
        #     check_list_str = commands.getoutput(
        #         '''awk '{{if($1~/^[0-9]/){{$1=$1/(1024*1024);printf "%.2f %s \\n",$1,$2;}}else{{print $0;}}}}' {0} | awk '{{if($1~/^[0-9]/){{printf "%-20s %s\\n",$1" M",$2}}else{{print $0}}}}' '''.format(checklist))
        # get datasize from checkSize.xls
        byte_size = commands.getoutput(
            '''cut -f 1 %s | awk '{sum += $1};END {print sum}' ''' % checklist)
        release_size = int(byte_size) / 1024 ** 3 * 1.15
        data_release = self.releaseType()
        data_release_dict = self.data_release_tran(data_release)
        tree_str = ""
        for i in ["RawData", "QC", "Mapping",
                         "Variation", "Advance", "Readme"]:
            if i in data_release_dict and i != "Advance":
                tree_str += (i + "\n")
            elif i in data_release_dict and i == "Advance":
                adv_detail = commands.getoutput(
                '''tree -dN -L 1 --noreport {0} | sed -E "s/ ->.+//g" | awk 'NR!=1 {{print}}' '''.format(self.outdir + "/" + data_release_dict["Advance"])
                )
                tree_str += i + "\n" + adv_detail + "\n"
        sample_count_dict = {}
        for i in data_release_dict:
            if i in ["RawData", "QC", "Mapping", "Variation"]:
                if i == "Variation":
                    sample_count_dict[i] = commands.getoutput(
                    ''' ls -l {0} | grep "^d" | wc -l '''.format(self.outdir +
                        "/" + data_release_dict[i]+"/SampleVariation"))
                else:
                    sample_count_dict[i] = commands.getoutput(
                    ''' ls -l {0} | grep "^d" | wc -l '''.format(self.outdir +
                        "/" + data_release_dict[i]))
                    
        sample_count_str = "\n".join([":".join([i, sample_count_dict[i]]) for i in sample_count_dict])
        email_path = self.proj_dir + '/email.txt'
        f = open(email_path, "w")
        msg_txt = u'''
Hi {0}：

一、信息自查（每个冒号后面写Yes 或 No）：
1）交付目录命名是否规范：
2）交付目录下是否有checkSize.xls文件：
3）交付目录结构是否正确：
4）是否不存在异常的空文件：
5）文件名不包含空格、#或中文等特殊字符：
--------------------总判由doublecheck填写--------------------
总判定：
备注：
（说明：1-5都是Yes，则总判定为Yes，否则为No，如果是No，需备注说明原因，返回信息分析负责人重新整理交付目录；如果总判定是Yes，则发送运营负责人。）

二、项目释放信息
1）数据所在地：{4}（选择天津集群、南京电信、南京扬子）
2）数据总大小：{5}（单位：Gbyte，是数据占用存储大小，非碱基数目）
3）样品数目(以文件夹判定, 如果有异常请进行检查):
{6}
4）数据路径：{1}
5）分期编号：{7}
6）客户名称：(如果没有请手动删除该条)

三、数据释放自动化检测结果
说明：如果数据自动化检查有问题, 会在此处显示警告信息
{2}

四、数据释放目录结构
{3}
'''.format((self.mail).split("@")[0], release_dir, check_result, tree_str,
    release_position, release_size, sample_count_str, args["pre"].strip())
        # change code to solve foxmail's problem
        msg_txt = msg_txt.encode('gb2312')
        f.write(msg_txt)
        f.close()
        # os.remove('email.txt')
        # print (msg_txt)
        # 原始邮件发送部分
        # main_text = MIMEText(msg_txt, 'plain', 'utf-8')
        # message = MIMEMultipart()
        # message.attach(main_text)
        # add attachment
        # if len(self.job_name) == 1:
        #     if os.path.exists(self.tar_checklist):
        #         file = MIMEApplication(open(self.tar_checklist, "rb").read())
        #     if os.path.exists(self.not_tar_checklist):
        #         file = MIMEApplication(
        #             open(self.not_tar_checklist, "rb").read())
        #     if os.path.exists(self.not_tar_raw_checklist):
        #         file = MIMEApplication(
        #             open(self.not_tar_raw_checklist, "rb").read())
        # else:
        #     file = MIMEApplication(open(self.multi_checklist, "rb").read())
        # file.add_header('Content-Disposition', 'attachment',
        #                 filename="checkSize.xls")
        # message.attach(file)

        # send email 原始邮件发送部分
        # message['Subject'] = Header(title, 'utf-8')
        # message['From'] = Header("disease_info@novogene.com", 'utf-8')
        # message['To'] = Header("disease-pm@novogene.com", 'utf-8')
        # smtpObj = smtplib.SMTP(mail_host)
        # smtpObj.login(mail_user, mail_pass)
        # smtpObj.sendmail(sender, receivers, message.as_string())

        # 从流程脚本中复制过来的邮件部分
        cmd = '\n. ~/.bash_profile &&'
        cmd += '\n/NJPROJ1/DISEASE/share/Software/bin/sendEmail -f humaninfo@novogene.com  -t %s -u %s -o message-file=%s message-content-type=text -o message-charset=GB2312 -xu humaninfo@novogene.com -xp DhumanB0206 -s 183.57.48.39 -o tls=no' % (receivers, title, email_path)
        rm = '&& rm -f {0}'.format(email_path) 
        os.system('ssh njlogin04 \'' + cmd + rm + '\'')
        # os.system('rm -f %s/mail.txt' % self.job_dir)        
        # print('Email has been sent successfully!')
        return checklist


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="For Data Release", formatter_class=RawTextHelpFormatter)
    parser.add_argument(
        '--pn', help="the file with info,project name.", required=True)
    parser.add_argument('--projdir', help="the path of project", required=True)
    parser.add_argument('--separate', help="tar the result separately by Rawdata QC and so on", default="N")##add 20180511
    # parser.add_argument('--seqsty', help="sequencing strategy",
    #                     choices=['WES_ag', 'WES_illu', 'WGS', 'TS'], default='WES_ag')
    parser.add_argument(
        '--mail', help="the email reciver", default="")
    parser.add_argument(
        '--yymail', help="the manager's mail", default="disease-pm@novogene.com")
    parser.add_argument('--analy_array', help="Choose the steps of ananlysis,the steps delimited by comma\n(default=1,2.1,3.1)\n"
                        "must given in analysis order\n"
                        "1\tquality_control\n"
                        "2\tmapping\n"
                        "\t2.1\tmapping_noGATK\n"
                        "\t2.2\tmapping_withGATK\n"
                        "3\tsnpindel_call\n"
                        "\t3.1\tsnpindel_call_samtools\n"
                        "\t3.2\tsnpindel_call_GATK\n"
                        "\t3.3\tsnpindel_call_samtoolsm\n"
                        "4\tsv_call\n"
                        "\t4.1\tsv_call_crest\n"
                        "\t4.2\tsv_call_breakdancer\n"
                        "5\tcnv_call\n"
                        "\t5.1\tcnv_call_freec\n"
                        "\t5.2\tcnv_call_cnvnator\n"
                        "\t5.3\tWES_cnv_CoNIFER\n"
                        "6\tAdvance\n"
                        "\t6.1\tfilterDB\n"
                        "\t6.2\tDamLevel\n"
                        "\t6.3\tFilterSV\n"
                        "\t6.4\tFilterCNV\n"
                        "\t6.5\tNoncoding\n"
                        "\t6.6\tDrug\n"
                        "7\tModelFilter\n"
                        "\t7.1\tModelFD\n"
                        "\t7.2\tModelFR\n"
                        "\t7.3\tShareCompare\n"
                        "8\tDenovo\n"
                        "\t8.1\tDenovos_samtools\n"
                        "\t8.2\tDenovos_denovogear\n"
                        "\t8.3\tDenovoF\n"
                        "\t8.4\tDenovoR\n"
                        "\t8.5\tDenovoSV\n"
                        "\t8.6\tDenovoCNV\n"
                        "9\tLinkage\n"
                        "\t9.1\tMerlinkage\n"
                        "10\tOther\n"
                        "\t10.1\tROH\n"
                        "\t10.2\tPhenolyzer\n"
                        "\t10.3\tPathway\n"
                        "\t10.4\tPPI\n"
                        "\t10.5\tHPA", default='1,2.1,3.1')
    parser.add_argument('--pre', help="The prefix of qc_list", required=True)
    parser.add_argument('--date', help="the job name",
                        default=time.strftime('%Y.%m.%d', time.localtime(time.time())) + '.job')
    parser.add_argument('--odir', help="data release directory")
    parser.add_argument('--release', type=str,
                        help='data to be released', default="")
    parser.add_argument('--force', type=str, choices=['tar', 'not_tar'], default="",
                        help="tar/not_tar means tar/not_tar release regardless of the datasize")
    parser.add_argument('--no-record', action='store_true',
                        help="do not record to database")
    #parser.add_argument('--multijob',type=str,help="release several jobs together",default="",required=False)
    args = vars(parser.parse_args())
    Release = Release()
    releasetype=Release.releaseType()
    #self.raw_tardir=self.outdir + "/data_" + self.pn + "_" + args["pre"].strip()+"_"+self.now
    #self.raw_release = self.raw_tardir + "/data_" + self.pn + "_" + args["pre"].strip()+"_"+self.now + ".tar"   
    #self.tardir = self.outdir + "/analysis_" + self.pn + "_" +args["pre"].strip()+"_"+ self.now
    #self.release = self.tardir + "/analysis_" + self.pn + "_" +args["pre"].strip()+"_"+ self.now + ".tar"
    #self.raw_release_dir = self.outdir + "/data_" + self.pn + "_" + self.now
    #self.release_dir = self.outdir + "/" + "analysis_" + self.pn + "_" + self.now
    # if --mail not given, set to path owner's email
    if Release.mail == "":
        ls_output = commands.getoutput('ls -ld {0}'.format(Release.job_dir))
        path_owner = ls_output.strip().split(' ')[2]
        Release.mail = path_owner + '@novogene.com'
        print(Release)
    checklist = Release.sendEmail()
    if args['no_record']:
        exit('done')

    if args["release"]:
        print('record to database')
        cmd = "\n. ~/.bash_profile &&"
        cmd += "\n python /NJPROJ2/DISEASE/share/Disease/Result/Version_4.6/Record_Data_releaseV2.0.py  --projdir {0} --analy_array {1} --odir {2} --date {3} --pre {4} --pn {5} --yymail {6} --release {7}".format(args["projdir"],args["analy_array"],args["odir"],args["date"],args["pre"],args["pn"],args["yymail"],args["release"])
        os.system('ssh njlogin04 \'' + cmd +'\'')
        #os.system("python /NJPROJ2/DISEASE/share/Disease/Result/Version_4.6/Record_Data_releaseV2.0.py  --projdir " + args["projdir"] + " --analy_array " + args["analy_array"] + " --odir " + args["odir"] +  " --date " + args["date"] + " --pre " + args["pre"] + "  --pn " + args["pn"] + " --yymail " + args["yymail"] + " --release " + args["release"])
    else:
        cmd = "\n. ~/.bash_profile &&"
        cmd += "\n python /NJPROJ2/DISEASE/share/Disease/Result/Version_4.6/Record_Data_releaseV2.0.py  --projdir {0} --analy_array {1} --odir {2} --date {3} --pre {4} --pn {5} --yymail {6} --release {7}".format(args["projdir"],args["analy_array"],args["odir"],args["date"],args["pre"],args["pn"],args["yymail"],",".join(Release.releaseType()))
        os.system('ssh njlogin04 \'' + cmd +'\'')
        #os.system("python /NJPROJ2/DISEASE/share/Disease/Result/Version_4.6/Record_Data_releaseV2.0.py  --projdir " + args["projdir"] + " --analy_array " + args["analy_array"] + " --odir " + args["odir"] +  " --date " + args["date"] + " --pre " + args["pre"] + "  --pn " + args["pn"] + " --yymail " + args["yymail"] + " --release " + ",".join(Release.releaseType()))
