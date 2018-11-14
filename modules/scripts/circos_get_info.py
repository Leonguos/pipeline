#!/usr/bin/env python
# -*- coding=utf-8 -*-
import os
import re
import sys
import math
import argparse

from collections import defaultdict

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(os.path.dirname(BASE_DIR)))

import utils
from config import config


TYPE = {
    'sv': ('crest', 'breakdancer', 'lumpy'),
    'cnv': ('freec', 'cnvnator'),
    'mutation': ('snp', 'indel'),
    'depth': ('bam',)
}


CHROM_LENGTH = {
    '1': 249250621,
    '2': 243199373,
    '3': 198022430,
    '4': 191154276,
    '5': 180915260,
    '6': 171115067,
    '7': 159138663,
    '8': 146364022,
    '9': 141213431,
    '10': 135534747,
    '11': 135006516,
    '12': 133851895,
    '13': 115169878,
    '14': 107349540,
    '15': 102531392,
    '16': 90354753,
    '17': 81195210,
    '18': 78077248,
    '19': 59128983,
    '20': 63025520,
    '21': 48129895,
    '22': 51304566,
    'X': 155270560,
    'Y': 59373566
}


class GetInfo(object):

    def __init__(self, *args, **kwargs):
        # self.__dict__.update(kwargs)
        self.infile = kwargs['infile']
        self.vtype = kwargs['vtype']
        self.outdir = kwargs['outdir']
        self.region_length = kwargs['region_length']
        self.window_size = kwargs['window_size']

        self.outsuffix = os.path.splitext(os.path.basename(self.infile))[0]

        self.igvtools = config.CONFIG.get('software', 'igvtools')

        self.normal_chrom = [str(n) for n in range(1, 23)] + ['X', 'Y']

        # print self.__dict__

    def start(self):

        print 'extracting informations for {}'.format(self.vtype)

        if self.vtype in TYPE['sv']:
            self.sv_info()

        elif self.vtype in TYPE['cnv']:
            self.cnv_info()

        elif self.vtype in TYPE['mutation']:
            self.mutation_info()

        elif self.vtype in TYPE['depth']:
            self.depth_info()

        else:
            exit('error vtype: {}'.format(self.vtype))

    def sv_info(self):

        sv_context = defaultdict(list)
        svid_list = []

        with utils.safe_open(self.infile) as f:
            print 'open file: {}'.format(self.infile)
            for line in f:
                linelist = line.strip().split('\t')
                if linelist[0] in ('Chr',):
                    headerlist = linelist
                    continue

                chrom = linelist[headerlist.index('Chr')].strip('Chr').strip('chr')
                start = linelist[headerlist.index('Start')]
                end = linelist[headerlist.index('End')]
                func = linelist[headerlist.index('Func')]
                tchr = linelist[headerlist.index('TCHR')].strip('Chr').strip('chr')
                tstart = linelist[headerlist.index('TSTART')]
                svid = linelist[headerlist.index('SVID')]
                svtype = linelist[headerlist.index('SVType')]

                # 对于同一个ID，只记录第一条
                if svid in svid_list:
                    continue

                # 只保留常染色体+XY
                if any(each not in self.normal_chrom + ['na'] for each in [chrom, tchr]):
                    # sys.stderr.write('skip a line of unnormal chrom: ' + line)
                    continue

                # 只保留exonic或splicing区的变异
                if not (func.startswith('exonic') or func.startswith('splicing')):
                    # sys.stderr.write('skip a line of not exonic or splicing: ' + line)
                    continue

                # 跳过breakpoint行
                if svtype == 'breakpoint':
                    continue

                # for lumpy
                if svtype == 'DUP':
                    svtype = 'INS'

                # for breakdancer
                if svtype not in ('CTX', 'ITX', 'DEL', 'INS', 'INV'):
                    svtype = linelist[headerlist.index('TX')][:3].upper()

                end1 = int(start) + 1
                if svtype in ('CTX', 'ITX'):
                    chrom2 = tchr
                    start2 = tstart
                    end2 = int(start2) + 1
                elif svtype in ('DEL', 'INS', 'INV'):
                    chrom2 = chrom
                    start2 = end
                    end2 = int(start2) + 1

                info = 'hs{chrom} {start} {end1} hs{chrom2} {start2} {end2}'.format(**locals())

                svid_list.append(svid)
                sv_context[svtype].append(info)

        for svtype in ('CTX', 'ITX', 'DEL', 'INS', 'INV'):
            outfile = '{outdir}/{vtype}_{svtype}_{outsuffix}'.format(**dict(self.__dict__, **locals()))
            with utils.safe_open(outfile, 'w') as out:
                for info in sv_context[svtype]:
                    out.write(info + '\n')
            print 'write file: {}'.format(outfile)

    def cnv_info(self):

        outfile = '{outdir}/{vtype}_{outsuffix}'.format(**self.__dict__)

        with utils.safe_open(self.infile) as f, utils.safe_open(outfile, 'w') as out:
            print 'open file: {}'.format(self.infile)
            for line in f:
                linelist = line.strip().split('\t')
                if linelist[0] in ('Chr', '#Chr'):
                    headerlist = linelist
                    continue
                chrom = linelist[0].strip('Chr').strip('chr')
                start = linelist[headerlist.index('Start')]
                end = linelist[headerlist.index('End')]

                if chrom not in self.normal_chrom:
                    continue

                if self.vtype == 'freec':
                    copynumber = int(linelist[headerlist.index('CopyNumber')]) - 2
                    copynumber = 6 if copynumber > 6 else copynumber
                elif self.vtype == 'cnvnator':
                    copynumber = float(linelist[headerlist.index('RD')])
                    copynumber = 3 if copynumber > 3 else copynumber

                line = 'hs{chrom} {start} {end} {copynumber}\n'.format(**locals())
                out.write(line)

        print 'write file: {}'.format(outfile)

    def depth_info(self):

        cmd = '{igvtools} count -w {window_size} {infile} {outdir}/depthRaw_{outsuffix}.wig novo37'.format(**self.__dict__)
        print 'run cmd:', cmd
        assert not os.system(cmd)

        out_wig = '{outdir}/depthRaw_{outsuffix}.wig'.format(**self.__dict__)
        outfile = '{outdir}/depth_{outsuffix}'.format(**self.__dict__)

        with utils.safe_open(out_wig) as f, utils.safe_open(outfile, 'w') as out:
            for line in f:
                if line.startswith('track'):
                    continue
                elif line.startswith('variableStep'):
                    chrom = re.findall(r'chrom=(.+?) ', line)[0].strip('chr')
                    continue
                linelist = line.strip().split('\t')
                start = int(linelist[0])
                depth = float(linelist[1])

                depth_log10 = math.log10(depth + 1)

                end = start + self.window_size - 1
                if start + self.window_size > CHROM_LENGTH[chrom]:
                    end = CHROM_LENGTH[chrom]

                line = 'hs{chrom}\t{start}\t{end}\t{depth_log10}\n'.format(**locals())
                out.write(line)
                
        print 'write file: {}'.format(outfile)

    def mutation_info(self):

        chrom_region = self._get_chrom_region()

        # print chrom_region['1'].items()[0]

        with utils.safe_open(self.infile) as f:
            print 'open file: {}'.format(self.infile)
            for line in f:
                linelist = line.strip().split('\t')
                if linelist[0] == '#CHROM':
                    headerlist = linelist
                if line.startswith('#'):
                    continue
                chrom = linelist[headerlist.index('#CHROM')]
                pos = int(linelist[headerlist.index('POS')])

                now_region = self._get_now_region(chrom, pos, chrom_region)
                # print chrom, now_region

                if chrom not in self.normal_chrom:
                    continue

                if self.vtype == 'snp':
                    genotype = self._get_genotype(headerlist, linelist)
                    if genotype == 'hom':
                        chrom_region[chrom][now_region][0] += 1
                    else:
                        chrom_region[chrom][now_region][1] += 1
                elif self.vtype == 'indel':
                    chrom_region[chrom][now_region][0] += 1

        density_outfile = '{outdir}/{vtype}_density_{outsuffix}'.format(**self.__dict__)
        if self.vtype == 'snp':
            snp_hom_het_ratio_outfile = '{outdir}/{vtype}_ratio_{outsuffix}'.format(**self.__dict__)
            snp_hom_het_ratio_out = utils.safe_open(snp_hom_het_ratio_outfile, 'w')

        chrom_order = map(str, range(1, 23)) + ['X', 'Y']
        with utils.safe_open(density_outfile, 'w') as density_out:
            for chrom, regions in sorted(chrom_region.iteritems(), key=lambda (k, v): chrom_order.index(k)):
                for start, end in sorted(regions):
                    site_number = sum(regions[start, end])
                    density = float(site_number) /  self.region_length
                    line = 'hs{chrom}\t{start}\t{end}\t{density}\n'.format(**locals())
                    density_out.write(line)

                    if self.vtype == 'snp':
                        hom_ratio = het_ratio = 0
                        if site_number:
                            hom_ratio = regions[start, end][0] / float(site_number)
                            het_ratio = regions[start, end][1] / float(site_number)
                        line = 'hs{chrom}\t{start}\t{end}\t{hom_ratio},{het_ratio}\n'.format(**locals())
                        snp_hom_het_ratio_out.write(line)

        print 'write file: {}'.format(density_outfile)

        if self.vtype == 'snp':
            snp_hom_het_ratio_out.close()
            print 'write file: {}'.format(snp_hom_het_ratio_outfile)
            

    def _get_genotype(self, headerlist, linelist):

        format_idx = headerlist.index('FORMAT')
        format = linelist[format_idx]
        sample_info = linelist[format_idx + 1]
        genotype_idx = format.split(':').index('GT')
        genotype = sample_info.split(':')[genotype_idx]
        # print genotype

        if genotype.split('/')[0] == genotype.split('/')[1]:
            return 'hom'

        return 'het'

    def _get_chrom_region(self):

        chrom_region = defaultdict(dict)

        for chrom, length in CHROM_LENGTH.iteritems():
            for start in range(1, length+1, self.region_length):
                end = start + self.region_length - 1
                end = length if end > length else end
                chrom_region[chrom][start, end] = [0, 0]  # means [hom_num, het_num]

        return chrom_region

    def _get_now_region(self, chrom, pos, chrom_region):

        for start, end in chrom_region[chrom]:
            if start <= pos <= end:
                return start, end
        print 'can not find a region from {chrom}:{pos}'.format(**locals())
        exit(1)


def main():

    GetInfo(**args).start()


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description='extract info for circos drawing')

    parser.add_argument('-i', '--infile', help='the input file', required=True)
    parser.add_argument('-o', '--outdir', help='the output directory[default=%(default)s]', default='.')
    parser.add_argument(
        '-vt',
        '--vtype',
        help='the input file type',
        choices=[x for sublist in TYPE.values() for x in sublist])
    
    parser.add_argument(
        '-len',
        '--region-length',
        help='the length of region to calculate mutation density[default=%(default)s]',
        default=1000000)
    
    parser.add_argument(
        '-w',
        '--window-size',
        help='the window size for igvtools count to calclate depth[default=%(default)s]',
        default=500000)

    args = vars(parser.parse_args())

    main()
