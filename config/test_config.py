#!/usr/bin/env python
# -*- coding=utf-8 -*-
import os
from ConfigParser import ConfigParser


BASE_DIR = os.path.dirname(os.path.abspath(__file__))

CONFIG = ConfigParser()

CONFIG.read(os.path.join(BASE_DIR, 'config.ini'))


print CONFIG.get('software', 'genemaniaDIR')
print CONFIG.get('genome_mm9', 'gatkRecalSnpVCF')
