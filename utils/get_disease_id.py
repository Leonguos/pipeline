#!/usr/bin/env python
# -*- coding=utf-8 -*-
import os
import json
import fuzzywuzzy.process


def get_disease_id(diseases, data):

    if not diseases:
        return None

    disease_ids = []

    database = os.path.basename(data).split('.')[0]
    print 'Use database: {} [ {} ]'.format(database, data)

    with open(data) as f:
        diseaseDict = json.load(f)

    for disease in diseases.split(';'):
        result = fuzzywuzzy.process.extractBests(disease, diseaseDict.keys())
        disease, score = result[0]
        disease_id = diseaseDict.get(disease)
        disease_ids.append(disease_id)
        print 'Matched disease: %s[ID:%s] (sorce:%s)' % (disease, disease_id, score)
        if int(score) < 90:
            print 'But the score is too low, maybe you should try again with another name'

        if int(score) < 95:
            print 'Other matched results: \n\t%s\n\t%s\n\t%s\n\t%s ' % tuple(
                map(lambda x:'{disease}[ID:{disease_id}] (score:{score})'.format(
                    disease=x[0],
                    disease_id=diseaseDict.get(x[0]),
                    score=x[1]), result[1:6]))

    return ';'.join(disease_ids)


if __name__ == "__main__":

    print get_disease_id('heart disease')
    print get_disease_id('heart disease;Neuromyelitis Optica')