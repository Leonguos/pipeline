#!/usr/bin/env python
# -*- coding=utf-8 -*-
from mongoengine import connect, DynamicDocument


connect(
    db='project',
    host='192.168.20.19',
    port=8888,
    username='disease',
    password='novogene')
# connect(db='project', host='192.168.20.19', port=8085, username='normal', password='novogene2017')


class pathlist(DynamicDocument):

    pass


class samplelist(DynamicDocument):

    pass

class projprepare(DynamicDocument):

    pass


if __name__ == "__main__":

    # test
    # print pathlist.objects()[0].to_mongo().to_dict()
    print pathlist.objects(
        projectid='P101SC18070637-01',
        novoid='ND180701860')[0].to_mongo().to_dict()
