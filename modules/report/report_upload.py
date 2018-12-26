#!/usr/bin/env python2
# -*- coding=utf-8 -*-
import glob
import argparse
import paramiko
from scp import SCPClient


def main():

    context = {
        'hostname': args['hostname'],
        'port': args['port'],
        'username': args['username'],
        'password': args['password'],
    }
    ssh = paramiko.SSHClient()
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    ssh.connect(**context)

    print 'connected to remote: {hostname}'.format(**args)

    ssh.exec_command('mkdir -p {remote_path}'.format(**args))

    print 'uploading reports to {hostname}:{remote_path}'.format(**args)

    with SCPClient(ssh.get_transport(), socket_timeout=60) as scp:
        for each in glob.glob('{files}/*'.format(**args)):
            print '> uploading {} ...'.format(each)
            scp.put(each, args['remote_path'], recursive=True)

    ssh.close()


if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument('files', help='a single path, or a list of paths to be transferred.')
    parser.add_argument('remote_path', help='the remote path to receive files', default='aaa')

    parser.add_argument('-host', '--hostname', help='the hostname of remote server', default='192.168.20.19')
    parser.add_argument('-port', help='the port to login', default=22, type=int)
    parser.add_argument('-u', '--username', help='the username to login', default='data')
    parser.add_argument('-p', '--password', help='the password to login', default='novogene2017')
    
    args = vars(parser.parse_args())

    # print args

    main()
