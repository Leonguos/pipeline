#!/usr/bin/env python2
# -*- coding=utf-8 -*-
import glob
import argparse
import paramiko
from scp import SCPClient, SCPException


def main():

    context = {
        'hostname': args['hostname'],
        'port': args['port'],
        'username': args['username'],
        'password': args['password'],
    }

    with paramiko.SSHClient() as ssh:
        ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        ssh.connect(**context)

        print 'connected to remote: {hostname}'.format(**args)

        stdin, stdout, stderr = ssh.exec_command('mkdir -p {remote_path}'.format(**args))

        print stdout.read()

        print 'uploading reports to {hostname}:{remote_path}'.format(**args)

        with SCPClient(ssh.get_transport(), socket_timeout=60) as scp:
            for each in glob.glob('{files}/*'.format(**args)):
                print '> uploading {} ...'.format(each)
                n = 0
                while n < 2:
                    try:
                        scp.put(each, args['remote_path'], recursive=True)
                        n += 1
                    except SCPException as e:
                        print e
                        n += 1


if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument('files', help='a single path, or a list of paths to be transferred.')
    parser.add_argument('remote_path', help='the remote path to receive files')

    # parser.add_argument('-host', '--hostname', help='the hostname of remote server', default='39.106.135.106')
    parser.add_argument('-host', '--hostname', help='the hostname of remote server', default='172.17.8.62')
    parser.add_argument('-port', help='the port to login', default=22, type=int)
    parser.add_argument('-u', '--username', help='the username to login', default='report')
    parser.add_argument('-p', '--password', help='the password to login', default='disease-report-2019')
    
    args = vars(parser.parse_args())

    # print args

    main()
