import os
import sys
import paramiko
import datetime


hostname = '192.168.20.19'
username = 'data'
password = 'novogene2017'
port = 22

def mkdir_p(sftp, remote_directory):
    """Change to this directory, recursively making new folders if needed.
    Returns True if any folders were created."""
    if remote_directory == '/':
        # absolute path so change directory to root
        sftp.chdir('/')
        return
    if remote_directory == '':
        # top-level relative directory must exist
        return
    try:
        sftp.chdir(remote_directory) # sub-directory exists
    except IOError:
        dirname, basename = os.path.split(remote_directory.rstrip('/'))
        mkdir_p(sftp, dirname) # make parent directories
        sftp.mkdir(basename) # sub-directory missing, so created it
        sftp.chdir(basename)
        return True


def upload(local_dir, remote_dir):

    try:
        t = paramiko.Transport((hostname, port))
        t.connect(username=username, password=password)
        sftp = paramiko.SFTPClient.from_transport(t)
        print 'upload file start %s ' % datetime.datetime.now()
        for root, dirs, files in os.walk(local_dir):
            for filespath in files:
                local_file = os.path.join(root, filespath)
                a = local_file.replace(local_dir, '').lstrip('/')
                remote_file = os.path.join(remote_dir, a)
                # print local_file
                # print remote_file
                try:
                    sftp.put(local_file, remote_file)
                except Exception, e:
                    print e
                    mkdir_p(sftp, os.path.split(remote_file)[0])
                    sftp.put(local_file, remote_file)
                # print "upload %s to remote %s" % (local_file, remote_file)
            for name in dirs:
                local_path = os.path.join(root, name)
                a = local_path.replace(local_dir, '')
                remote_path = os.path.join(remote_dir, a)
                try:
                    mkdir_p(sftp, remote_path)
                    print "mkdir path %s" % remote_path
                except Exception, e:
                    print e
        print 'upload file success %s ' % datetime.datetime.now()
        t.close()
    except Exception, e:
        print e

upload(sys.argv[1], sys.argv[2])
