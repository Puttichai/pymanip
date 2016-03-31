from setuptools import setup
import os

setup(name='pymanip',
      packages=['pymanip', 
                'pymanip.planners', 
                'pymanip.utils'],
      include_package_data=True)

# Create a pymanip directory to store various files
home = os.path.expanduser('~')
pymanipdir = os.path.join(home, '.pymanip')
if not os.path.exists(pymanipdir):
    os.makedirs(pymanipdir)

# Change permission to pymanip directory
uid = int(os.environ.get('SUDO_UID'))
gid = int(os.environ.get('SUDO_GID'))
os.chown(pymanipdir, uid, gid)
