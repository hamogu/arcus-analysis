# -*- coding: utf-8 -*-
'''Compile page showing simulations for ARCUS from various input directories.
'''
from __future__ import print_function

from glob import glob
import os
from os.path import join as pjoin
import shutil
import json
import yaml
from copy import copy

from jinja2 import Environment, FileSystemLoader, select_autoescape

try:
    import configparser  # Py 3
except ImportError:
    import ConfigParser as configparser  # Py 2

cfgpath = [os.path.join(os.path.dirname(__file__), '..', 'site.cfg')]
'Path list to search for configuration files.'
conf = configparser.ConfigParser()
cfgfile = conf.read(cfgpath)
basepath = conf.get("Website", "htmlbase")

# Generate html
env = Environment(loader=FileSystemLoader([pjoin(basepath, 'htmltemplates')]),
                  autoescape=select_autoescape(['html']))

outpath = conf.get("Website", "outpath")
if not os.path.exists(outpath):
    os.makedir(outpath)

with open(conf.get("Website", 'toc'), 'r') as f:
    pagelist = yaml.load(f)


# The page list is nested to represent the structure of the menu.
# Make a flattened copy that can be iterated over:
def add_children(flatlist, page):
    flatlist.append(page)
    if 'children' in page:
        for p in page['children']:
            add_children(flatlist, p)

flatlist = []
add_children(flatlist, pagelist)

for page in flatlist:
    if 'title' not in page:
        continue
    print("Working on {0}".format(page['title']))
    kwargs = copy(page)
    if page['type'] == 'included':
        template = env.get_template(page['href'] + '.html')
    else:
        template = env.get_template(page['type'] + '.html')
        if page['type'] == 'x3d':
            x3dpath = conf.get("Website", "x3dpath")
            with open(pjoin(x3dpath, page['href'] + '.json'), 'r') as f:
                kwargs['x3d'] = json.load(f)
            kwargs['x3d']['path'] = pjoin('x3d', page['href'])
            outdir = pjoin(outpath, 'x3d')
            if not os.path.exists(outdir):
                os.makedirs(outdir)
            shutil.copy(pjoin(x3dpath, page['href'] + '.x3d'),
                        pjoin(outdir, page['href'] + '.x3d'))
        elif page['type'] == 'notebook':
            nbpath = conf.get("Website", "notebookpath")
            outdir = pjoin(outpath, 'notebooks')
            if not os.path.exists(outdir):
                os.makedirs(outdir)
            shutil.copy(pjoin(nbpath, page['href'] + '.html'),
                        pjoin(outdir, page['href'] + '.html'))
            kwargs['path'] = pjoin('notebooks', page['href'])

    with open(pjoin(outpath, page['href'] + '.html'), "w") as f:
        f.write(template.render(pages=pagelist['children'],
                                active_page=page['href'],
                                **kwargs))


# copy several directories verbatim
for d in ['css', 'images']:
    outdir = pjoin(outpath, d)
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    filelist = glob(pjoin(basepath, d, '*'))
    for f in filelist:
        shutil.copy(f, outdir)

print("Done. Website is in directory: {}.".format(outpath))
print("If this version is meant for publication, take the following steps:")
print(" cd {}".format(outpath))
print(" git add any_new_files")
print(" git commit -am'My message here'")
print(" scp -r * space:/space/web/home/guenther/ARCUS")
