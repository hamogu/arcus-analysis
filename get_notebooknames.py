from __future__ import print_function
import yaml
import argparse


def get_nbnames(namelist, element):
    if 'children' not in element:
        if element['type'] == 'notebook':
            namelist.append(element['href'])
    else:
        for c in element['children']:
            get_nbnames(namelist, c)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='''
    Run plotting scripts to make output for website.

    For example:
    > python get_notebooknames.py "../html/" ".html"
    will list the names of all html files for the website that are made by
    rendering notebooks.
    ''')
    parser.add_argument('prename',
                        help='directory patter to be prefixed')
    parser.add_argument('postname',
                        help='file extension to be post-fixed')
    args = parser.parse_args()

    with open('website/toc.yaml', 'r') as f:
        pagelist = yaml.load(f)

    names = []
    get_nbnames(names, pagelist)
    names = [args.prename + n + args.postname for n in names]
    print(' '.join(names))
