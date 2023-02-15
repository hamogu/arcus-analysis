'''This module collects utilities for the Arcus notebooks.'''
import sys
import os
import subprocess
import configparser
import re
import logging

from IPython.display import display, Markdown, HTML, Image

__all__ = ['display_header', 'get_path', 'cfgpath', 'display_codetoggle']

cfgpath = [os.path.join(os.path.dirname(sys.modules[__name__].__file__), '..', 'site.cfg')]
'Path list to search for configuration files.'


class DisableLogger():
    def __enter__(self):
        logging.disable(logging.CRITICAL)

    def __exit__(self, a, b, c):
        logging.disable(logging.NOTSET)


codetoggle = HTML('''<script>
code_show=true;
function code_toggle() {
 if (code_show){
 $('div.input').hide();
 $('div.prompt').hide();
 } else {
 $('div.input').show();
 $('div.prompt').show();
}
 code_show = !code_show
}
$( document ).ready(code_toggle);
</script>
<form action="javascript:code_toggle()"><input type="submit" value="Click here to toggle on/off the display of raw code."></form>''')

logo = Image('logos/logo.png', height=80, width=80)

reexp = re.compile(r"(?P<version>[\d.dev]+[\d]+)[+]?(g(?P<gittag>\w+))?[.]?(d(?P<dirtydate>[\d]+))?")


def parse_git_scm_version(version):
    ver = reexp.match(version)
    v = ver.group('version')
    if not ver.group('gittag') is None:
        v += ' (commit hash: ' + ver.group('gittag') + ')'
    if not ver.group('dirtydate') is None:
        v += ' (Date of dirty version: ' + ver.group('dirtydate') + ')'
    return v


def get_marxs_status():
    try:
        import marxs.version as mvers
    except ImportError:
        return 'MARXS cannot be imported. No version information is available.'
    return 'MARXS ray-trace code version ' + parse_git_scm_version(mvers.version)


def get_caldb_status():
    try:
        with DisableLogger():
            from marxs.missions.arcus import utils
    except ImportError:
        return 'ARCUS CALDB cannot be imported. No version information is available.'
    return 'ARCUS CALDB version ' + utils.git_hash.decode() + ' (' + utils.git_info.decode()[:19] + ')'


def get_nb_status(filename):
    try:
        gitlog = subprocess.check_output(['git',  'log', '-1', '--use-mailmap',
                                          '--format=medium',  '--', filename])
        gitlog = gitlog.decode(sys.stdout.encoding)
    except subprocess.CalledProcessError:
        return '''git is not installed or notebook was run outside of git version control.
        No versioning information can be displayed.'''

    if len(gitlog) == 0:
        return '''file: {} not found in repository (path missing or new file not yet commited?).
        No versioning information can be displayed.'''.format(filename)
    else:
        gitlog = gitlog.split('\n')
        out = '''Last revision in version control:

- {1}
- {0}
- {2}
'''.format(gitlog[0], gitlog[1], gitlog[2])
        modifiedfiles = subprocess.check_output(['git', 'ls-files', '-m'])
        modifiedfiles = modifiedfiles.decode(sys.stdout.encoding)
        if filename in modifiedfiles:
            out = out + '''
**The version shown here is modified compared to the last commited revision.**

            '''
    return out


def revision_status(filename, status=None):
    if status is None:
        statusstring = ''
    else:
        statusstring = ': *{}*'.format(status)
    out='### Revision status{}\n'.format(statusstring)

    out = out + '''{nbstatus}

This document is git version controlled. The repository is available at https://github.com/hamogu/arcus.
See git commit log for full revision history.

Code was last run with:

- {marxs}
- {caldb}
'''.format(nbstatus=get_nb_status(filename),
           marxs=get_marxs_status(),
           caldb=get_caldb_status())
    return Markdown(out)


def display_header(filename, status=None):
    display(logo)
    display(revision_status(filename, status=status))


def get_path(name):
    '''Get path name info from site.cfg file in root directory.

    If a path does not exist, it will be created.
    '''
    conf = configparser.ConfigParser()
    cfgfile = conf.read(cfgpath)
    if not cfgfile:
        raise Exception("No config file with path specifications found. File must be called 'site.py' and be located in one of the following directories: {}".format(cfgpath))
    path = conf.get("Path", name)
    if not os.path.exists(path):
        os.makedirs(path)
    return path


def display_codetoggle():
    '''Display button "show code on/off".

    Calling this function toggles off. Call this at the end of a notebook after
    all processing steps have run.  Otherwise, newly run cell will have
    prompts.
    '''
    display(codetoggle)
