import subprocess
from IPython.display import display, Markdown, HTML, Image

__all__ = ['display_header']

import logging
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
 } else {
 $('div.input').show();
 }
 code_show = !code_show
}
$( document ).ready(code_toggle);
</script>
<form action="javascript:code_toggle()"><input type="submit" value="Click here to toggle on/off the display of raw code."></form>''')

logo = Image('logos/logo.png', height=80, width=80)

def get_marxs_status():
    try:
        import marxs.version
    except ImportError:
        return 'MARXS cannot be imported. No version information is available.'
    return 'MARXS ray-trace code version {} (commit hash: {} from {})'.format(marxs.version.version,
                                                                             marxs.version.githash[:10],
                                                                             marxs.version.timestamp.date())

def get_arcus_status():
    try:
        with DisableLogger():
            import arcus.version
    except ImportError:
        return 'ARCUS cannot be imported. No version information is available.'
    return 'ARCUS python code version {} (commit hash: {} from {})'.format(arcus.version.version,
                                                                          arcus.version.githash[:10],
                                                                          arcus.version.timestamp.date())

def get_caldb_status():
    try:
        with DisableLogger():
            from arcus import load_csv
    except ImportError:
        return 'ARCUS CALDB cannot be imported. No version information is available.'
    return 'ARCUS CALDB version ' + load_csv.string_git_info()


def get_nb_status(filename):
    try:
        gitlog = subprocess.check_output(['git',  'log', '-1', '--use-mailmap',
                                          '--format=medium',  '--', filename])
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
        modified = filename in subprocess.check_output(['git', 'ls-files', '-m'])
        if modified:
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
- {arcus}
- {caldb}
'''.format(nbstatus=get_nb_status(filename),
           marxs=get_marxs_status(),
           arcus=get_arcus_status(),
           caldb=get_caldb_status())
    return Markdown(out)

def display_header(filename, status=None):
    display(logo)
    display(revision_status(filename, status=status))
    display(codetoggle)
