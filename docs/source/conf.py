import re
import sys
from os.path import join, abspath
sys.path.insert(0, abspath("../../"))


def get_version():
    VERSIONFILE = join('..', '..', 'Cargo.toml')
    with open(VERSIONFILE, 'rt') as f:
        lines = f.readlines()
    vgx = '^version = \"[0-9+.0-9+.0-9+]*[a-zA-Z0-9]*\"'
    for line in lines:
        mo = re.search(vgx, line, re.M)
        if mo:
            return mo.group().split('"')[1]
    raise RuntimeError('Unable to find version in %s.' % (VERSIONFILE,))


project = 'Polymers'
version = get_version()
release = version
author = 'Michael R. Buche'
copyright = '2022 National Technology & Engineering Solutions of Sandia, LLC'

add_module_names = False
bibtex_bibfiles = ['main.bib']
bibtex_default_style = 'plain'
extensions = [
    'matplotlib.sphinxext.plot_directive',
    'nbsphinx',
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
    'sphinxcontrib.bibtex',
    'sphinx_copybutton'
]
html_css_files = [
    'custom.css'
]
html_show_sphinx = False
html_show_sourcelink = False
html_static_path = ['_static']
html_theme = 'sphinx_rtd_theme'
html_theme_options = {
    'navigation_depth': 8
}
latex_engine = 'xelatex'
nbsphinx_allow_errors = True
plot_html_show_formats = False
plot_html_show_source_link = False
plot_include_source = True
plot_rcparams = {'font.size': 10}
plot_formats = [('png', 300)]
templates_path = ['_templates']
