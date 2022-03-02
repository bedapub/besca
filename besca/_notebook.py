"""Helper functions for working with notebooks"""

from IPython.display import display, Javascript
from glob import glob
import os
import subprocess


def save_notebook():
    """Use JavaScript to simulate saving notebook, which makes sure
       that the current notebook is the last one in the current directory
       that is modified"""

    display(Javascript(
        "document.body.dispatchEvent("
        "new KeyboardEvent('keydown', {key:'s', keyCode: 83, ctrlKey: true}"
        "))"
    ))


def save_notebook_return_path():
    """Save the current notebook and return its full path"""

    save_notebook()
    ipynbs = glob("*.ipynb")
    curr_dir = os.getcwd()
    max_mtime = 0
    for fname in ipynbs:
        full_path = os.path.join(curr_dir, fname)
        mtime = os.stat(full_path).st_mtime
        if mtime > max_mtime:
            max_mtime = mtime
            max_file = full_path
    return max_file


def convert_notebook_to_HTML():
    """Convert the current notebook to HTML"""

    current = save_notebook_return_path()
    res = subprocess.run(['jupyter', 'nbconvert',
                          '--to', 'html', current],
                         shell=False, capture_output=True)
    return(res)
