# -*- coding: utf-8 -*-
"""
Created on Mon Jun  6 12:22:14 2016

@author: cornkle
"""

import re
import os
import glob

def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    '''
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    (See Toothy's implementation in the comments)
    '''
    return [ atoi(c) for c in re.split('(\d+)', text) ]


def dir_names(data_path, dir_only=False, file_only=False):
    """
    In datapath, looks for all dir and filenames and returns the basenames
    """
    dir_basename = []
    for d in glob.glob(data_path + '/*'):
        if dir_only:
            if not os.path.isdir:
                continue
        if file_only:
            if not os.path.isfile:
                continue
                
        dir_basename.append(os.path.basename(d))
    return dir_basename



def tolist(arg, length=None):
    """Makes sure that arg is a list."""

    if isinstance(arg, str):
        return [arg]

    try:
        (e for e in arg)
    except TypeError:
        arg = [arg]

    arg = list(arg)

    if length is not None:

        if len(arg) == 1:
            arg *= length
        elif len(arg) == length:
            pass
        else:
            raise ValueError('Cannot broadcast len {} '.format(len(arg)) +
                             'to desired length: {}.'.format(length))

    return arg
