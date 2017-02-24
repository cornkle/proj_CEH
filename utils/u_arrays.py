# -*- coding: utf-8 -*-
"""
Created on Tue Mar 15 17:27:16 2016

@author: cornkle
"""

import os, fnmatch
import numpy as np
from math import sqrt

def locate(pattern, root_path):
    strg = []
    llist = os.listdir(root_path)
    llist.sort()
    for file in llist:
        if file.endswith(pattern):
            strg.append(os.path.join(root_path, file))
    return strg


def distance(x1, y1, x2, y2):
    return np.sqrt((x1 - x2) *(x1 - x2) + (y1 - y2) * (y1 - y2))


def closest_point(point, points):
    dist_2 = np.sum((points - point) * (points - point), axis=1)
    return np.argmin(dist_2)


"""create one unique integer from two positive integers
 Cantor pairing function"""
def unique_of_pair(x,y):

    uni = (x + y) * (x + y + 1) / 2 + y
    return uni


"""
Find all indices within the local circle of radius
Input:
x: x index of center point
y: y index of center point
radius: radius in pixels, floats are handled including the farthest point
Returns a tuple of (y index, x index)
"""
def draw_circle(x, y, radius):

    xloc1 = np.arange(x - radius, x + radius + 1)
    yloc1 = np.arange(y - radius, y + radius + 1)
    xloc, yloc = np.meshgrid(xloc1, yloc1)
    distloc = ( ((xloc - x) * (xloc - x)) + ((yloc - y) * (yloc - y)) )**.5

    indloc = (distloc <= radius).nonzero()
    ycirc = indloc[0] - radius + y
    xcirc = indloc[1] - radius + x

    return (ycirc, xcirc)


"""
Find all indices within the local circle of radius but remove indeces that
are out of an area box, specified with an 2d array.
Input:
x: x index of center point
y: y index of center point
radius: radius in pixels, floats are handled including the farthest point
Returns a tuple of (y index, x index)
"""
def draw_cut_circle(x, y, radius, array):

    ycirc, xcirc = draw_circle(x, y, radius)
    noky = np.where(ycirc >= array.shape[0])  # if the circle is off the edge
    if noky[0].size > 0:
        ycirc = np.delete(ycirc, noky)
        xcirc = np.delete(xcirc, noky)

    nokx = np.where(xcirc >= array.shape[1])
    if nokx[0].size > 0:
        ycirc = np.delete(ycirc, nokx)
        xcirc = np.delete(xcirc, nokx)

    return (ycirc, xcirc)

def draw_ring(x, y, inrad, outrad, array):

    in_ycirc, in_xcirc = draw_cut_circle(x, y, inrad, array)
    out_ycirc, out_xcirc = draw_cut_circle(x, y, outrad, array)

    in_uni=unique_of_pair(in_xcirc, in_ycirc)
    out_uni = unique_of_pair(out_xcirc, out_ycirc)

    inter = np.in1d(out_uni, in_uni, assume_unique=True)

    if np.sum(inter) != 0:
        nok = np.where(inter)
        out_ycirc = np.delete(out_ycirc, nok)
        out_xcirc = np.delete(out_xcirc, nok)

    return (out_ycirc, out_xcirc)







