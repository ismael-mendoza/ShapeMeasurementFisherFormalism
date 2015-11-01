import os
os.chdir("/Users/Ismael/code/research/repo/")

import fisher
import galfun
import galsim
import sys
import copy
import numpy as np
import math
import os
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import pylab

def shearEllipticity(g, e):
    """Changes given ellipticity to a sheared ellpticity according to the
    formula (14) of paper: http://arxiv.org/abs/1409.6273.
    Ellipticity should be given by inputting the two components.
    Returns both sheared components of ellpicity.
    Uses complex numbers.
    """
    e_s = (e + g)/(1 + g.conjugate()*e)
    return e_s

def shearBias(fish, g):
    """Returns the value of the bias of the given lensing shear for the
    particular galaxy analyzed by using the ring test. Only works with
    a single galaxy profile.
    Assume galaxy is parametrized with e1,e2.
    """
    angle_range = (0, 2*math.pi)
    steps = 7 #6 points on the ring. excluding theta=2pi
    id_params = copy.deepcopy(fish.g_parameters.id_params)
    ids = id_params.keys()
    id1 = ids[0]
    snr = fish.snr
    angles = list(np.linspace(angle_range[0], angle_range[1], steps))
    angles.pop() #remove 2pi angle redundancy.
    biases = []
    orig_e1 = id_params[id1]['e1']
    orig_e2 = id_params[id1]['e2']
    orig_e = complex(orig_e1,orig_e2) #unsheared ellipticity
    for angle in angles:
        e = complex(abs(orig_e)*math.cos(angle),abs(orig_e)*math.sin(angle))
        e_s = shearEllipticity(g, e) #get sheared components
        id_params[id1]['e1'] = e_s.real
        id_params[id1]['e2'] = e_s.imag
        g_parameters = galfun.GParameters(id_params=id_params)
        fish = fisher.Fisher(g_parameters, snr)
        bias = complex(fish.biases['e1_'+id1],fish.biases['e2_'+id1])
        biases.append(bias)

    #sanity check.
    #print np.mean(ellipticities_s)

    #return bias(g) which is average of ellipticity bias
    return np.mean(biases)

def biasesEllipticities(fish, g):
    """Returns a lists of the biases of both ellipticies that are in the process of being sheared.
    """
    angle_range = (0, 2*math.pi)
    steps = 7 #6 points on the ring. excluding theta=2pi
    id_params = copy.deepcopy(fish.g_parameters.id_params)
    ids = id_params.keys()
    id1 = ids[0]
    snr = fish.snr
    angles = list(np.linspace(angle_range[0], angle_range[1], steps))
    angles.pop() #remove 2pi angle redundancy.
    biases_e1 = []
    biases_e2 = []
    orig_e1 = id_params[id1]['e1']
    orig_e2 = id_params[id1]['e2']
    orig_e = complex(orig_e1,orig_e2) #unsheared ellipticity
    for angle in angles:
        e = complex(abs(orig_e)*math.cos(angle),abs(orig_e)*math.sin(angle))
        e_s = shearEllipticity(g, e) #get sheared components
        ellipticities_s.append(e_s)
        id_params[id1]['e1'] = e_s.real
        id_params[id1]['e2'] = e_s.imag
        g_parameters = galfun.GParameters(id_params=id_params)
        fish = fisher.Fisher(g_parameters, snr)
        biases_e1.append(fish.biases['e1_'+id1])
        biases_e2.append(fish.biases['e2_'+id1])


    #sanity check.
    #print np.mean(ellipticities_s)

    #return bias(g) which is average of ellipticity bias
    return biases_e1, biases_e2

#get multiplicative and additive biases of shear by doing a linear fit.
#we let g2 be a constant  = 0
def getMultAddBiasG1(fish):
    g_range = (.0001, .02)
    steps = 10 #number of points to plot.
    g_ranges = list(np.linspace(g_range[0], g_range[1], steps))
    gs = []
    bs = []
    g2 = 0
    for g1 in g_ranges:
        g = complex(g1,g2)
        b = shearBias(fish, g)
        gs.append(g)
        bs.append(b)
    gs_real = np.array([g.real for g in gs])
    gs_imag = np.array([g.imag for g in gs])
    bs_real = np.array([b.real for b in bs])
    bs_imag = np.array([b.imag for b in bs])
    m1, c1 = np.polyfit(gs_real, bs_real,1)
    m2, c2 = np.polyfit(gs_imag, bs_imag, 1)
    return (m1, c1, m2, c2)

#get multiplicative and additive biases of shear by doing a linear fit.
#we let g1 be a constant  = 0
def getMultAddBiasG2(fish):
    g_range = (.0001, .02)
    steps = 10 #number of points to plot.
    g_ranges = list(np.linspace(g_range[0], g_range[1], steps))
    gs = []
    bs = []
    g1 = 0
    for g2 in g_ranges:
        g = complex(g1,g2)
        b = shearBias(fish, g)
        gs.append(g)
        bs.append(b)
    gs_real = np.array([g.real for g in gs])
    gs_imag = np.array([g.imag for g in gs])
    bs_real = np.array([b.real for b in bs])
    bs_imag = np.array([b.imag for b in bs])

    m1, c1 = np.polyfit(gs_real, bs_real,1)
    m2, c2 = np.polyfit(gs_imag, bs_imag, 1)
    return (m1, c1, m2, c2)
