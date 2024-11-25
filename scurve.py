import ROOT
import yaml
import numpy as np
import argparse
import matplotlib.pyplot as plt
import parameters
import math
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.colors import Normalize
from omegaconf import OmegaConf
import csv
import os
import mplhep as hep

parser = argparse.ArgumentParser(description='Plot a scurve map from a given root file')
parser.add_argument('cfg', help='YAML file with all the analysis parameters', type=str)
args = parser.parse_args()

def get_rootpath(hybrid, chip):
    return f'Detector/Board_0/OpticalGroup_0/Hybrid_{hybrid}/Chip_{chip}'

def get_rootprefix(hybrid):
    return f'D_B(0)_O(0)_H({hybrid})'

def TH2F_to_matrix(hist):
    nx, ny = hist.GetNbinsX(), hist.GetNbinsY()
    matrix = np.ndarray((nx, ny))
    for i in range(nx):
        for j in range(ny):
            matrix[i,j] = hist.GetBinContent(i+1, j+1)
    return matrix

def TH1F_to_matrix(hist):
    nx = hist.GetNbinsX() + 2
    xdata = np.ndarray(nx)
    ydata = np.ndarray(nx)
    for i in range(nx):
        xdata[i] = hist.GetXaxis().GetBinCenter(i)
        ydata[i] = hist.GetBinContent(i)
    return np.vstack((xdata, ydata))

# takes the difference between 2d plots and puts in a matrix the new 2d plot
def diff_TH2F_to_matrix(hist1, hist2):
    nx, ny = hist1.GetNbinsX(), hist1.GetNbinsY()
    matrix = np.ndarray((nx, ny))
    for i in range(nx):
        for j in range(ny):
            matrix[i,j] = hist1.GetBinContent(i+1, j+1) - hist2.GetBinContent(i+1, j+1)
    return matrix


def get_hist(conf,filename,hybrid,prefix,scan_list):
    hist = []
    for root_scan in scan_list:
        for chip in conf.ph2acf.chip_id:
            path = get_rootpath(hybrid,chip)
            try:
                get_hist = filename.Get(f"{path}/{prefix}_" + root_scan + f"({chip})").GetPrimitive(f"{prefix}_" + root_scan + f"({chip})")
                matrix = TH2F_to_matrix(get_hist)
            except:
                print('Directory ', f"{path}/{prefix}_" + root_scan + f"({chip})", 'does not exist -> set to 0') 
                matrix = np.zeros((432,336)).astype('int')   
            hist.append(matrix)
    return hist

def quadMap(conf, occupancy_matrix):
    occupancyMap = []
    for chip in conf.ph2acf.chip_id:
        if chip in range(0,2):
            occupancy_matrix[chip] = np.rot90(occupancy_matrix[chip])
            occupancy_matrix[chip] = np.flip(occupancy_matrix[chip],1)
            occupancy_matrix[chip] = np.rot90(occupancy_matrix[chip])
        else:
            occupancy_matrix[chip] = np.flip(occupancy_matrix[chip],0)
        occupancyMap.append(occupancy_matrix[chip]) 
    hitmap = np.zeros((864,672)).astype('int')
    hitmap [432:864, 0:336]   = occupancyMap[0]
    hitmap [0:432, 0:336]     = occupancyMap[1]
    hitmap [0:432, 336:672]   = occupancyMap[2]
    hitmap [432:864, 336:672] = occupancyMap[3]
    hitmap = np.rot90(hitmap)
    hitmap = np.flip(hitmap,0)
    return hitmap
    
def dualMap(conf, occupancy_matrix):
    hitmap = np.zeros((864,336)).astype('int')
    hitmap [432:864, 0:336] = occupancy_matrix[0]
    hitmap [0:432, 0:336]   = occupancy_matrix[1]
    hitmap = np.rot90(hitmap)
    return hitmap

def save_cfg(conf):
    out_folder = conf['output']['output_folder']+ 'config/'
    if not os.path.isdir(out_folder):
        os.makedirs(out_folder)
    out_yaml = out_folder+conf['output']['module_name']+".yaml"
    OmegaConf.save(config=conf, f=out_yaml)
    print('\n'+ 'Configuration file saved in', out_yaml + '\n')

################################################################################################################

base_conf = parameters.get_default_parameters()
second_conf = parameters.get_parameters(args.cfg)
conf = parameters.merge_parameters(base_conf, second_conf)

'''
Some plotting parameters from the config
'''
module        = conf['output']['module_name']
maxvcal_thr   = conf['param']['maxVcal_thr']
minvcal_thr   = conf['param']['minVcal_thr']
maxvcal_noise = conf['param']['maxVcal_noise']
minvcal_noise = conf['param']['minVcal_noise']
outfolder     = conf['output']['output_folder']
dual          = conf['ph2acf']['is_dual']
chip_n = len(conf.ph2acf.chip_id)

if not os.path.isdir(outfolder):
    os.makedirs(outfolder)

#print(chip_n)

hybrid = conf['ph2acf']['hybrid']
prefix = get_rootprefix(hybrid)
f = ROOT.TFile.Open(conf['input']['input_folder']+conf['input']['input_file'])

scans_list = [conf['scans']['threshold2d'],conf['scans']['noise2d']]

thr2d   = []
noise2d = []
hist    = []

hist = get_hist(conf, f, hybrid, prefix, scans_list)

#print(len(hist))

for chip_id in conf.ph2acf.chip_id:
    thr2d_chip   = hist[chip_id]
    noise2d_chip = hist[chip_id+chip_n]
    thr2d.append(thr2d_chip)
    noise2d.append(noise2d_chip)

if dual == True:
    figsize_min = 13
    figsize_max = 6
    thrMap   = dualMap(conf, thr2d)
    noiseMap = dualMap(conf, noise2d)
    col_labels = ["0","100","200","300","432","100","200","300","432"]
    row_labels = ["0","100","200","336"]
    col_x = np.array([0,100,200,300,432,532,632,732,864])
    row_y = np.array([0,100,200,336])

else:
    figsize_min = 13
    figsize_max = 10
    thrMap   = quadMap(conf, thr2d)
    noiseMap = quadMap(conf, noise2d)
    col_labels = ["0","100","200","300","432","100","200","300","432"]
    row_labels = ["0","100","200","336","100","200","336"]
    col_x = np.array([0,100,200,300,432,532,632,732,864])
    row_y = np.array([0,100,200,336,436,536,672])
    roc23_labels = ["ROC2","ROC3"]

rocs  = np.array([216,648])
roc01_labels = ["ROC1","ROC0"]

'''
THRESHOLD MAP 
2d histogram containing the threshold distributions 
'''
plt.style.use([hep.cms.style.ROOT])
f0,ax0 = plt.subplots(figsize=(figsize_min,figsize_max))
f0.tight_layout(pad=3)
a0 = ax0.pcolor(thrMap, cmap=plt.cm.viridis, vmin=0, vmax=500)
ax0.axvline(x=432, color='black', linestyle='--', linewidth=2)
cbar = f0.colorbar(a0, ax=ax0)
cbar.set_label('$\Delta$VCAL[VCAL]', labelpad=20)
ax0.set_ylabel('row')
ax0.set_xlabel('column')
ax0.set_xticks(ticks=col_x, labels=col_labels)
ax0.set_yticks(ticks=row_y, labels=row_labels)
ax0.tick_params(axis='x', labelsize = 16)
ax0.tick_params(axis='y', labelsize = 16)
secax = ax0.secondary_xaxis(0)
secax.set_xticks(ticks=rocs, labels=roc01_labels)
secax.tick_params(axis='x', pad=30)
if dual == False:
    ax0.axvline(y=336, color='black', linestyle='--', linewidth=2)
    thirdax = ax0.secondary_xaxis(1)
    thirdax.set_xticks(ticks=rocs, labels=roc23_labels)
ax0.set_title('2D Threshold Map', y=1.02)
f0.savefig(outfolder + '2DThreshold_' + module, dpi=300)
plt.show()


'''
NOISE MAP 
2d histogram containing the noise distributions 
'''
f1,ax1 = plt.subplots(figsize=(figsize_min,figsize_max))
f1.tight_layout(pad=3)
a1 = ax1.pcolor(noiseMap, cmap=plt.cm.viridis, vmin=0, vmax=50)
ax1.axvline(x=432, color='black', linestyle='--', linewidth=2)
cbar = f1.colorbar(a1, ax=ax1)
cbar.set_label('$\Delta$ENC', labelpad=20)
ax1.set_ylabel('row')
ax1.set_xlabel('column')
ax1.set_xticks(ticks=col_x, labels=col_labels)
ax1.set_yticks(ticks=row_y, labels=row_labels)
ax1.tick_params(axis='x', labelsize = 16)
ax1.tick_params(axis='y', labelsize = 16)
secax = ax1.secondary_xaxis(0)
secax.set_xticks(ticks=rocs, labels=roc01_labels)
secax.tick_params(axis='x', pad=30)
if dual == False:
    ax1.axvline(y=336, color='black', linestyle='--', linewidth=2)
    thirdax = ax1.secondary_xaxis(1)
    thirdax.set_xticks(ticks=rocs, labels=roc23_labels)
ax1.set_title('2D Noise Map', y=1.02)
f1.savefig(outfolder + '2DNoise_' + module, dpi=300)
plt.show()

'''
THRESHOLD HIST 
1d histogram containing the threshold distributions 
'''
f2,ax2 = plt.subplots(figsize=(13,10))
a2 = ax2.hist(thrMap.flatten(), bins=(maxvcal_thr-minvcal_thr), range=(minvcal_thr, maxvcal_thr),density=False, histtype='step', color = 'blue')
ax2.set_xlabel('$\Delta$VCAL[VCAL]')
ax2.set_ylabel('counts')
ax2.set_title('1D Threshold Histogram', y=1.02)
f2.savefig(outfolder + '1DThreshold_' + module, dpi=300)
plt.show()

'''
NOISE HIST 
1d histogram containing the noise distributions 
'''
f3,ax3 = plt.subplots(figsize=(13,10))
a3 = ax3.hist(noiseMap.flatten(), bins=(maxvcal_noise-minvcal_noise), range=(minvcal_noise, maxvcal_noise),density=False, histtype='step', color = 'red')
ax3.set_xlabel('$\Delta$ENC')
ax3.set_ylabel('counts')
ax3.set_title('1D Noise Histogram', y=1.02)
f3.savefig(outfolder + '1DNoise_' + module, dpi=300)
plt.show()