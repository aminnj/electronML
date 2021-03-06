from __future__ import print_function


import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_auc_score
from sklearn.metrics import roc_curve
# from sklearn.ensemble import (RandomTreesEmbedding, RandomForestClassifier, GradientBoostingClassifier, AdaBoostClassifier)

import plottery.plottery as ply
import plottery.utils as plu
# import xgboost as xgb
import sys

# data = np.load("todump.npa")
data = np.load("todump_prime.npa")
# data = np.load("todump_prime_test.npa")
y_test = data[:,0]
y_pred = data[:,1]
pt_test = data[:,2]
mva_test = data[:,3]
weights_test = data[:,4]
sieie_test = data[:,5]
xgb_test = data[:,6]
xgbprime_test = data[:,7]
fname = "plots/roc.pdf"

# y_test = y_test[sieie_test < 0.011]
# y_pred = y_pred[sieie_test < 0.011]
# pt_test = pt_test[sieie_test < 0.011]
# mva_test = mva_test[sieie_test < 0.011]
# weights_test = weights_test[sieie_test < 0.011]
# sieie_test = sieie_test[sieie_test < 0.011]

use_sieie_as_disc = False
use_xgb_as_disc = True
use_xgbprime_as_disc = False

# grd = GradientBoostingClassifier(n_estimators=3,verbose=1, max_depth=1,loss="exponential") 
# print(grd)
# x_train = np.c_[sieie_test]
# x_test = x_train
# grd.fit(x_test,y_test)
# sieie_disc = grd.predict_proba(x_test)[:,1]
# print(sieie_disc)



# sys.exit()

# print(y_pred)
# print(mva_test)

print(y_test)
print(y_pred)
print("AUC total",roc_auc_score(y_test,y_pred))
print("AUC total",roc_auc_score(y_test,xgb_test))
print("AUC total",roc_auc_score(y_test,xgbprime_test))

def culled_indices(thresholds, precision=0.0002):
    """
    input a list of discriminant values and a precision
    and this returns a list of indices on the values for which
    the difference wrt previous values is >= precision
    """
    prev = -999.
    to_keep = []
    for ithr,thr in enumerate(thresholds):
        diff = abs(thr - prev)
        if diff > precision:
            to_keep.append(ithr)
            prev = thr
        else:
            continue
    return np.array(to_keep)

vals = []
legend_labels = []
tot = 0.
bininfo = [
        (0., 10., "%s: 0 < p_{T} < 10"),
        (10., 30., "%s: 10 < p_{T} < 30"),
        (30., 80., "%s: 30 < p_{T} < 80"),
        (80., 1000., "%s: 80 < p_{T} < 1e3"),
        ]
for ptlow,pthigh,ptlabel in bininfo:
    y_test_sub = y_test[(pt_test > ptlow) & (pt_test < pthigh)]
    if use_sieie_as_disc:
        ptlabel = ptlabel % "#sigma_{i#etai#eta}"
        y_pred_sub = 1.0-sieie_test[(pt_test > ptlow) & (pt_test < pthigh)]
    elif use_xgb_as_disc:
        ptlabel = ptlabel % "XGB"
        y_pred_sub = xgb_test[(pt_test > ptlow) & (pt_test < pthigh)]
    elif use_xgbprime_as_disc:
        ptlabel = ptlabel % "XGBCNN"
        y_pred_sub = xgbprime_test[(pt_test > ptlow) & (pt_test < pthigh)]
    else:
        ptlabel = ptlabel % "CMS"
        y_pred_sub = mva_test[(pt_test > ptlow) & (pt_test < pthigh)]
    fpr, tpr, thresholds = roc_curve(y_test_sub, y_pred_sub) #, pos_label=2)
    # if ptlow == 10.:
    #     for bkgeff,sigeff,thresh in zip(fpr,tpr,thresholds):
    #         if 0.985 < thresh < 0.990:
    #             print(bkgeff,sigeff,thresh)
    auc = np.trapz(tpr,fpr)
    to_keep = culled_indices(thresholds)
    fpr = fpr[to_keep]
    tpr = tpr[to_keep]
    vals.append( (fpr,tpr) )
    ptlabel = "{} [{:.3f}]".format(ptlabel,auc)
    legend_labels.append(ptlabel)

for ptlow,pthigh,ptlabel in bininfo:
    ptlabel = ptlabel % "XGBCNN"
    y_test_sub = y_test[(pt_test > ptlow) & (pt_test < pthigh)]
    # y_pred_sub = y_pred[(pt_test > ptlow) & (pt_test < pthigh)]
    y_pred_sub = xgbprime_test[(pt_test > ptlow) & (pt_test < pthigh)]
    fpr, tpr, thresholds = roc_curve(y_test_sub, y_pred_sub)
    auc = np.trapz(tpr,fpr)
    to_keep = culled_indices(thresholds)
    fpr = fpr[to_keep]
    tpr = tpr[to_keep]
    vals.append( (fpr,tpr) )
    ptlabel = "{} [{:.3f}]".format(ptlabel,auc)
    legend_labels.append(ptlabel)

# # colors = []
# # colors.extend(plu.interpolate_colors_rgb( (1.,0.,0.),(0.,0.,1.), len(legend_labels)/2 ))
# # colors.extend(plu.interpolate_colors_rgb( (1.,0.,0.),(0.,0.,1.), len(legend_labels)/2 ))
# colors = []
# colors.extend(plu.get_brightdefault_colors()[:len(legend_labels)/2])
# colors.extend(plu.get_brightdefault_colors()[:len(legend_labels)/2])
# draw_styles = [0]*(len(legend_labels)/2) + [2]*(len(legend_labels)/2)
colors = []
colors.extend(plu.get_brightdefault_colors()[:len(legend_labels)/2])
colors.extend(plu.get_brightdefault_colors()[:len(legend_labels)/2])
draw_styles = [0]*(len(legend_labels)/2) + [2]*(len(legend_labels)/2)

temp = plu.get_brightdefault_colors()[:len(legend_labels)/2]
vals = sum([[v1,v2] for v1,v2 in zip(vals[:len(vals)/2],vals[len(vals)/2:])],[])
legend_labels = sum([[l1,l2] for l1,l2 in zip(legend_labels[:len(vals)/2],legend_labels[len(vals)/2:])],[])
colors = sum([[t,t] for t in temp],[])
draw_styles = [0,2]*(len(legend_labels)/2)

ply.plot_graph(
        vals,
        colors = colors,
        legend_labels = legend_labels,
        draw_styles = draw_styles,
        options = {
            "legend_alignment": "bottom right",
            "legend_scalex": 1.45,
            "legend_scaley": 1.5,
            "xaxis_label": "bkg. eff.",
            "yaxis_label": "sig. eff.",
            "title": "DY: prompt vs unmatched e & b/c#rightarrow e",
            "xaxis_log": True,
            "yaxis_log": True,
            "xaxis_range": [0.005,0.8],
            "yaxis_range": [0.7,1.0],
            # "xaxis_range": [0.005,0.8],
            # "yaxis_range": [0.5,1.0],
            # "output_name": "plots/roc.pdf",
            "output_name": fname,
            "output_ic": True,
        }
)

