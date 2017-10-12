from __future__ import print_function


import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_auc_score
from sklearn.metrics import roc_curve
from sklearn.ensemble import (RandomTreesEmbedding, RandomForestClassifier, GradientBoostingClassifier, AdaBoostClassifier)

import plottery.plottery as ply
import plottery.utils as plu
import xgboost as xgb
import sys

data = np.load("todump.npa")
y_test = data[:,0]
y_pred = data[:,1]
pt_test = data[:,2]
mva_test = data[:,3]
weights_test = data[:,4]
sieie_test = data[:,5]
xgb_test = data[:,6]

# y_test = y_test[sieie_test < 0.011]
# y_pred = y_pred[sieie_test < 0.011]
# pt_test = pt_test[sieie_test < 0.011]
# mva_test = mva_test[sieie_test < 0.011]
# weights_test = weights_test[sieie_test < 0.011]
# sieie_test = sieie_test[sieie_test < 0.011]

use_sieie_as_disc = False
use_xgb_as_disc = True

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
for ptlow,pthigh,ptlabel in [
        (0., 10., "CMS: 0 < p_{T} < 10"),
        (10., 30., "CMS: 10 < p_{T} < 30"),
        # (21., 22., "CMS: 21 < p_{T} < 22"),
        (30., 80., "CMS: 30 < p_{T} < 80"),
        (80., 1000., "CMS: 80 < p_{T} < 1000"),
        ]:
    y_test_sub = y_test[(pt_test > ptlow) & (pt_test < pthigh)]
    if use_sieie_as_disc:
        y_pred_sub = 1.0-sieie_test[(pt_test > ptlow) & (pt_test < pthigh)]
    elif use_xgb_as_disc:
        y_pred_sub = xgb_test[(pt_test > ptlow) & (pt_test < pthigh)]
    else:
        y_pred_sub = mva_test[(pt_test > ptlow) & (pt_test < pthigh)]
    fpr, tpr, thresholds = roc_curve(y_test_sub, y_pred_sub) #, pos_label=2)
    # if ptlow == 10.:
    #     for bkgeff,sigeff,thresh in zip(fpr,tpr,thresholds):
    #         if 0.985 < thresh < 0.990:
    #             print(bkgeff,sigeff,thresh)
    auc = np.trapz(tpr,fpr)
    print(auc)
    to_keep = culled_indices(thresholds)
    fpr = fpr[to_keep]
    tpr = tpr[to_keep]
    vals.append( (fpr,tpr) )
    ptlabel = "{} [{:.3f}]".format(ptlabel,auc)
    legend_labels.append(ptlabel)

for ptlow,pthigh,ptlabel in [
        (0., 10., "Me: 0 < p_{T} < 10"),
        (10., 30., "Me: 10 < p_{T} < 30"),
        # (21., 22., "Me: 21 < p_{T} < 22"),
        (30., 80., "Me: 30 < p_{T} < 80"),
        (80., 1000., "Me: 80 < p_{T} < 1000"),
        ]:
    y_test_sub = y_test[(pt_test > ptlow) & (pt_test < pthigh)]
    y_pred_sub = y_pred[(pt_test > ptlow) & (pt_test < pthigh)]
    fpr, tpr, thresholds = roc_curve(y_test_sub, y_pred_sub)
    auc = np.trapz(tpr,fpr)
    to_keep = culled_indices(thresholds)
    fpr = fpr[to_keep]
    tpr = tpr[to_keep]
    vals.append( (fpr,tpr) )
    ptlabel = "{} [{:.3f}]".format(ptlabel,auc)
    legend_labels.append(ptlabel)
colors = []
colors.extend(plu.interpolate_colors_rgb( (1.,0.,0.),(0.,0.,1.), 4 ))
colors.extend(plu.interpolate_colors_rgb( (1.,0.,0.),(0.,0.,1.), 4 ))
draw_styles = [0]*4 + [2]*4

ply.plot_graph(
        vals,
        colors = colors,
        legend_labels = legend_labels,
        draw_styles = draw_styles,
        options = {
            "legend_alignment": "bottom right",
            "xaxis_label": "bkg. eff.",
            "yaxis_label": "sig. eff.",
            "title": "DY: prompt vs unmatched e & b/c#rightarrow e",
            "xaxis_log": False,
            "yaxis_log": False,
            # "xaxis_range": [0.005,0.8],
            # "yaxis_range": [0.7,1.0],
            "xaxis_range": [0.005,0.8],
            "yaxis_range": [0.2,1.0],
            "output_name": "plots/roc.pdf",
            "output_ic": True,
        }
)

