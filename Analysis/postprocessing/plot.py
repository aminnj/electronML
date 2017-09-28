from __future__ import print_function


import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_auc_score
from sklearn.metrics import roc_curve

import plottery.plottery as ply
import plottery.utils as plu

data = np.load("todump.npa")
y_test = data[:,0]
y_pred = data[:,1]
pt_test = data[:,2]
mva_test = data[:,3]

print(y_test)
print(y_pred)
print("AUC total",roc_auc_score(y_test,y_pred))

def culled_indices(thresholds, precision=0.002):
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
        (30., 80., "CMS: 30 < p_{T} < 80"),
        (80., 1000., "CMS: 80 < p_{T} < 1000"),
        ]:
    y_test_sub = y_test[(pt_test > ptlow) & (pt_test < pthigh)]
    y_pred_sub = mva_test[(pt_test > ptlow) & (pt_test < pthigh)]
    fpr, tpr, thresholds = roc_curve(y_test_sub, y_pred_sub) #, pos_label=2)
    auc = np.trapz(tpr,fpr)
    to_keep = culled_indices(thresholds)
    fpr = fpr[to_keep]
    tpr = tpr[to_keep]
    vals.append( (fpr,tpr) )
    ptlabel = "{} [{:.3f}]".format(ptlabel,auc)
    legend_labels.append(ptlabel)

for ptlow,pthigh,ptlabel in [
        (0., 10., "Me: 0 < p_{T} < 10"),
        (10., 30., "Me: 10 < p_{T} < 30"),
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
            "xaxis_log": True,
            "yaxis_log": True,
            "xaxis_range": [0.005,0.8],
            "yaxis_range": [0.7,1.0],
            "output_name": "plots/roc.pdf",
            "output_ic": True,
        }
)

