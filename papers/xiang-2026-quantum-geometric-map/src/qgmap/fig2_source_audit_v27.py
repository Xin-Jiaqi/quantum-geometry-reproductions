"""Independent PDF-curve audit of Fig. 2(b) (Xiang et al. PRB 113 L201406).

Source points are *approximate pixel-extracted* observations, NOT authors' raw data.
No PDF image is included in the public output package. Compare normalized shapes
at 20/40 K with a parameter-free Fermi tanh shape and literal SI amplitude.
"""
from pathlib import Path
import csv
import numpy as np
from scipy.optimize import curve_fit
from scipy import constants as c
from .si_fig2_v25 import FigureInputs,evaluate

SOURCE_CALIBRATION={
  'pdf_page_1indexed':4,
  'render_dpi_equivalent':576,
  'pdf_crop_points':[175,51,279,160],
  'plot_x_center_px':438,
  'plot_x_left_px':97,
  'plot_x_right_px':779,
  'plot_y_zero_px':450,
  'plot_y_upper_px':110,
  'plot_y_lower_px':792,
  'pixel_selection':'color-filtered curve medians, positive branch, 454..535 px, sampled every 3 px',
}

def px_to_values(x,y):
    a=SOURCE_CALIBRATION
    return (np.asarray(x)-a['plot_x_center_px'])*80/(a['plot_x_right_px']-a['plot_x_left_px']), \
           (a['plot_y_zero_px']-np.asarray(y))*4/(a['plot_y_lower_px']-a['plot_y_upper_px'])

def thermal_shape(mu_meV,temp_K):
    mu=np.asarray(mu_meV,float)
    if temp_K<=0:raise ValueError('for fitted finite-temperature curves T must be positive')
    return np.tanh(mu*1.e-3*c.e/(2*c.k*temp_K))

def load_points(csvpath):
    with Path(csvpath).open(encoding='utf-8',newline='') as handle:
        rows=list(csv.DictReader(handle))
    return {k:[(float(r['mu_meV']),float(r['ordinate'])) for r in rows if r['curve']==k]
            for k in ('T20','T40')}

def fit_branch(pts):
    pts=np.asarray(pts,float); mu,y=pts.T
    def model(mu,A,T):return A*thermal_shape(mu,T)
    pars,cov=curve_fit(model,mu,y,p0=[1.75,30],bounds=([1.,5.],[2.5,80.]),maxfev=10000)
    dif=y-model(mu,*pars)
    return dict(amplitude=float(pars[0]),effective_T_K=float(pars[1]),n_points=len(pts),
                rms_ordinate=float(np.sqrt(np.mean(dif**2))),max_abs_ordinate=float(np.max(np.abs(dif))))

def sample_pdf_crop(pdf_file,out_csv):
    """Optional offline PDF-only step; PyMuPDF is not a runtime dependency for frozen tests.

    Uses published PDF on user's local disk. Pixels are approximate and kept as
    machine-readable points, not embedded publisher images.
    """
    import fitz
    page=fitz.open(pdf_file)[3]
    pix=page.get_pixmap(matrix=fitz.Matrix(8,8),clip=fitz.Rect(175,51,279,160),alpha=False)
    from PIL import Image
    import io
    rgb=np.asarray(Image.open(io.BytesIO(pix.tobytes('png'))).convert('RGB')).astype(np.float64)
    red,green,blue=rgb[:,:,0],rgb[:,:,1],rgb[:,:,2]
    masks={
        'T20':(green>75)&(green>red*1.4)&(green>blue*1.4)&(red<100),
        'T40':(blue>150)&(blue>red*1.4)&(blue>green*1.4)&(red<90),
    }
    records=[]
    for name,mask in masks.items():
        for x in range(454,536,3):
            yy=np.flatnonzero(mask[:,x] & (np.arange(mask.shape[0])>137) & (np.arange(mask.shape[0])<455))
            if len(yy)>2:
                y=float(np.median(yy));mu,val=px_to_values(float(x),y)
                records.append((name,int(x),y,float(mu),float(val)))
    with open(out_csv,'w',newline='',encoding='utf-8') as f:
        w=csv.writer(f);w.writerow(['curve','pixel_x','pixel_y','mu_meV','ordinate']);w.writerows(records)
    return records

def ratio_to_literal(plateau):
    paper_lit=evaluate(FigureInputs(),mu_mev=20,temp_K=0)['minus_1e4_sigma_over_e2_hbar']
    return plateau/paper_lit, paper_lit

def calibration_sensitivity(csvpath):
    """Uncertainty envelope from +/-0.5-pixel read-off of the two central axes."""
    with Path(csvpath).open(encoding='utf-8',newline='') as handle:
        rows=list(csv.DictReader(handle))
    res={}
    for label in ('T20','T40'):
        selected=[r for r in rows if r['curve']==label]
        x=np.array([float(r['pixel_x']) for r in selected]); y=np.array([float(r['pixel_y']) for r in selected])
        models=[]
        for dx in (-.5,0,.5):
            for dy in (-.5,0,.5):
                mu=(x-438-dx)*80/682
                ordinate=(450+dy-y)*4/682
                models.append(fit_branch(np.column_stack([mu,ordinate])))
        res[label]={
            'amplitude_range':[min(r['amplitude'] for r in models),max(r['amplitude'] for r in models)],
            'T_K_range':[min(r['effective_T_K'] for r in models),max(r['effective_T_K'] for r in models)],
            'calibration_model':'plot origin +/-0.5 raster pixels in x and y, all other anchors fixed',
        }
    return res
