import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
import astropy.units as u
from scipy.optimize import curve_fit
import warnings
warnings.filterwarnings('ignore')

def king_profile(r, f0, rc, bg):
    return f0 / (1 + (r/rc)**2) + bg

def analyze_cluster(name, members_file, ax):
    df = pd.read_csv(members_file)
    ra_cen = df['ra'].median()
    dec_cen = df['dec'].median()
    center = SkyCoord(ra=ra_cen*u.deg, dec=dec_cen*u.deg)
    stars = SkyCoord(ra=df['ra'].values*u.deg, dec=df['dec'].values*u.deg)
    
    seps = center.separation(stars).arcminute
    max_r = np.max(seps)
    
    # Dynamic binning for sparse clusters
    n_bins = max(6, int(np.sqrt(len(df))))
    bins = np.linspace(0, max_r, n_bins)
    r_centers = 0.5 * (bins[1:] + bins[:-1])
    
    counts, _ = np.histogram(seps, bins=bins)
    areas = np.pi * (bins[1:]**2 - bins[:-1]**2)
    density = counts / areas
    
    try:
        popt, _ = curve_fit(king_profile, r_centers, density, p0=[np.max(density), np.median(seps), 0.0], bounds=(0, np.inf))
        f0, rc, bg = popt
        r_fit = np.linspace(0, max_r, 100)
        ax.plot(r_fit, king_profile(r_fit, *popt), 'r-', lw=2.5, label=f'King Fit (rc={rc:.2f}\')')
    except Exception as e:
        print(f"Fit failed for {name}: {e}")
        pass
        
    ax.bar(r_centers, density, width=(bins[1]-bins[0])*0.8, color='silver', edgecolor='black', alpha=0.6, label='Observed Density', zorder=1)
    ax.scatter(r_centers, density, color='black', s=50, zorder=2)
    
    ax.set_xlabel('Distance from Center (arcmin)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Surface Density (stars / arcmin2$)', fontsize=12, fontweight='bold')
    ax.set_title(f'{name} Spatial Profile (N={len(df)})', fontsize=14, fontweight='bold')
    ax.grid(True, linestyle='--', alpha=0.5)
    ax.legend(fontsize=11)
    
fig, axes = plt.subplots(1, 2, figsize=(14, 6))
analyze_cluster('Berkeley 87', r'C:\Users\raffe\Desktop\Dissertation\GaussCluster_Software\src\v8\Results\Be87\Be87_GaussCluster_Members.csv', axes[0])
analyze_cluster('NGC 957', r'C:\Users\raffe\Desktop\Dissertation\GaussCluster_Software\src\v8\Results\NGC957\NGC957_GaussCluster_Members.csv', axes[1])

plt.tight_layout()
out_plot = r'C:\Users\raffe\Desktop\Dissertation\Appendix\King_Profiles.png'
plt.savefig(out_plot, dpi=300)
print(f"Plot saved to {out_plot}")
