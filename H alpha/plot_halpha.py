import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
from astropy.coordinates import SkyCoord
import astropy.units as u

OUT_DIR = r"C:\Users\raffe\Desktop\Dissertation\Appendix\h alpha"
MEMBERS_PATH = r"C:\Users\raffe\Desktop\Dissertation\Appendix\h alpha\Be87_GaussCluster_Members.csv"
HALPHA_CSV = r"C:\Users\raffe\Desktop\Dissertation\Observations\Observation 8 H Alpha Be87\Be87_halpha_results.csv"

df_clean = pd.read_csv(MEMBERS_PATH)
df_halpha = pd.read_csv(HALPHA_CSV)

# Match spatially between Members and Halpha
c_clean = SkyCoord(ra=df_clean['ra'].values*u.deg, dec=df_clean['dec'].values*u.deg)
c_halpha = SkyCoord(ra=df_halpha['gaia_ra'].values*u.deg, dec=df_halpha['gaia_dec'].values*u.deg)
idx, d2d, _ = c_clean.match_to_catalog_sky(c_halpha)
mask = d2d < 1*u.arcsec

# Build matched dataframe
df_matched = df_clean[mask].copy()
df_matched['Halpha_inst_mag'] = df_halpha.iloc[idx[mask]]['Halpha_inst_mag'].values

sns.set_theme(style="white", context="paper")
fig, ax = plt.subplots(figsize=(7, 7))

rp = df_matched['phot_rp_mean_mag']
ha = df_matched['Halpha_inst_mag']
mask_ha = df_matched['Halpha_inst_mag'].notna()
med_diff = np.nanmedian((rp - ha)[mask_ha]) if mask_ha.any() else 0.0

if not np.isnan(med_diff) and not mask_ha.empty:
    min_val = np.nanmin([np.nanmin(rp), np.nanmin(ha + med_diff)]) - 0.5
    max_val = np.nanmax([np.nanmax(rp), np.nanmax(ha + med_diff)]) + 0.5

    x_line = np.linspace(min_val, max_val, 100)
    y_line = x_line - med_diff

    ax.plot(x_line, y_line, 'r--', lw=2.5, zorder=1, label="1:1 Emission Threshold")
    ax.scatter(rp, ha, color='darkblue', s=60, edgecolor='black', linewidth=0.8, alpha=0.9, zorder=2, label="3D Members")
    
    # -------------------------------------------------------------
    # MATCH V439 CYG SPATIALLY (RA = 305.313998, Dec = 37.408612)
    # -------------------------------------------------------------
    v439_coord = SkyCoord(ra=[305.313998]*u.deg, dec=[37.408612]*u.deg)
    c_matched = SkyCoord(ra=df_matched['ra'].values*u.deg, dec=df_matched['dec'].values*u.deg)
    idx_v439, d2d_v439, _ = v439_coord.match_to_catalog_sky(c_matched)
    
    if d2d_v439[0].arcsec < 1.0:
        # We found it spatially!
        v439_rp = df_matched.iloc[idx_v439[0]]['phot_rp_mean_mag']
        v439_ha = df_matched.iloc[idx_v439[0]]['Halpha_inst_mag']
        ax.scatter(v439_rp, v439_ha, color='gold', marker='*', s=350, edgecolor='black', linewidth=1.2, zorder=4, label="V439 Cyg (Spatial Match)")
        ax.annotate("V439 Cyg", (v439_rp, v439_ha), xytext=(12, 12), textcoords="offset points", fontsize=12, fontweight='bold', color='black', bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="black", lw=1, alpha=0.8))

    ax.set_xlim(min_val, max_val)
    ax.set_ylim(min_val - med_diff, max_val - med_diff)
    ax.set_aspect('equal', 'box')

ax.set_xlabel("Gaia RP (mag)", fontsize=14, fontweight='bold')
ax.set_ylabel("Instrumental H-alpha (mag)", fontsize=14, fontweight='bold')
ax.set_title("Berkeley 87: H-alpha Photometry", fontsize=16, fontweight='bold', pad=10)
for spine in ax.spines.values():
    spine.set_edgecolor('black')
    spine.set_linewidth(1.5)
ax.tick_params(axis='both', labelsize=12)
ax.legend(fontsize=12, loc="upper left", edgecolor='black')

plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "Be87_Halpha_TCD_Single_3D.png"), dpi=300)
print("Saved single H-alpha TCD with strictly spatial V439 Cyg labeling")
