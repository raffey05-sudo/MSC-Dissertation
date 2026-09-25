import argparse
import os
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--gc', required=True, help="GaussCluster members CSV")
    parser.add_argument('--storm', required=True, help="STORM members CSV")
    parser.add_argument('--extname', default="STORM", help="External catalog name")
    parser.add_argument('--field', required=True, help="GaussCluster field CSV")
    parser.add_argument('--outimg', required=True, help="Path to save the output diagnostic image")
    parser.add_argument('--title', required=True, help="Cluster Title")
    parser.add_argument('--spatial-cut', type=float, default=15.0, help="Cutoff in arcmin from GC center")
    args = parser.parse_args()

    gc = pd.read_csv(args.gc)
    storm = pd.read_csv(args.storm)
    field = pd.read_csv(args.field)
    
    # Safely get variables for generic catalogs
    s_ra = 'RA_ICRS' if 'RA_ICRS' in storm.columns else 'ra'
    s_dec = 'DE_ICRS' if 'DE_ICRS' in storm.columns else 'dec'
    s_pmra = 'pmRA' if 'pmRA' in storm.columns else 'pmra'
    s_pmdec = 'pmDE' if 'pmDE' in storm.columns else 'pmdec'
    s_plx = 'Plx' if 'Plx' in storm.columns else 'parallax'
    s_g = 'Gmag' if 'Gmag' in storm.columns else 'phot_g_mean_mag'
    s_bp = 'BPmag' if 'BPmag' in storm.columns else 'phot_bp_mean_mag'
    s_rp = 'RPmag' if 'RPmag' in storm.columns else 'phot_rp_mean_mag'

    # Apply spatial cut if requested
    if args.spatial_cut is not None:
        # The GC catalog is ALREADY cut by GaussCluster.py, so we just use its center
        med_ra = gc['ra'].median()
        med_dec = gc['dec'].median()
        center = SkyCoord(ra=med_ra*u.deg, dec=med_dec*u.deg)
        
        # We DO NOT cut gc again here! Re-cutting it based on the median of the already-cut 
        # catalog causes edge-cases to drop out due to the center shifting slightly.
        
        # Cut External Catalog
        storm_coords = SkyCoord(ra=storm[s_ra].values*u.deg, dec=storm[s_dec].values*u.deg)
        storm_valid = center.separation(storm_coords).to(u.arcmin).value <= args.spatial_cut
        storm = storm[storm_valid]

        print(f"Applied {args.spatial_cut} arcmin cut. GC remaining: {len(gc)}, {args.extname} remaining: {len(storm)}")

    # Coordinate crossmatch
    gc_sc = SkyCoord(ra=gc['ra'].values*u.deg, dec=gc['dec'].values*u.deg)
    storm_sc = SkyCoord(ra=storm[s_ra].values*u.deg, dec=storm[s_dec].values*u.deg)
    idx, d2d, _ = gc_sc.match_to_catalog_sky(storm_sc)
    
    overlap_mask = d2d < 1.0 * u.arcsec
    overlap_count = overlap_mask.sum()
    
    n_gc = len(gc)
    n_storm = len(storm)
    
    print(f"[{args.title}] Overlap Analysis:")
    print(f"GaussCluster Members: {n_gc}")
    print(f"{args.extname} Members: {n_storm}")
    print(f"Intersection: {overlap_count}")

    # Set up the 3-panel plot
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))

    # Panel 0: Spatial Distribution
    # axes[0].scatter(field['ra'], field['dec'], s=1, color='grey', alpha=0.05)
    axes[0].scatter(storm[s_ra], storm[s_dec], facecolors='none', edgecolors='red', s=40, alpha=0.5, label=f'{args.extname} ({n_storm})')
    axes[0].scatter(gc['ra'], gc['dec'], s=15, marker='+', color='blue', label=f'GaussCluster ({n_gc})')
    axes[0].set_xlabel('RA (deg)', fontsize=15)
    axes[0].set_ylabel('Dec (deg)', fontsize=15)
    # axes[0].set_xlim(field['ra'].min(), field['ra'].max())
    # axes[0].set_ylim(field['dec'].min(), field['dec'].max())
    axes[0].tick_params(axis='both', which='major', labelsize=12)
    axes[0].legend(fontsize=12, loc='best')
    axes[0].set_title('Spatial Distribution', fontsize=20)

    # Panel 1: Parallax Distribution
    field_plx = field['parallax_corrected'] if 'parallax_corrected' in field.columns else field['parallax']
    gc_plx = gc['parallax_corrected'] if 'parallax_corrected' in gc.columns else gc['parallax']
    storm_plx = storm[s_plx]
    
    # Range centered on the actual plx being plotted (gc_plx)
    plx_mean = gc_plx.mean()
    axes[1].hist(field_plx.dropna(), bins=50, alpha=0.3, color='grey', label='Field', density=True, range=(plx_mean - 3, plx_mean + 3))
    axes[1].hist(storm_plx.dropna(), bins=15, alpha=0.5, color='red', label=args.extname, density=True, histtype='step', linewidth=2)
    axes[1].hist(gc_plx.dropna(), bins=20, alpha=0.5, color='blue', label='GaussCluster', density=True)
    axes[1].set_xlim(plx_mean - 1.5, plx_mean + 1.5)
    axes[1].set_xlabel('Parallax (mas)', fontsize=15)
    axes[1].set_ylabel('Density', fontsize=15)
    axes[1].tick_params(axis='both', which='major', labelsize=12)
    axes[1].legend(fontsize=12, loc='best')
    axes[1].set_title('Parallax Distribution', fontsize=20)

    # Panel 2: CMD
    storm_bp_rp = storm[s_bp] - storm[s_rp]
    axes[2].scatter(storm_bp_rp, storm[s_g], facecolors='none', edgecolors='red', s=40, alpha=0.5, label=args.extname)
    
    # GaussCluster CMD using phot_quality_ok if available
    gc_cmd = gc[gc['phot_quality_ok']] if 'phot_quality_ok' in gc.columns else gc
    axes[2].scatter(gc_cmd['bp_rp'], gc_cmd['phot_g_mean_mag'], s=15, marker='+', color='blue', label='GaussCluster')
    axes[2].invert_yaxis()
    axes[2].set_xlabel('BP-RP (mag)', fontsize=15)
    axes[2].set_ylabel('G (mag)', fontsize=15)
    axes[2].tick_params(axis='both', which='major', labelsize=12)
    axes[2].legend(fontsize=12, loc='best')
    axes[2].set_title('Colour-Magnitude Diagram', fontsize=20)

    plt.tight_layout()
    os.makedirs(os.path.dirname(os.path.abspath(args.outimg)), exist_ok=True)
    plt.savefig(args.outimg, dpi=150)
    plt.close()

if __name__ == "__main__":
    main()



