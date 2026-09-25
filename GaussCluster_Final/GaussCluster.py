import os
import sys
import argparse
import warnings

import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
import astropy.units as u
from astroquery.gaia import Gaia
from scipy.stats import gaussian_kde, median_abs_deviation
import matplotlib.pyplot as plt

try:
    from zero_point import zpt
    ZPT_AVAILABLE = True
except ImportError:
    ZPT_AVAILABLE = False

GLOBAL_ZPT_FALLBACK_MAS = -0.017  # Lindegren+21 global mean offset (mas), fallback only


def sigma_clip_3d(x, y, z, n_sigma=3, max_iter=5, floor_x=0.05, floor_y=0.05, floor_z=0.02):
    """MAD-based, floored, iterative 3D sigma clip (robust to single outliers)."""
    mask = np.ones(len(x), dtype=bool)
    for _ in range(max_iter):
        if mask.sum() < 3:
            break
        mx, my, mz = np.median(x[mask]), np.median(y[mask]), np.median(z[mask])
        sx = max(median_abs_deviation(x[mask], scale='normal'), floor_x)
        sy = max(median_abs_deviation(y[mask], scale='normal'), floor_y)
        sz = max(median_abs_deviation(z[mask], scale='normal'), floor_z)
        new_mask = (np.abs(x - mx) < n_sigma * sx) & (np.abs(y - my) < n_sigma * sy) & (np.abs(z - mz) < n_sigma * sz)
        if (new_mask == mask).all():
            break
        mask = new_mask
    return mask


def fit_kinematics(seed_pmra, seed_pmdec, seed_plx):
    mask = sigma_clip_3d(seed_pmra, seed_pmdec, seed_plx)
    if mask.sum() < 3:
        raise RuntimeError("Fewer than 3 seed members survive sigma-clipping -- "
                            "check the seed catalogue / crossmatch radius.")
    mu = np.array([np.median(seed_pmra[mask]), np.median(seed_pmdec[mask]), np.median(seed_plx[mask])])
    sigma = np.array([
        max(np.std(seed_pmra[mask]), 0.05),
        max(np.std(seed_pmdec[mask]), 0.05),
        max(np.std(seed_plx[mask]), 0.02),
    ])
    return mu, sigma

# Data acquisition + zero-point correction (official Gaia TAP archive) #

def query_gaia_dr3(ra_cen, dec_cen, radius_deg, row_limit=-1):
    query = f"""
    SELECT source_id, ra, dec, pmra, pmra_error, pmdec, pmdec_error,
           parallax, parallax_error, phot_g_mean_mag, phot_bp_mean_mag,
           phot_rp_mean_mag, bp_rp, ruwe, phot_bp_rp_excess_factor,
           nu_eff_used_in_astrometry, pseudocolour, ecl_lat,
           astrometric_params_solved,
           pmra_pmdec_corr, parallax_pmra_corr, parallax_pmdec_corr
    FROM gaiadr3.gaia_source
    WHERE 1 = CONTAINS(
        POINT('ICRS', ra, dec),
        CIRCLE('ICRS', {ra_cen}, {dec_cen}, {radius_deg}))
    """
    Gaia.ROW_LIMIT = row_limit
    import time
    for attempt in range(5):
        try:
            job = Gaia.launch_job_async(query)
            return job.get_results().to_pandas()
        except Exception as e:
            print(f"TAP error (attempt {attempt+1}/5): {e}")
            time.sleep(2 ** attempt)
    raise RuntimeError("Failed to query Gaia TAP after 5 attempts.")


def apply_zero_point_correction(df):
    """Add 'parallax_corrected': parallax minus the Lindegren+21 zero-point.
    Source-dependent correction where the package + columns + valid solution
    type (31 or 95) are available; documented global fallback otherwise --
    never a hardcoded solution type."""
    df['parallax_corrected'] = df['parallax'] - GLOBAL_ZPT_FALLBACK_MAS

    needed = {'nu_eff_used_in_astrometry', 'pseudocolour', 'ecl_lat',
              'astrometric_params_solved', 'phot_g_mean_mag'}
    if not (ZPT_AVAILABLE and needed.issubset(df.columns)):
        warnings.warn(
            "gaiadr3-zeropoint package or required columns not available -- "
            f"using the fixed global offset of {GLOBAL_ZPT_FALLBACK_MAS * 1000:.0f} uas "
            "for all sources. Install with pip install gaiadr3-zeropoint for the "
            "source-dependent correction."
        )
        return df

    try:
        zpt.load_tables()
        valid_31 = (df['astrometric_params_solved'] == 31) & df['nu_eff_used_in_astrometry'].notna()
        valid_95 = (df['astrometric_params_solved'] == 95) & df['pseudocolour'].notna()
        valid = (valid_31 | valid_95) & df['ecl_lat'].notna() & df['phot_g_mean_mag'].notna()
        
        if valid.any():
            zp = zpt.get_zpt(
                np.asarray(df.loc[valid, 'phot_g_mean_mag'].values, dtype=np.float64),
                np.asarray(df.loc[valid, 'nu_eff_used_in_astrometry'].fillna(0).values, dtype=np.float64),
                np.asarray(df.loc[valid, 'pseudocolour'].fillna(0).values, dtype=np.float64),
                np.asarray(df.loc[valid, 'ecl_lat'].values, dtype=np.float64),
                np.asarray(pd.to_numeric(df.loc[valid, 'astrometric_params_solved']).fillna(0), dtype=np.int64),
            )
            zp = np.asarray(zp, dtype=float)
            df.loc[valid, 'parallax_corrected'] = df.loc[valid, 'parallax'] - zp
            print(f"Applied source-dependent Lindegren+21 zero-point correction to "
                  f"{valid.sum()}/{len(df)} sources (median = {np.nanmedian(zp) * 1000:.1f} uas); "
                  f"remaining sources use the global {GLOBAL_ZPT_FALLBACK_MAS * 1000:.0f} uas fallback.")
        else:
            warnings.warn("No sources had valid inputs for the source-dependent zero-point "
                           "correction -- using the global fallback offset for all sources.")
    except Exception as e:
        warnings.warn(f"gaiadr3-zeropoint correction failed ({e}) -- using the global "
                       f"{GLOBAL_ZPT_FALLBACK_MAS * 1000:.0f} uas fallback offset instead.")
    return df



# Main #

def main():
    parser = argparse.ArgumentParser(
        description="GaussCluster: Iterative 3D Bayesian Mixture Model "
                    "(full 3x3 covariance, zero-point corrected, KDE field model, Mahalanobis Truncated)")
    parser.add_argument('--seed', required=True)
    parser.add_argument('--outdir', required=True)
    parser.add_argument('--ra-col', default='RA')
    parser.add_argument('--dec-col', default='Dec')
    parser.add_argument('--name-col', default=None)
    parser.add_argument('--cluster-name', default=None)
    parser.add_argument('--radius', type=float, default=0.5)
    parser.add_argument('--threshold', type=float, default=0.9)
    parser.add_argument('--max-iter', type=int, default=10)
    parser.add_argument('--tol', type=float, default=0.005, help="Convergence tol on mean+sigma shift.")
    parser.add_argument('--max-pm-error', type=float, default=0.5)
    parser.add_argument('--max-ruwe', type=float, default=1.4)
    parser.add_argument('--excl-sigma', type=float, default=4.0,
                         help="Core/field split radius, in normalised sigma units.")
    parser.add_argument('--kde-cap', type=int, default=50000)
    parser.add_argument('--spatial-cut', type=float, default=15.0, help="Maximum distance (in arcmin) from the median RA/Dec to keep a member.")
    parser.add_argument('--storm', default=None, help="STORM catalogue for comparison.")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    print("=== GaussCluster Initialization ===")

    # 1. Load seed data
    try:
        seed_df = pd.read_csv(args.seed, low_memory=False)
        if args.name_col and args.cluster_name:
            seed_df = seed_df[seed_df[args.name_col] == args.cluster_name]
        if args.ra_col not in seed_df.columns and 'RA_ICRS' in seed_df.columns:
            args.ra_col = 'RA_ICRS'
        if args.dec_col not in seed_df.columns and 'DE_ICRS' in seed_df.columns:
            args.dec_col = 'DE_ICRS'
        if len(seed_df) < 3:
            print("Error: Too few seed members (< 3) to initialize cluster covariance. Exiting.")
            sys.exit(1)
        ra_cen = seed_df[args.ra_col].median()
        dec_cen = seed_df[args.dec_col].median()
        print(f"Loaded {len(seed_df)} seed members. Cluster Center: RA={ra_cen:.4f}, Dec={dec_cen:.4f}")
    except Exception as e:
        print(f"Error loading seed file: {e}")
        sys.exit(1)

    # 2. Query Gaia DR3 (official TAP archive)
    print(f"\nQuerying Gaia DR3 (TAP+) within {args.radius} degrees...")
    try:
        df = query_gaia_dr3(ra_cen, dec_cen, args.radius)
    except Exception as e:
        print(f"Error querying Gaia archive: {e}")
        sys.exit(1)
    if len(df) == 0:
        print("Error: No sources returned from Gaia archive.")
        sys.exit(1)

    quality_mask = (
        df['pmra'].notna() & df['pmdec'].notna() & df['parallax'].notna() &
        (df['pmra_error'] <= args.max_pm_error) &
        (df['pmdec_error'] <= args.max_pm_error) &
        (df['ruwe'] < args.max_ruwe)
    )
    df = df[quality_mask].copy().reset_index(drop=True)
    print(f"Retrieved {len(df)} astrometrically clean sources.")
    if len(df) < 200:
        print("Warning: very few sources survive quality cuts -- consider a larger --radius "
              "or relaxing --max-pm-error / --max-ruwe.")

    bp_rp2 = df['bp_rp'] ** 2
    c_upper = 1.3 + 0.06 * bp_rp2
    c_lower = 1.0 + 0.015 * bp_rp2
    df['phot_quality_ok'] = (
        df['phot_bp_rp_excess_factor'].isna()
        | df['phot_bp_rp_excess_factor'].between(c_lower, c_upper)
    )

    # 2b. Parallax zero-point correction
    df = apply_zero_point_correction(df)

    # 3. Match seeds to Gaia, get an initial kinematic estimate
    seed_coords = SkyCoord(ra=seed_df[args.ra_col].values * u.deg, dec=seed_df[args.dec_col].values * u.deg)
    gaia_coords = SkyCoord(ra=df['ra'].values * u.deg, dec=df['dec'].values * u.deg)
    idx, d2d, _ = seed_coords.match_to_catalog_sky(gaia_coords)
    good = d2d < 1.0 * u.arcsec
    if good.sum() < 3:
        print(f"Error: only {good.sum()} seed members matched to a Gaia source within 1 arcsec -- "
              "check --ra-col/--dec-col and that the seed catalogue is at the Gaia DR3 (J2016.0) epoch.")
        sys.exit(1)
    matched_gaia = df.iloc[idx[good]].copy()

    cluster_mean, cluster_sigma = fit_kinematics(
        matched_gaia['pmra'].values, matched_gaia['pmdec'].values, matched_gaia['parallax_corrected'].values)

    print("\nInitial cluster kinematics (zero-point corrected, MAD-clipped seeds):")
    print(f"pmRA = {cluster_mean[0]:.3f} +- {cluster_sigma[0]:.3f} mas/yr")
    print(f"pmDE = {cluster_mean[1]:.3f} +- {cluster_sigma[1]:.3f} mas/yr")
    print(f"Plx  = {cluster_mean[2]:.3f} +- {cluster_sigma[2]:.3f} mas")

    # 4. Per-source 3x3 error covariance (correlations included)
    N = len(df)
    pmra_e = df['pmra_error'].values
    pmdec_e = df['pmdec_error'].values
    plx_e = df['parallax_error'].values
    c_pmra_pmdec = df['pmra_pmdec_corr'].fillna(0).values
    c_plx_pmra = df['parallax_pmra_corr'].fillna(0).values
    c_plx_pmdec = df['parallax_pmdec_corr'].fillna(0).values

    err_covs = np.zeros((N, 3, 3))
    err_covs[:, 0, 0] = pmra_e ** 2
    err_covs[:, 1, 1] = pmdec_e ** 2
    err_covs[:, 2, 2] = plx_e ** 2
    err_covs[:, 0, 1] = err_covs[:, 1, 0] = pmra_e * pmdec_e * c_pmra_pmdec
    err_covs[:, 0, 2] = err_covs[:, 2, 0] = pmra_e * plx_e * c_plx_pmra
    err_covs[:, 1, 2] = err_covs[:, 2, 1] = pmdec_e * plx_e * c_plx_pmdec
    err_var = np.stack([pmra_e ** 2, pmdec_e ** 2, plx_e ** 2], axis=1)  # for dispersion deconvolution

    all_data = np.vstack([df['pmra'].values, df['pmdec'].values, df['parallax_corrected'].values]).T

    # 5. Static background KDE (fit once at the initial estimate
    r3_init = np.sqrt(np.sum(((all_data - cluster_mean) / cluster_sigma) ** 2, axis=1))
    kde_field_mask = r3_init > args.excl_sigma
    kde_field_df = df[kde_field_mask]
    if len(kde_field_df) < 100:
        print("Error: too few field sources to build a background model -- try a larger --radius.")
        sys.exit(1)
    kde_sample = kde_field_df.sample(args.kde_cap, random_state=42) if len(kde_field_df) > args.kde_cap else kde_field_df
    print(f"\nComputing background KDE from {len(kde_sample)} field sources "
          f"(this may take 1-2 minutes)...")
    field_data = np.vstack([kde_sample['pmra'].values, kde_sample['pmdec'].values,
                             kde_sample['parallax_corrected'].values])
    kde = gaussian_kde(field_data, bw_method='scott')
    phi_f = kde.evaluate(all_data.T) + 1e-10

    # 6. Vectorized E-M: mean AND covariance refit every iteration
    print(f"\nStarting E-M iterative refinement (max {args.max_iter} iterations)...")
    floors = np.array([0.05, 0.05, 0.02])
    P_new = None
    for it in range(1, args.max_iter + 1):
        r3 = np.sqrt(np.sum(((all_data - cluster_mean) / cluster_sigma) ** 2, axis=1))
        core_mask = r3 <= args.excl_sigma
        n_in_core = int(core_mask.sum())
        n_field = N - n_in_core
        max_r = np.max(r3)
        if max_r <= args.excl_sigma:
            print("Error: no sources found outside the core exclusion radius -- check --radius.")
            sys.exit(1)
        background_density = n_field / (max_r ** 3 - args.excl_sigma ** 3)
        expected_background = background_density * args.excl_sigma ** 3
        true_cluster_N = max(len(seed_df), n_in_core - expected_background)
        nc = min(max(true_cluster_N / N, 1e-6), 0.999999)
        nf = 1.0 - nc

        cluster_cov = np.diag(cluster_sigma ** 2)
        total_covs = cluster_cov[np.newaxis, :, :] + err_covs
        inv_covs = np.linalg.inv(total_covs)
        det_covs = np.linalg.det(total_covs)
        norm_consts = 1.0 / np.sqrt((2 * np.pi) ** 3 * det_covs)
        diff = all_data - cluster_mean
        
        # Calculate full Mahalanobis Distance Squared for each source
        D_M_sq = np.einsum('ni,nij,nj->n', diff, inv_covs, diff)
        D_M = np.sqrt(D_M_sq)
        
        exponent = -0.5 * D_M_sq
        phi_c = norm_consts * np.exp(exponent)
        
        # MAHALANOBIS TRUNCATION: Any star further than 3 sigma from the 
        # cluster mean in covariance space has its cluster probability forced to exactly 0. 
        phi_c[D_M > 3.0] = 0.0

        P_new = (nc * phi_c) / (nc * phi_c + nf * phi_f)
        df['P_memb_3D'] = P_new

        weights = P_new
        sum_w = max(weights.sum(), 1.0)
        new_mean = np.average(all_data, weights=weights, axis=0)

        # Error-deconvolved dispersion: raw weighted scatter is intrinsic +
        # measurement error, so subtract the weighted mean error variance
        # before flooring, to avoid inflating sigma iteration-over-iteration.
        var_obs = np.average((all_data - new_mean) ** 2, weights=weights, axis=0)
        mean_err_var = np.average(err_var, weights=weights, axis=0)
        new_sigma = np.sqrt(np.maximum(var_obs - mean_err_var, floors ** 2))

        mean_shift = np.max(np.abs(new_mean - cluster_mean))
        sigma_shift = np.max(np.abs(new_sigma - cluster_sigma))
        n_members_now = int((P_new >= args.threshold).sum())
        print(f"Iter {it}: mean_shift={mean_shift:.5f}, sigma_shift={sigma_shift:.5f}, "
              f"N_expected={true_cluster_N:.1f}, N_members(P>={args.threshold})={n_members_now}, "
              f"plx_c={new_mean[2]:.3f}")

        cluster_mean, cluster_sigma = new_mean, new_sigma
        if mean_shift < args.tol and sigma_shift < args.tol:
            print(f"Converged after {it} iterations.")
            break

    members = df[df['P_memb_3D'] >= args.threshold].copy()
    
    if args.spatial_cut is not None:
        med_ra = members['ra'].median()
        med_dec = members['dec'].median()
        center = SkyCoord(ra=med_ra*u.deg, dec=med_dec*u.deg)
        members_coords = SkyCoord(ra=members['ra'].values*u.deg, dec=members['dec'].values*u.deg)
        seps = center.separation(members_coords).to(u.arcmin).value
        valid = seps <= args.spatial_cut
        members = members[valid]
        print(f"Applied spatial cut of {args.spatial_cut} arcmin. Remaining members: {len(members)}")

    out_csv = os.path.join(args.outdir, "GaussCluster_Members.csv")
    members.to_csv(out_csv, index=False)
    df.to_csv(os.path.join(args.outdir, 'GaussCluster_Field.csv'), index=False)
    print(f"\nSUCCESS: Extracted {len(members)} true 3D members (P >= {args.threshold}).")
    print(f"Catalog saved to: {out_csv}")
    print(f"Final: pmRA={cluster_mean[0]:.3f}+-{cluster_sigma[0]:.3f}, "
          f"pmDE={cluster_mean[1]:.3f}+-{cluster_sigma[1]:.3f}, "
          f"Plx={cluster_mean[2]:.3f}+-{cluster_sigma[2]:.3f} mas"
          + (f" (-> {1000.0 / cluster_mean[2]:.1f} pc)" if cluster_mean[2] > 0 else ""))

    # 7. Diagnostic plot
    fig, axes = plt.subplots(1, 4, figsize=(24, 6))

    axes[0].scatter(df['pmra'], df['pmdec'], s=1, color='grey', alpha=0.1, label='Field')
    axes[0].scatter(members['pmra'], members['pmdec'], s=10, marker='+', color='blue',
                     label=f'GaussCluster ({len(members)})')

    storm_matched = None
    if args.storm is not None:
        try:
            storm_df = pd.read_csv(args.storm, low_memory=False)
            s_ra = 'RA_ICRS' if 'RA_ICRS' in storm_df.columns else 'RA' if 'RA' in storm_df.columns else 'ra'
            s_dec = 'DE_ICRS' if 'DE_ICRS' in storm_df.columns else 'Dec' if 'Dec' in storm_df.columns else 'dec'
            sc_storm = SkyCoord(ra=storm_df[s_ra].values * u.deg, dec=storm_df[s_dec].values * u.deg)
            sidx, sd2d, _ = sc_storm.match_to_catalog_sky(gaia_coords)
            storm_matched = df.iloc[sidx[sd2d < 1.0 * u.arcsec]].copy()
            axes[0].scatter(storm_matched['pmra'], storm_matched['pmdec'], facecolors='none',
                             edgecolors='red', s=40, label=f'STORM ({len(storm_matched)})')
        except Exception as e:
            print(f"Failed to overlay STORM: {e}")

    axes[0].set_xlim(cluster_mean[0] - 5, cluster_mean[0] + 5)
    axes[0].set_ylim(cluster_mean[1] - 5, cluster_mean[1] + 5)
    axes[0].set_xlabel('pmRA (mas/yr)', fontsize=15)
    axes[0].set_ylabel('pmDE (mas/yr)', fontsize=15)
    axes[0].tick_params(axis='both', which='major', labelsize=12)
    axes[0].legend(fontsize=12, loc='best')
    axes[0].set_title('Vector Point Diagram', fontsize=20)

    axes[1].scatter(df['ra'], df['dec'], s=1, color='grey', alpha=0.05)
    axes[1].scatter(members['ra'], members['dec'], s=10, marker='+', color='blue', label='GaussCluster')
    if storm_matched is not None:
        axes[1].scatter(storm_matched['ra'], storm_matched['dec'], facecolors='none', edgecolors='red',
                         s=40, label='STORM')
    axes[1].set_xlabel('RA (deg)', fontsize=15)
    axes[1].set_ylabel('Dec (deg)', fontsize=15)
    axes[1].tick_params(axis='both', which='major', labelsize=12)
    axes[1].legend(fontsize=12, loc='best')
    axes[1].set_title('Spatial Distribution', fontsize=20)

    axes[2].hist(df['parallax_corrected'], bins=50, range=(cluster_mean[2] - 3, cluster_mean[2] + 3), alpha=0.3,
                 color='grey', label='Field', density=True)
    axes[2].hist(members['parallax_corrected'], bins=20, alpha=0.5, color='blue',
                 label='GaussCluster', density=True)
    if storm_matched is not None:
        axes[2].hist(storm_matched['parallax_corrected'], bins=15, alpha=0.5, color='red',
                     label='STORM', density=True, histtype='step', linewidth=2)
    axes[2].set_xlim(cluster_mean[2] - 1.5, cluster_mean[2] + 1.5)
    axes[2].set_xlabel('Parallax, zero-point corrected (mas)', fontsize=15)
    axes[2].set_ylabel('Density', fontsize=15)
    axes[2].tick_params(axis='both', which='major', labelsize=12)
    axes[2].legend(fontsize=12, loc='best')
    axes[2].set_title('Parallax Distribution', fontsize=20)

    # CMD panel uses the phot_quality_ok flag (both C* bounds) -- for
    # display/isochrone-fitting purposes only, never for membership itself.
    cmd_members = members[members['phot_quality_ok']]
    cmd_storm = storm_matched[storm_matched['phot_quality_ok']] if storm_matched is not None else None
    axes[3].scatter(cmd_members['bp_rp'], cmd_members['phot_g_mean_mag'], s=10, marker='+', color='blue',
                     label='GaussCluster')
    if cmd_storm is not None:
        axes[3].scatter(cmd_storm['bp_rp'], cmd_storm['phot_g_mean_mag'], facecolors='none',
                         edgecolors='red', s=40, label='STORM')
    axes[3].invert_yaxis()
    axes[3].set_xlabel('BP-RP (mag)', fontsize=15)
    axes[3].set_ylabel('G (mag)', fontsize=15)
    axes[3].tick_params(axis='both', which='major', labelsize=12)
    axes[3].legend(fontsize=12, loc='best')
    axes[3].set_title('Colour-Magnitude Diagram', fontsize=20)

    plt.tight_layout()
    out_img = os.path.join(args.outdir, 'GaussCluster_Diagnostics_Cosmetic.png')
    plt.savefig(out_img, dpi=150)
    print(f"Diagnostic plot saved to: {out_img}")
    print("=== Run Complete ===")


if __name__ == "__main__":
    main()
