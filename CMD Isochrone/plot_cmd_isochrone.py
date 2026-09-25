import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import seaborn as sns

def plot_cluster(cluster_name, members_path, iso_path, out_dir, e_bv_adopted, best_logage, best_dm, color, xlim, label_val=3):
    print(f"Plotting {cluster_name}...")
    
    A_G_COEFF = 2.27 * e_bv_adopted
    A_BP_COEFF = 2.98 * e_bv_adopted
    A_RP_COEFF = 1.72 * e_bv_adopted
    E_BPRP = A_BP_COEFF - A_RP_COEFF
    
    df = pd.read_csv(members_path)
    if 'bp_rp' not in df.columns:
        df['bp_rp'] = df['phot_bp_mean_mag'] - df['phot_rp_mean_mag']

    with open(iso_path, "r") as f:
        for line in f:
            stripped = line.lstrip("#").strip()
            if "logAge" in stripped and "Gmag" in stripped:
                col_names = stripped.split()
                break
                
    iso = pd.read_csv(iso_path, comment="#", sep=r"\s+", header=None, names=col_names)
    
    if "G_BPmag" in iso.columns:
        iso.rename(columns={"G_BPmag": "BPmag", "G_RPmag": "RPmag"}, inplace=True)
    iso["BP_RP"] = iso["BPmag"] - iso["RPmag"]
    
    ages = np.sort(iso["logAge"].unique())
    closest_age = ages[np.argmin(np.abs(ages - best_logage))]
    best_track = iso[iso["logAge"] == closest_age].copy()
    
    print(f"Fitting Isochrone: log(Age) = {closest_age} ({10**closest_age / 1e9:.2f} Gyr / {10**closest_age / 1e6:.2f} Myr)")
    print(f"Distance Modulus: {best_dm} ({10**(best_dm/5 + 1):.0f} pc)")
    print(f"Reddening E(B-V): {e_bv_adopted}")
    
    best_track = best_track[best_track["label"] <= label_val].sort_values("Mini")
    brightest_idx = best_track["Gmag"].idxmin()
    max_mini = best_track.loc[brightest_idx, "Mini"]
    best_track = best_track[best_track["Mini"] <= max_mini].sort_values("Mini")
    
    sns.set_theme(style="white", context="paper")
    fig, ax = plt.subplots(figsize=(8, 10))
    
    # Scatter of cluster members
    ax.scatter(df["bp_rp"], df["phot_g_mean_mag"], s=45, color=color,
               edgecolors='none', alpha=0.9, zorder=2, label=f"{cluster_name} Members")
               
    # Plot the Isochrone track
    ax.plot(best_track["BP_RP"] + E_BPRP, best_track["Gmag"] + best_dm + A_G_COEFF,
            color="crimson", linewidth=2.5, zorder=3, alpha=1.0,
            label=f"PARSEC Isochrone\nAge ~ {10**closest_age / 1e6:.2f} Myr, d ~ {10**(best_dm/5 + 1):.0f} pc")
            
    ax.invert_yaxis()
    ax.set_xlim(xlim[0], xlim[1])
    ymin, ymax = ax.get_ylim()
    ax.set_ylim(20.0, ymax)
    ax.set_xlabel("BP - RP (mag)", fontsize=15, fontweight='bold')
    ax.set_ylabel("Apparent G Magnitude (mag)", fontsize=15, fontweight='bold')
    ax.tick_params(axis='both', labelsize=12)
    ax.set_title(f"Color-Magnitude Diagram of {cluster_name}", fontsize=18, fontweight='bold', pad=15)
    
    for spine in ax.spines.values():
        spine.set_edgecolor('black')
        spine.set_linewidth(1.5)
        
    ax.legend(fontsize=12, loc='lower left', frameon=True, edgecolor='black')
    
    plt.tight_layout()
    out_file = os.path.join(out_dir, f"{cluster_name.replace(' ', '')}_CMD_Isochrone_Fit.png")
    plt.savefig(out_file, dpi=300)
    print(f"Saved CMD plot to {out_file}\n")
    plt.close()

if __name__ == '__main__':
    # Be87
    plot_cluster(
        cluster_name="Be87",
        members_path=r"C:\Users\raffe\Desktop\Dissertation\GaussCluster_Software\src\v8\Results\Be87\Be87_GaussCluster_Members.csv",
        iso_path=r"C:\Users\raffe\Desktop\Dissertation\CMD Isochrone\v11\Be87\BE87.txt",
        out_dir=r"C:\Users\raffe\Desktop\Dissertation\CMD Isochrone\v11\Be87",
        e_bv_adopted=1.68,
        best_logage=7.0,
        best_dm=11.25,
        color="black",
        xlim=(1, 3),
        label_val=3
    )
    
    # NGC 957
    plot_cluster(
        cluster_name="NGC 957",
        members_path=r"C:\Users\raffe\Desktop\Dissertation\GaussCluster_Software\src\v8\Results\NGC957\NGC957_GaussCluster_Members.csv",
        iso_path=r"C:\Users\raffe\Desktop\Dissertation\CMD Isochrone\v11\NGC957\NGC957.txt",
        out_dir=r"C:\Users\raffe\Desktop\Dissertation\CMD Isochrone\v11\NGC957",
        e_bv_adopted=0.87,
        best_logage=7.0,
        best_dm=11.71,
        color="darkblue",
        xlim=(0, 3),
        label_val=3
    )









