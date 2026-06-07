import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import math
from glob import glob

# ---------------------------
# PARAMETERS
# ---------------------------
CSV_FILE = "../../Data/Raw/MSD_experiment.csv"
PX_PER_CM = 20.0               # pixels per cm (for unit conversion if needed)
CENTER = (1601 / PX_PER_CM, 1588 / PX_PER_CM)             # circle center (x0, y0)
RADIUS = 370 / PX_PER_CM               # circle radius
DT = 1/6.0                       # time per frame (units per frame).
MIN_SEGMENT_LENGTH = 120          # ignore segments shorter than this (too noisy)
MAX_LAG = 120                  # maximum lag (in frames) to compute MSD for a segment
FIT_LAGS = 6                    # number of small-lag points to use when fitting quadratic
# ---------------------------

# Quadratic model for MSD: a*tau^2 + b*tau + c
def quad_model(tau, a, b, c):
    return a * tau ** 2 + b * tau + c

def compute_msd_segment(x, y, max_lag=None):
    """
    Compute MSD for a single segment with coordinates x, y (1D arrays).
    Returns taus (array of lag times in seconds) and msd (array).
    """
    n = len(x)
    if n < 2:
        return np.array([]), np.array([])

    if max_lag is None:
        max_lag = n - 1
    else:
        max_lag = min(max_lag, n - 1)

    taus = np.arange(1, max_lag + 1) * DT
    msd = np.empty(len(taus), dtype=float)
    for idx, lag in enumerate(range(1, max_lag + 1)):
        dx = x[lag:] - x[:-lag]
        dy = y[lag:] - y[:-lag]
        sqdisp = dx * dx + dy * dy
        msd[idx] = np.mean(sqdisp) if sqdisp.size > 0 else np.nan
    return taus, msd

def split_into_inside_out_segments(df_particle, center, radius):
    """
    Given a DataFrame for one particle (with columns frame,x,y),
    return a list of segments: each is dict {inside: bool, frames, x, y}.
    Only contiguous runs with same inside/outside status are returned.
    """
    x = df_particle['x'].values  / PX_PER_CM
    y = df_particle['y'].values  / PX_PER_CM
    frames = df_particle['frame'].values

    # Smooth values
    x = np.convolve(x, np.ones(7)/7, mode='valid')
    y = np.convolve(y, np.ones(7)/7, mode='valid')

    # compute inside mask
    dx = x - center[0]
    dy = y - center[1]
    dist2 = dx * dx + dy * dy
    inside_mask = dist2 <= radius * radius  # boolean array

    if len(inside_mask) == 0:
        return []

    # find contiguous runs
    runs = []
    start = 0
    current = inside_mask[0]
    for i in range(1, len(inside_mask)):
        if inside_mask[i] != current:
            # close run [start, i)
            runs.append({
                'inside': bool(current),
                'frames': frames[start:i],
                'x': x[start:i],
                'y': y[start:i]
            })
            start = i
            current = inside_mask[i]
    # final run
    runs.append({
        'inside': bool(current),
        'frames': frames[start:len(inside_mask)],
        'x': x[start:len(inside_mask)],
        'y': y[start:len(inside_mask)]
    })
    return runs

def fit_velocity_from_msd(taus, msd, fit_points=6):
    """
    Fit a quadratic MSD(τ) = a τ^2 + b τ + c over the first fit_points.
    Returns estimated velocity v = sqrt(a) (or np.nan if fit fails).
    taus in seconds, msd in distance^2 units.
    """
    if len(taus) < 2 or len(msd) < 2:
        return np.nan, None
    n_fit = min(fit_points, len(taus))
    x = taus[:n_fit]
    y = msd[:n_fit]
    # Guard: need finite y
    mask = np.isfinite(y)
    if np.sum(mask) < 2:
        return np.nan, None
    try:
        popt, pcov = curve_fit(quad_model, x[mask], y[mask], p0=[0.0, 0.0, y[0] if len(y)>0 else 0.0])
        a = popt[0]
        if a < 0:
            # Negative a -> can't take sqrt; return NaN
            return np.nan, popt
        v = math.sqrt(a)
        return v, popt
    except Exception as e:
        return np.nan, None

def main():
    print(f"Processing file: {CSV_FILE}")
    # Load CSV
    df = pd.read_csv(CSV_FILE)
    required_cols = {'frame', 'x', 'y', 'particle'}
    if not required_cols.issubset(df.columns):
        raise ValueError(f"CSV must contain columns: {required_cols}")

    # Sort for safety
    df = df.sort_values(['particle', 'frame']).reset_index(drop=True)

    # storage
    msd_segments_inside = []   # list of (taus, msd) for inside segments
    msd_segments_outside = []  # same for outside segments
    msd_segments_all = []  # largest segments

    # iterate particles
    for pid, group in df.groupby('particle'):
        # ensure sorted by frame
        group = group.sort_values('frame')
        segments = split_into_inside_out_segments(group, CENTER, RADIUS)
        for seg in segments:
            seg_len = len(seg['x'])
            if seg_len < MIN_SEGMENT_LENGTH:
                continue
            taus, msd = compute_msd_segment(np.array(seg['x']), np.array(seg['y']), max_lag=MAX_LAG)
            if len(taus) == 0:
                continue
            if seg['inside']:
                msd_segments_inside.append((taus, msd))
            else:
                msd_segments_outside.append((taus, msd))
        
        large_segments = split_into_inside_out_segments(group, CENTER, 1400)  # large radius to get full segments
        for seg in large_segments:
            seg_len = len(seg['x'])
            if seg_len < MIN_SEGMENT_LENGTH:
                continue
            taus, msd = compute_msd_segment(np.array(seg['x']), np.array(seg['y']), max_lag=MAX_LAG)
            if len(taus) == 0:
                continue
            msd_segments_all.append((taus, msd))

    # If no segments found, exit
    if len(msd_segments_inside) + len(msd_segments_outside) == 0:
        raise RuntimeError("No segments found (maybe MIN_SEGMENT_LENGTH too large?). Exiting.")

    # Plotting: first plot all individual msds in black translucent
    plt.figure(figsize=(6,8))
    ax = plt.gca()

    # Determine a common tau axis for computing mean: choose a reasonable max common lag
    # Find minimal max lag among segments to avoid NaNs when averaging.
    def common_tau_axis(segments):
        if len(segments) == 0:
            return np.array([])
        min_maxlag = min([seg[0][-1] if len(seg[0])>0 else 0 for seg in segments])
        if min_maxlag <= 0:
            return np.array([])
        # We'll sample integer lags in frames: taus are multiples of DT; re-create index
        max_index = int(min_maxlag / DT)
        taus_common = (np.arange(1, max_index + 1) * DT)
        return taus_common

    taus_in_common = common_tau_axis(msd_segments_inside)
    taus_out_common = common_tau_axis(msd_segments_outside)
    # to compute means, we'll use the shorter of the two commons (so inside/outside means use same tau axis)
    if len(taus_in_common) == 0 and len(taus_out_common) == 0:
        raise RuntimeError("No valid common taus for averaging.")
    if len(taus_in_common) == 0:
        taus_common = taus_out_common
    elif len(taus_out_common) == 0:
        taus_common = taus_in_common
    else:
        taus_common = taus_in_common if len(taus_in_common) <= len(taus_out_common) else taus_out_common

    # Helper to interpolate msd onto taus_common
    def interp_msd_to_common(taus, msd, taus_common):
        if len(taus_common) == 0:
            return np.array([])
        # we assume taus is monotonically increasing
        # only consider overlap range
        from numpy import interp
        # interp requires x increasing; ensure float arrays
        return interp(taus_common, taus, msd, left=np.nan, right=np.nan)

    # Plot individual segments
    all_inside_msds_on_common = []
    for taus, msd in msd_segments_inside:
        ax.plot(taus, msd, color="#400000", alpha=0.02, linewidth=0.8)
        if len(taus_common) > 0:
            all_inside_msds_on_common.append(interp_msd_to_common(taus, msd, taus_common))

    all_outside_msds_on_common = []
    for taus, msd in msd_segments_outside:
        ax.plot(taus, msd, color="#000040", alpha=0.02, linewidth=0.8)
        if len(taus_common) > 0:
            all_outside_msds_on_common.append(interp_msd_to_common(taus, msd, taus_common))
    

    all_large_msds_on_common = []
    for taus, msd in msd_segments_all:
        if len(taus_common) > 0:
            all_large_msds_on_common.append(interp_msd_to_common(taus, msd, taus_common))

    # Convert to arrays and compute mean across segments (ignoring NaNs)
    mean_inside = None
    mean_outside = None
    v_inside = np.nan
    v_outside = np.nan

    if len(all_inside_msds_on_common) > 0:
        arr = np.array(all_inside_msds_on_common)  # shape (n_segments, n_taus)
        mean_inside = np.nanmean(arr, axis=0)
        # fit velocity on mean
        v_inside, popt_in = fit_velocity_from_msd(taus_common, mean_inside, fit_points=FIT_LAGS)
    if len(all_outside_msds_on_common) > 0:
        arr = np.array(all_outside_msds_on_common)
        mean_outside = np.nanmean(arr, axis=0)
        v_outside, popt_out = fit_velocity_from_msd(taus_common, mean_outside, fit_points=FIT_LAGS)
    
    estimated_velocity, _ = fit_velocity_from_msd(taus_common, np.nanmean(np.array(all_large_msds_on_common), axis=0), fit_points=FIT_LAGS)

    # Plot mean curves
    if mean_inside is not None:
        ax.plot(taus_common, mean_inside, color='blue', linewidth=2.2,
                label="$v_{\\circ}\\approx$ " + f"{v_inside:.2f} cm/s")
        # Plot line showing fit and a rotated annotation for velocity
        if not np.isnan(v_inside):
            fit_line = quad_model(taus_common, *(popt_in)) if popt_in is not None else None
            if fit_line is not None:
                ax.plot(taus_common[1:6], fit_line[1:6]/2, color='darkblue', linestyle='--', linewidth=1.5)
                # Annotate velocity
                tau_annot = taus_common[min(FIT_LAGS-1, len(taus_common)-1)]
                msd_annot = quad_model(tau_annot, *(popt_in))
                angle = math.degrees(math.atan2(1, 0.5))
                ax.annotate(f"$\propto \\tau^2$",
                            xy=(0.4, 0.04),
                            xytext=(10, 10),
                            textcoords='offset points',
                            rotation=angle,
                            color='darkblue',
                            fontsize=20)
        # Plot line showing fit and a rotated annotation for velocity
        if not np.isnan(v_inside):
            fit_line = quad_model(taus_common, *(popt_out)) if popt_out is not None else None
            if fit_line is not None:
                ax.plot(taus_common[1:6], fit_line[1:6]*2, color='darkred', linestyle='--', linewidth=1.5)
                # Annotate velocity
                tau_annot = taus_common[min(FIT_LAGS-1, len(taus_common)-1)]
                msd_annot = quad_model(tau_annot, *(popt_out))
                angle = math.degrees(math.atan2(1, 0.5))
                ax.annotate(f"$\propto \\tau^2$",
                            xy=(0.17, 2.0),
                            xytext=(10, 10),
                            textcoords='offset points',
                            rotation=angle,
                            color='darkred',
                            fontsize=20)
    if mean_outside is not None:
        ax.plot(taus_common, mean_outside, color='red', linewidth=2.2,
                label="$v_{\\bullet}\\approx$ " + f"{v_outside:.2f} cm/s")
    ax.set_xscale('log' if len(taus_common)>0 and taus_common[-1]/taus_common[0]>50 else 'linear')
    ax.set_yscale('log' if (mean_inside is not None and np.nanmax(mean_inside)/np.nanmin(mean_inside[np.nonzero(mean_inside)])>50) or
                (mean_outside is not None and np.nanmax(mean_outside)/np.nanmin(mean_outside[np.nonzero(mean_outside)])>50)
                else 'linear')
    plt.ylim(10e-3, 2*10e2)
    ax.grid(True, which='both', linestyle='--', linewidth=0.5, alpha=0.6)
    plt.gca().set_aspect('equal')
    plt.xlabel('$\\tau$ (s)')
    plt.ylabel('MSD ($cm^2$)')  
    plt.tight_layout()
    OUTPUT_FIG = f"../out/MSD_plot.png"

    handles, labels = ax.get_legend_handles_labels()
    # sort both labels and handles by labels
    labels, handles = zip(*sorted(zip(labels, handles), key=lambda t: t[0]))
    ax.legend(handles, labels, fontsize=19, loc='upper left')
    plt.savefig(OUTPUT_FIG, dpi=200)
    print(f"Saved MSD plot to: {OUTPUT_FIG}")
    print(f"Segments used: inside={len(msd_segments_inside)}, outside={len(msd_segments_outside)}")
    print(f"Estimated mean velocity (from largest segments): {estimated_velocity:.6g} cm / s")
    if not np.isnan(v_inside):
        print(f"Estimated mean velocity inside (from fit): {v_inside:.6g} cm / s")
    else:
        print("Estimated mean velocity inside: nan (fit failed or not enough data)")
    if not np.isnan(v_outside):
        print(f"Estimated mean velocity outside (from fit): {v_outside:.6g} cm / s")
    else:
        print("Estimated mean velocity outside: nan (fit failed or not enough data)")
    plt.show()

if __name__ == "__main__":
    plt.rcParams.update({
        "text.usetex": True,
        "font.family": "Times New Roman",
        "font.size": 24
    })
    main()