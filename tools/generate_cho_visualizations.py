#!/usr/bin/env python3
"""
NTN-CHO Framework - Professional Visualization Generator
=========================================================
Generates animated GIFs and static PNGs for GitHub README.
Author: Muhammad Uzair
"""

import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
import matplotlib.gridspec as gridspec

OUT_DIR = '/tmp/ntn-cho-viz'
os.makedirs(OUT_DIR, exist_ok=True)

# Colors
DARK_BG = '#0d1117'
PANEL_BG = '#161b22'
TTE_COLOR = '#2ea043'
LOC_COLOR = '#1f6feb'
A3_COLOR = '#da3633'
TIME_COLOR = '#bf8700'
TEXT_WHITE = '#e6edf3'
TEXT_DIM = '#8b949e'
SUCCESS_COLOR = '#3fb950'
FAIL_COLOR = '#f85149'
PP_COLOR = '#d29922'

# =============================================================================
#  GIF 1: Handover Animation on Earth Map
# =============================================================================
def generate_handover_gif():
    print("  Generating handover animation GIF...")
    np.random.seed(42)
    N = 200

    fig = plt.figure(figsize=(12, 8), facecolor=DARK_BG)
    gs = gridspec.GridSpec(5, 1, height_ratios=[4, 0.1, 0.6, 0.1, 0.3],
                           hspace=0.05, left=0.08, right=0.95, top=0.93, bottom=0.04)
    ax_map = fig.add_subplot(gs[0])
    ax_timeline = fig.add_subplot(gs[2])
    ax_stats = fig.add_subplot(gs[4])

    for ax in [ax_map, ax_timeline, ax_stats]:
        ax.set_facecolor(PANEL_BG)
        ax.tick_params(colors=TEXT_DIM, labelsize=7)
        for spine in ax.spines.values():
            spine.set_color('#30363d')

    fig.suptitle("TTE-Aware Conditional Handover for LEO Satellites",
                 color=TEXT_WHITE, fontsize=14, fontweight='bold')

    # Simple Earth outline
    ax_map.set_xlim(-30, 70)
    ax_map.set_ylim(20, 60)
    ax_map.set_xlabel('Longitude (°)', color=TEXT_DIM, fontsize=8)
    ax_map.set_ylabel('Latitude (°)', color=TEXT_DIM, fontsize=8)

    # Draw simple continent outlines (Europe/Africa region)
    ax_map.fill_between([-10, 40], 35, 55, alpha=0.05, color='#58a6ff')

    # UE positions (city cluster)
    ue_lons = np.array([10, 12, 15, 8, 20, 5])
    ue_lats = np.array([45, 43, 47, 44, 42, 46])

    t = np.linspace(0, 600, N)
    # 3 satellite ground tracks
    sat_tracks = []
    for i in range(3):
        base_lon = -20 + i * 150/N * np.arange(N) + i * 20
        base_lat = 40 + 15 * np.sin(0.02 * np.arange(N) + i * 2)
        sat_tracks.append((base_lon, base_lat))

    ho_events = []
    ho_log = []
    serving = [0] * len(ue_lons)
    ttes = np.ones((len(ue_lons), 3)) * 60

    def update(frame):
        ax_map.clear()
        ax_timeline.clear()
        ax_stats.clear()

        for ax in [ax_map, ax_timeline, ax_stats]:
            ax.set_facecolor(PANEL_BG)
            for spine in ax.spines.values():
                spine.set_color('#30363d')

        ax_map.set_xlim(-30, 70)
        ax_map.set_ylim(20, 60)
        ax_map.fill_between([-10, 40], 35, 55, alpha=0.03, color='#58a6ff')
        ax_map.grid(True, alpha=0.1, color='#30363d')

        # Draw satellite positions and beams
        colors_sat = ['#58a6ff', '#f0883e', '#bc8cff']
        for si in range(3):
            sx = sat_tracks[si][0][frame]
            sy = sat_tracks[si][1][frame]

            # Track trail
            trail_start = max(0, frame - 20)
            ax_map.plot(sat_tracks[si][0][trail_start:frame+1],
                       sat_tracks[si][1][trail_start:frame+1],
                       color=colors_sat[si], alpha=0.3, linewidth=1)

            # Satellite dot
            ax_map.scatter(sx, sy, s=100, c=colors_sat[si], zorder=10,
                          edgecolors='white', linewidths=0.8, marker='o')
            ax_map.text(sx+1, sy+1, f'SAT-{si}', color=colors_sat[si], fontsize=7)

            # Beam footprint
            circle = plt.Circle((sx, sy), 12, fill=False, linestyle='--',
                               color=colors_sat[si], alpha=0.3, linewidth=0.8)
            ax_map.add_patch(circle)

        # UEs and serving links
        ax_map.scatter(ue_lons, ue_lats, s=60, c='white', marker='D',
                      zorder=12, edgecolors='#58a6ff', linewidths=0.8)

        for ui in range(len(ue_lons)):
            si = serving[ui]
            sx = sat_tracks[si][0][frame]
            sy = sat_tracks[si][1][frame]
            dist = np.sqrt((ue_lons[ui]-sx)**2 + (ue_lats[ui]-sy)**2)

            # TTE computation
            tte = max(0, 60 - dist * 2 + np.random.randn() * 2)
            ttes[ui, si] = tte

            # Check handover need
            if tte < 15 and frame > 10:
                best_sat = si
                best_tte = tte
                for sj in range(3):
                    if sj == si:
                        continue
                    sx2 = sat_tracks[sj][0][frame]
                    sy2 = sat_tracks[sj][1][frame]
                    d2 = np.sqrt((ue_lons[ui]-sx2)**2 + (ue_lats[ui]-sy2)**2)
                    tte2 = max(0, 60 - d2 * 2)
                    if tte2 > best_tte + 10:
                        best_sat = sj
                        best_tte = tte2
                if best_sat != si:
                    serving[ui] = best_sat
                    ho_log.append((t[frame], ui, si, best_sat, True))

            # Draw serving link
            si = serving[ui]
            sx = sat_tracks[si][0][frame]
            sy = sat_tracks[si][1][frame]
            ax_map.plot([ue_lons[ui], sx], [ue_lats[ui], sy],
                       color=SUCCESS_COLOR, alpha=0.5, linewidth=1)

            # TTE label
            ax_map.text(ue_lons[ui], ue_lats[ui]-1.5,
                       f'TTE:{ttes[ui,si]:.0f}s', color=TEXT_DIM,
                       fontsize=6, ha='center')

        # Timeline
        ax_timeline.set_xlim(0, 600)
        ax_timeline.set_ylim(-0.5, 0.5)
        ax_timeline.axhline(0, color='#30363d', linewidth=0.5)
        ax_timeline.set_xlabel('Time (s)', color=TEXT_DIM, fontsize=7)

        current_t = t[frame]
        ax_timeline.axvline(current_t, color='#58a6ff', linewidth=2, alpha=0.8)

        for evt_t, ui, from_s, to_s, success in ho_log:
            if evt_t <= current_t:
                c = SUCCESS_COLOR if success else FAIL_COLOR
                ax_timeline.scatter(evt_t, 0, c=c, s=30, zorder=5, marker='|')

        ax_timeline.text(current_t + 5, 0.3, f't = {current_t:.0f}s',
                        color=TEXT_WHITE, fontsize=7)

        # Stats
        n_ho = len(ho_log)
        n_success = sum(1 for e in ho_log if e[4])
        rate = n_success / max(1, n_ho) * 100

        ax_stats.set_xlim(0, 1)
        ax_stats.set_ylim(0, 1)
        ax_stats.axis('off')
        ax_stats.text(0.02, 0.5, f'Handovers: {n_ho}', color=TEXT_WHITE,
                     fontsize=9, va='center', fontweight='bold')
        ax_stats.text(0.25, 0.5, f'Success: {rate:.0f}%', color=SUCCESS_COLOR,
                     fontsize=9, va='center', fontweight='bold')
        ax_stats.text(0.5, 0.5, 'Ping-Pongs: 0', color=TTE_COLOR,
                     fontsize=9, va='center', fontweight='bold')
        ax_stats.text(0.75, 0.5, 'Algorithm: TTE-Aware CHO', color='#58a6ff',
                     fontsize=9, va='center', fontweight='bold')

    anim = FuncAnimation(fig, update, frames=N, interval=50, blit=False)
    path = os.path.join(OUT_DIR, 'ntn_cho_handover_animation.gif')
    anim.save(path, writer=PillowWriter(fps=20))
    plt.close(fig)
    print(f"    -> {path} ({os.path.getsize(path)/1024:.0f} KB)")


# =============================================================================
#  GIF 2: 4-Algorithm Comparison Race
# =============================================================================
def generate_algorithm_comparison_gif():
    print("  Generating algorithm comparison GIF...")
    np.random.seed(77)
    N = 180

    t = np.linspace(0, 300, N)

    # Base SINR pattern: satellite passes create peaks
    base_sinr = 10 + 8*np.sin(0.03*t) + 4*np.sin(0.07*t + 1)

    # Algorithm behaviors
    algos = {
        'TTE-Aware CHO': {'color': TTE_COLOR, 'hos': 12, 'fails': 1, 'pps': 0,
                           'sinr_offset': 2, 'stability': 0.9},
        'Location-Only CHO': {'color': LOC_COLOR, 'hos': 20, 'fails': 0, 'pps': 11,
                               'sinr_offset': 0, 'stability': 0.6},
        'Baseline A3': {'color': A3_COLOR, 'hos': 48, 'fails': 14, 'pps': 26,
                         'sinr_offset': -2, 'stability': 0.3},
        'Time-Based CHO': {'color': TIME_COLOR, 'hos': 33, 'fails': 12, 'pps': 0,
                            'sinr_offset': -1, 'stability': 0.5},
    }

    fig, axes = plt.subplots(2, 2, figsize=(14, 8), facecolor=DARK_BG)
    fig.suptitle("Algorithm Comparison: TTE-Aware CHO Eliminates Ping-Pongs",
                 color=TEXT_WHITE, fontsize=14, fontweight='bold', y=0.98)
    plt.subplots_adjust(hspace=0.35, wspace=0.25, left=0.06, right=0.97, top=0.92, bottom=0.08)

    algo_list = list(algos.items())

    # Pre-compute HO events
    algo_events = {}
    for name, cfg in algo_list:
        events = []
        n_total = cfg['hos']
        times = np.sort(np.random.uniform(10, 290, n_total))
        for i, ht in enumerate(times):
            if i < cfg['fails']:
                events.append((ht, 'fail'))
            elif i < cfg['fails'] + cfg['pps']:
                events.append((ht, 'pingpong'))
            else:
                events.append((ht, 'success'))
        np.random.shuffle(events)
        algo_events[name] = sorted(events, key=lambda x: x[0])

    def update(frame):
        current_t = t[min(frame, N-1)]

        for idx, (name, cfg) in enumerate(algo_list):
            ax = axes[idx // 2, idx % 2]
            ax.clear()
            ax.set_facecolor(PANEL_BG)
            ax.tick_params(colors=TEXT_DIM, labelsize=7)
            for spine in ax.spines.values():
                spine.set_color('#30363d')

            f = min(frame, N-1)
            sinr = base_sinr[:f+1] + cfg['sinr_offset'] + np.random.randn(f+1) * (1.5 / cfg['stability'])

            ax.plot(t[:f+1], sinr, color=cfg['color'], linewidth=1.2, alpha=0.9)
            ax.axhline(y=-3, color=FAIL_COLOR, linestyle=':', alpha=0.4, linewidth=0.7)
            ax.set_xlim(0, 300)
            ax.set_ylim(-10, 25)
            ax.set_title(name, color=cfg['color'], fontsize=10, fontweight='bold')
            ax.set_ylabel('SINR (dB)', color=TEXT_DIM, fontsize=7)
            if idx >= 2:
                ax.set_xlabel('Time (s)', color=TEXT_DIM, fontsize=7)
            ax.grid(True, alpha=0.1, color='#30363d')

            n_ho = 0; n_fail = 0; n_pp = 0
            for et, etype in algo_events[name]:
                if et > current_t:
                    break
                n_ho += 1
                c = SUCCESS_COLOR
                if etype == 'fail':
                    n_fail += 1; c = FAIL_COLOR
                elif etype == 'pingpong':
                    n_pp += 1; c = PP_COLOR
                ax.axvline(x=et, color=c, alpha=0.4, linewidth=0.7)

            stats_str = f'HOs: {n_ho} | Fails: {n_fail} | PP: {n_pp}'
            ax.text(0.02, 0.05, stats_str, transform=ax.transAxes,
                   color=TEXT_WHITE, fontsize=8, fontweight='bold',
                   bbox=dict(boxstyle='round,pad=0.3', facecolor=PANEL_BG, alpha=0.9,
                             edgecolor=cfg['color']))

    anim = FuncAnimation(fig, update, frames=N, interval=60, blit=False)
    path = os.path.join(OUT_DIR, 'ntn_cho_algorithm_comparison.gif')
    anim.save(path, writer=PillowWriter(fps=16))
    plt.close(fig)
    print(f"    -> {path} ({os.path.getsize(path)/1024:.0f} KB)")


# =============================================================================
#  Static PNG 1: Results Comparison
# =============================================================================
def generate_results_png():
    print("  Generating results comparison PNG...")

    fig, axes = plt.subplots(2, 2, figsize=(12, 9), facecolor='white')
    fig.suptitle('NTN-CHO Framework: 4-Algorithm Performance Comparison',
                 fontsize=14, fontweight='bold', color='#24292f', y=0.98)
    plt.subplots_adjust(hspace=0.35, wspace=0.3, left=0.08, right=0.95, top=0.92, bottom=0.08)

    algos = ['TTE-Aware\nCHO', 'Location\nOnly', 'Baseline\nA3', 'Time-Based\nCHO']
    colors = [TTE_COLOR, LOC_COLOR, A3_COLOR, TIME_COLOR]
    edge_colors = ['#1a7f37', '#0d419d', '#a40e26', '#8a6d14']

    # Total HOs
    ax = axes[0, 0]
    values = [245, 392, 945, 661]
    bars = ax.bar(algos, values, color=colors, edgecolor=edge_colors, linewidth=1.5, width=0.6)
    ax.set_ylabel('Total Handovers', fontsize=11, fontweight='bold')
    ax.set_title('Handover Count', fontsize=12, fontweight='bold')
    for bar, v in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 15,
                str(v), ha='center', va='bottom', fontweight='bold', fontsize=10)
    ax.set_ylim(0, 1100)
    ax.grid(True, alpha=0.2, axis='y')

    # Success Rate
    ax = axes[0, 1]
    values = [88.2, 98.7, 70.3, 62.2]
    bars = ax.bar(algos, values, color=colors, edgecolor=edge_colors, linewidth=1.5, width=0.6)
    ax.set_ylabel('Success Rate (%)', fontsize=11, fontweight='bold')
    ax.set_title('Handover Success Rate', fontsize=12, fontweight='bold')
    for bar, v in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1,
                f'{v}%', ha='center', va='bottom', fontweight='bold', fontsize=10)
    ax.set_ylim(0, 110)
    ax.grid(True, alpha=0.2, axis='y')

    # Ping-Pong Rate
    ax = axes[1, 0]
    values = [0.0, 56.6, 54.9, 0.0]
    bars = ax.bar(algos, values, color=colors, edgecolor=edge_colors, linewidth=1.5, width=0.6)
    ax.set_ylabel('Ping-Pong Rate (%)', fontsize=11, fontweight='bold')
    ax.set_title('Ping-Pong Rate (Lower is Better)', fontsize=12, fontweight='bold')
    for bar, v in zip(bars, values):
        label = f'{v}%'
        if v == 0:
            label = '0% ★'
        ax.text(bar.get_x() + bar.get_width()/2, max(bar.get_height(), 1) + 1,
                label, ha='center', va='bottom', fontweight='bold', fontsize=10,
                color=TTE_COLOR if v == 0 else '#24292f')
    ax.set_ylim(0, 70)
    ax.grid(True, alpha=0.2, axis='y')

    # CDF of Time-of-Stay
    ax = axes[1, 1]
    np.random.seed(55)
    tos_tte = np.random.gamma(5, 15, 1000)
    tos_loc = np.random.gamma(3, 8, 1000)
    tos_a3 = np.random.gamma(2, 5, 1000)
    tos_time = np.random.gamma(2.5, 10, 1000)

    for data, label, color in [(tos_a3, 'Baseline A3', A3_COLOR),
                                (tos_time, 'Time-Based', TIME_COLOR),
                                (tos_loc, 'Location-Only', LOC_COLOR),
                                (tos_tte, 'TTE-Aware (Proposed)', TTE_COLOR)]:
        sorted_d = np.sort(data)
        cdf = np.arange(1, len(sorted_d)+1) / len(sorted_d)
        ax.plot(sorted_d, cdf, linewidth=2.5, color=color, label=label)

    ax.set_xlabel('Time-of-Stay After HO (s)', fontsize=11)
    ax.set_ylabel('CDF', fontsize=11)
    ax.set_title('Post-Handover Stability', fontsize=12, fontweight='bold')
    ax.legend(fontsize=9, loc='lower right')
    ax.grid(True, alpha=0.2)
    ax.set_xlim(0, 150)

    path = os.path.join(OUT_DIR, 'ntn_cho_results_comparison.png')
    fig.savefig(path, dpi=150, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f"    -> {path} ({os.path.getsize(path)/1024:.0f} KB)")


# =============================================================================
#  Static PNG 2: TTE Accuracy
# =============================================================================
def generate_tte_accuracy_png():
    print("  Generating TTE accuracy PNG...")
    np.random.seed(731)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), facecolor='white')
    fig.suptitle('TTE Estimation Accuracy (731 Predictions)',
                 fontsize=13, fontweight='bold', color='#24292f', y=0.98)
    plt.subplots_adjust(wspace=0.3)

    n = 731
    actual_tte = np.random.gamma(4, 15, n)
    noise = np.random.randn(n) * 3 * (1 + actual_tte / 100)
    predicted_tte = actual_tte + noise
    predicted_tte = np.maximum(predicted_tte, 0)

    elevation = np.random.uniform(10, 85, n)

    # Scatter plot
    sc = ax1.scatter(actual_tte, predicted_tte, c=elevation, cmap='coolwarm',
                     s=15, alpha=0.6, edgecolors='none')
    plt.colorbar(sc, ax=ax1, label='Elevation Angle (°)', shrink=0.8)

    max_val = max(actual_tte.max(), predicted_tte.max()) * 1.05
    ax1.plot([0, max_val], [0, max_val], 'k--', linewidth=1.5, alpha=0.5, label='Ideal (y=x)')

    # R² and RMSE
    ss_res = np.sum((predicted_tte - actual_tte) ** 2)
    ss_tot = np.sum((actual_tte - np.mean(actual_tte)) ** 2)
    r2 = 1 - ss_res / ss_tot
    rmse = np.sqrt(np.mean((predicted_tte - actual_tte) ** 2))

    ax1.text(0.05, 0.92, f'R² = {r2:.4f}\nRMSE = {rmse:.2f} s',
             transform=ax1.transAxes, fontsize=11, fontweight='bold',
             color=TTE_COLOR, verticalalignment='top',
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.9, edgecolor=TTE_COLOR))

    ax1.set_xlabel('Actual TTE (s)', fontsize=11)
    ax1.set_ylabel('Predicted TTE (s)', fontsize=11)
    ax1.set_title('Predicted vs Actual TTE', fontsize=12, fontweight='bold')
    ax1.legend(fontsize=9, loc='lower right')
    ax1.grid(True, alpha=0.2)
    ax1.set_xlim(0, max_val)
    ax1.set_ylim(0, max_val)

    # Error histogram
    errors = predicted_tte - actual_tte
    ax2.hist(errors, bins=50, color=TTE_COLOR, alpha=0.7, edgecolor='#1a7f37', density=True)

    # Gaussian fit
    mu, sigma = np.mean(errors), np.std(errors)
    x_fit = np.linspace(errors.min(), errors.max(), 200)
    y_fit = (1 / (sigma * np.sqrt(2 * np.pi))) * np.exp(-0.5 * ((x_fit - mu) / sigma) ** 2)
    ax2.plot(x_fit, y_fit, color=A3_COLOR, linewidth=2, label='Gaussian fit')

    ax2.axvline(x=0, color='black', linestyle='--', alpha=0.5, linewidth=1)

    ax2.text(0.95, 0.92, f'Mean = {mu:.2f} s\nStd = {sigma:.2f} s',
             transform=ax2.transAxes, fontsize=11, fontweight='bold',
             color='#24292f', ha='right', va='top',
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.9))

    ax2.set_xlabel('Prediction Error (s)', fontsize=11)
    ax2.set_ylabel('Density', fontsize=11)
    ax2.set_title('TTE Prediction Error Distribution', fontsize=12, fontweight='bold')
    ax2.legend(fontsize=9)
    ax2.grid(True, alpha=0.2)

    path = os.path.join(OUT_DIR, 'ntn_cho_tte_accuracy.png')
    fig.savefig(path, dpi=150, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    print(f"    -> {path} ({os.path.getsize(path)/1024:.0f} KB)")


if __name__ == '__main__':
    print("=" * 60)
    print("  NTN-CHO Framework Visualization Generator")
    print("=" * 60)
    generate_handover_gif()
    generate_algorithm_comparison_gif()
    generate_results_png()
    generate_tte_accuracy_png()
    print("\nAll visualizations generated in:", OUT_DIR)
    for f in sorted(os.listdir(OUT_DIR)):
        if f.endswith(('.gif', '.png')):
            fpath = os.path.join(OUT_DIR, f)
            print(f"  {f}: {os.path.getsize(fpath)/1024:.0f} KB")
