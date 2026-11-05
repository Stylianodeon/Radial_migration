#!/usr/bin/env python3
"""Animate snapshots exported by LagrangeSnapshotsExec (no orbit integration)."""
import argparse
from pathlib import Path
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.lines import Line2D
from matplotlib.transforms import Affine2D
import numpy as np


def read_points(path):
    rows = []
    ring = None
    for line in path.read_text().splitlines():
        match = re.search(r"equilibrium ring: R=([^ ]+)", line)
        if match:
            ring = float(match.group(1))
        if line.strip() and not line.startswith("#"):
            rows.append([float(value) for value in line.split()])
    data = np.array(rows, dtype=float).reshape(-1, 3)
    if not np.isfinite(data).all():
        raise ValueError(f"Non-finite equilibrium coordinates or residuals in {path}")
    return data[:, :2], ring


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path, help="Directory exported by LagrangeSnapshotsExec")
    parser.add_argument("--fps", type=int, default=10, help="Animation frames per second (default: 10)")
    args = parser.parse_args()
    if args.fps <= 0:
        parser.error("--fps must be positive")
    directory = args.directory
    omega, vc, extent, spacing = np.loadtxt(directory / "metadata.dat")
    frames = np.atleast_2d(np.loadtxt(directory / "frames.dat"))
    samples = [read_points(directory / f"points_{int(frame):04d}.dat") for frame in frames[:, 0]]
    grid = np.loadtxt(directory / "spiral.dat")
    nx = len(np.unique(grid[:, 0]))
    ny = len(np.unique(grid[:, 1]))
    x, y, phi = (grid[:, i].reshape(ny, nx) for i in range(3))
    phi = np.ma.masked_where(~np.isfinite(phi) | (np.hypot(x, y) > extent), phi)
    vmax = max(float(np.ma.max(np.ma.abs(phi))), np.finfo(float).eps)
    rcor = abs(vc / omega) if omega else np.nan
    shifts = np.full((len(frames), 2), np.nan)
    for i, (points, _) in enumerate(samples):
        if len(points) and np.isfinite(rcor):
            radii = np.hypot(points[:, 0], points[:, 1])
            shifts[i] = (np.array([radii.min(), radii.max()]) - rcor) * 1000.0

    plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False})
    fig = plt.figure(figsize=(12.5, 5.8), layout="constrained")
    gs = fig.add_gridspec(2, 3, height_ratios=[1, .06], width_ratios=[1, 1, .85])
    axes = [fig.add_subplot(gs[0, i]) for i in range(3)]
    images, markers, rings = [], [], []
    angles = np.linspace(0, 2*np.pi, 361)
    bounds = [-extent-spacing/2, extent+spacing/2, -extent-spacing/2, extent+spacing/2]
    for ax, title in zip(axes[:2], ["Inertial frame", "Frame rotating with the spiral"]):
        images.append(ax.imshow(phi, origin="lower", extent=bounds, cmap="RdBu_r",
                                vmin=-vmax, vmax=vmax, interpolation="bilinear"))
        markers.append(ax.scatter([], [], s=48, c="#ffda48", edgecolors="#151515", linewidths=1, zorder=4))
        ring, = ax.plot([], [], "--", color="#555555", linewidth=1.1, zorder=3)
        rings.append(ring)
        ax.set(xlim=(-extent, extent), ylim=(-extent, extent), xlabel="x [kpc]", ylabel="y [kpc]", title=title)
        ax.set_aspect("equal")
    colorbar = fig.colorbar(images[1], cax=fig.add_subplot(gs[1, :2]), orientation="horizontal")
    colorbar.set_label(r"Spiral contribution to $\Phi_{\rm eff}$ [kpc$^2$/Myr$^2$] — fixed colour scale")
    colorbar.formatter.set_powerlimits((0, 0))
    colorbar.update_ticks()
    ax = axes[2]
    ax.plot(frames[:, 1], shifts[:, 0], color="#2166ac", label="Innermost equilibria")
    ax.plot(frames[:, 1], shifts[:, 1], color="#b2182b", label="Outermost equilibria")
    ax.axhline(0, color="0.65", linestyle="--", linewidth=.8)
    cursor = ax.axvline(frames[0, 1], color="0.25", linewidth=1)
    dots, = ax.plot([], [], "o", color="#222222", markersize=5)
    ax.set(xlim=(frames[0, 1], frames[-1, 1]), xlabel="Time [Myr]",
           ylabel=r"$R - R_{\rm cor}$ [pc]", title="Small radial shifts")
    ax.legend(loc="upper right", fontsize=8)
    ax.grid(alpha=.15)
    handles = [Line2D([], [], marker="o", color="none", markerfacecolor="#ffda48",
                      markeredgecolor="#151515", label="Instantaneous equilibrium"),
               Line2D([], [], color="#555555", linestyle="--", label="Corotation radius / ring when spiral is off")]
    axes[0].legend(handles=handles, loc="lower left", fontsize=7)
    caption = fig.add_subplot(gs[1, 2])
    caption.axis("off")
    caption.text(.5, .4, "Equilibrium locations, not stellar trajectories", ha="center", va="center", fontsize=8)
    heading = fig.suptitle("")

    def update(index):
        _, time, scale = frames[index]
        points, ring = samples[index]
        angle = (omega * time) % (2*np.pi)
        c, s = np.cos(angle), np.sin(angle)
        # x_rot = c*x + s*y; y_rot = -s*x + c*y.
        rotated = points @ np.array([[c, -s], [s, c]])
        markers[0].set_offsets(points)
        markers[1].set_offsets(rotated)
        images[0].set_data(phi * scale)
        images[0].set_transform(Affine2D().rotate(angle) + axes[0].transData)
        images[1].set_data(phi * scale)
        radius = ring if ring is not None else rcor
        for artist in rings:
            artist.set_data(radius*np.cos(angles), radius*np.sin(angles))
            artist.set_linewidth(1.8 if ring is not None else .8)
        cursor.set_xdata([time, time])
        dots.set_data([time, time], shifts[index])
        state = "equilibrium ring" if ring is not None else f"{len(points)} isolated points"
        heading.set_text(f"Spiral equilibria  |  t = {time:.0f} Myr  |  amplitude = {scale:.2f}  |  {state}")
        return [*markers, *images, *rings, cursor, dots, heading]

    peak = int(np.argmax(frames[:, 2]))
    update(peak)
    preview = directory / "lagrange_peak.png"
    fig.savefig(preview, dpi=140)
    animation = FuncAnimation(fig, update, frames=len(frames), interval=1000/args.fps, blit=False)
    output = directory / "lagrange_motion.gif"
    animation.save(output, writer=PillowWriter(fps=args.fps), dpi=90)
    plt.close(fig)
    print(f"Saved {output} and {preview}")


if __name__ == "__main__":
    main()
