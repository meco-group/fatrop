#!/usr/bin/env python3
"""
Visualize the closed-loop quadrotor MPC recovery logged by quadrotor_mpc_example.

Run the C++ example first (it writes quadrotor_mpc.csv in the working directory), then:

    python3 quadrotor_mpc_visualize.py [path/to/quadrotor_mpc.csv]

Produces:
  * quadrotor_mpc_summary.png  -- time histories (tilt, body rates, position, velocity, thrusts)
  * quadrotor_mpc_3d.png       -- 3D flight path with body-frame triads sampled along it
  * quadrotor_mpc.gif          -- animation of the recovering body frame (if pillow is available)
"""
import sys
import os
import csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401


def quat_to_R(q):
    """Hamilton [w,x,y,z] unit quaternion -> rotation matrix (body->world)."""
    w, x, y, z = q
    return np.array([
        [1 - 2 * (y * y + z * z), 2 * (x * y - z * w),     2 * (x * z + y * w)],
        [2 * (x * y + z * w),     1 - 2 * (x * x + z * z), 2 * (y * z - x * w)],
        [2 * (x * z - y * w),     2 * (y * z + x * w),     1 - 2 * (x * x + y * y)],
    ])


def load(path):
    rows = []
    with open(path) as f:
        for r in csv.DictReader(f):
            rows.append(r)
    goal = None
    data = {k: [] for k in
            ["t", "px", "py", "pz", "vx", "vy", "vz",
             "qw", "qx", "qy", "qz", "wx", "wy", "wz",
             "f0", "f1", "f2", "f3", "tilt_deg", "omega_norm"]}
    for r in rows:
        if int(r["step"]) == -1:
            goal = r
            continue
        for k in data:
            v = r[k]
            data[k].append(float(v) if v != "" else np.nan)
    for k in data:
        data[k] = np.array(data[k])
    return data, goal


G = 9.81
HOVER_PER_ROTOR = G / 4.0       # ~2.4525 N
F_MAX_PER_ROTOR = 3.0 * G / 4.0 # matches f_max in quadrotor_mpc_example.cpp (TWR 3.0)


def plot_summary(d, out):
    t = d["t"]
    fig, ax = plt.subplots(2, 3, figsize=(15, 8))
    fig.suptitle("Quadrotor SO(3) MPC — recovery from a heavy disturbance", fontsize=14)

    ax[0, 0].plot(t, d["tilt_deg"], color="crimson", lw=2)
    ax[0, 0].set_title("Attitude error |Log(q_t⁻¹q)|")
    ax[0, 0].set_ylabel("tilt [deg]"); ax[0, 0].axhline(0, color="k", lw=.5)

    ax[0, 1].plot(t, d["wx"], label="ωx")
    ax[0, 1].plot(t, d["wy"], label="ωy")
    ax[0, 1].plot(t, d["wz"], label="ωz")
    ax[0, 1].plot(t, d["omega_norm"], "k--", lw=1.2, label="|ω|")
    ax[0, 1].set_title("Body angular rate"); ax[0, 1].set_ylabel("[rad/s]")
    ax[0, 1].legend(ncol=2, fontsize=8)

    # per-rotor thrusts with the actuator box [0, f_max]
    for k, lab in zip(("f0", "f1", "f2", "f3"), ("f₀", "f₁", "f₂", "f₃")):
        ax[0, 2].plot(t, d[k], label=lab)
    ax[0, 2].axhline(HOVER_PER_ROTOR, color="gray", ls=":", lw=1, label="hover")
    ax[0, 2].axhline(F_MAX_PER_ROTOR, color="red", ls="--", lw=1, label="bounds")
    ax[0, 2].axhline(0.0, color="red", ls="--", lw=1)
    ax[0, 2].set_title("Commanded rotor thrusts (bounded)"); ax[0, 2].set_ylabel("[N]")
    ax[0, 2].legend(ncol=3, fontsize=8)

    pnorm = np.sqrt(d["px"] ** 2 + d["py"] ** 2 + d["pz"] ** 2)
    ax[1, 0].plot(t, d["px"], label="x")
    ax[1, 0].plot(t, d["py"], label="y")
    ax[1, 0].plot(t, d["pz"], label="z")
    ax[1, 0].plot(t, pnorm, "k--", lw=1.0, label="|p|")
    ax[1, 0].set_title("Position"); ax[1, 0].set_ylabel("[m]"); ax[1, 0].axhline(0, color="k", lw=.5)
    ax[1, 0].legend(ncol=2, fontsize=8)

    ax[1, 1].plot(t, d["vx"], label="vx")
    ax[1, 1].plot(t, d["vy"], label="vy")
    ax[1, 1].plot(t, d["vz"], label="vz")
    ax[1, 1].set_title("Velocity"); ax[1, 1].set_ylabel("[m/s]"); ax[1, 1].axhline(0, color="k", lw=.5)
    ax[1, 1].legend(fontsize=8)

    # control input magnitude: total commanded thrust T = Σf_i and Euclidean norm ‖f‖₂
    T_total = d["f0"] + d["f1"] + d["f2"] + d["f3"]
    f_norm = np.sqrt(d["f0"] ** 2 + d["f1"] ** 2 + d["f2"] ** 2 + d["f3"] ** 2)
    ax[1, 2].plot(t, T_total, color="navy", lw=2, label="total thrust ΣfᵢN")
    ax[1, 2].plot(t, f_norm, color="darkorange", lw=1.5, label="‖f‖₂")
    ax[1, 2].axhline(G, color="gray", ls=":", lw=1, label="weight (mg)")
    ax[1, 2].axhline(4 * F_MAX_PER_ROTOR, color="red", ls="--", lw=1, label="max ΣfᵢN")
    ax[1, 2].set_title("Control input magnitude"); ax[1, 2].set_ylabel("[N]")
    ax[1, 2].legend(fontsize=8)

    for a in ax.ravel():
        a.set_xlabel("time [s]"); a.grid(alpha=.3)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(out, dpi=130)
    print("wrote", out)


# Body-frame rotor positions for the "+ configuration" used by the model:
#   rotor 1 at -x_b, rotor 2 at +y_b, rotor 3 at +x_b, rotor 4 at -y_b
# (so tau_x = L(f2-f4), tau_y = L(f3-f1)).  Index order [f0,f1,f2,f3] -> [1,2,3,4].
_ARM_DIRS = np.array([[-1, 0, 0], [0, 1, 0], [1, 0, 0], [0, -1, 0]], dtype=float)
_ROTOR_COLORS = ["#1f77b4", "#1f77b4", "#d62728", "#d62728"]  # front pair red, rear pair blue


def draw_quad(ax, p, R, arm, thrusts=None, alpha=1.0, base_lw=2.5):
    """Draw a + configuration quadrotor: 4 arms, rotor disks, and a thrust-up arrow.

    If `thrusts` (4 rotor forces) is given, each rotor disk is scaled by its thrust so the
    control-input magnitude is visible directly on the body."""
    tips = [p + R @ (arm * dir_) for dir_ in _ARM_DIRS]
    # arms (cross)
    for tip in tips:
        ax.plot([p[0], tip[0]], [p[1], tip[1]], [p[2], tip[2]],
                color="0.25", lw=base_lw, alpha=alpha, solid_capstyle="round")
    # rotor disks (circles in the body xy-plane), radius encodes thrust magnitude
    th = np.linspace(0, 2 * np.pi, 24)
    bx, by = R[:, 0], R[:, 1]
    for i, tip in enumerate(tips):
        if thrusts is not None:
            frac = np.clip(thrusts[i] / F_MAX_PER_ROTOR, 0.0, 1.0)
            r = arm * (0.28 + 0.42 * frac)
        else:
            r = arm * 0.38
        circ = np.array([tip + r * (np.cos(a) * bx + np.sin(a) * by) for a in th])
        ax.plot(circ[:, 0], circ[:, 1], circ[:, 2],
                color=_ROTOR_COLORS[i], lw=1.8, alpha=alpha)
    # body +z (thrust direction) arrow
    up = R[:, 2] * arm * 1.4
    ax.plot([p[0], p[0] + up[0]], [p[1], p[1] + up[1]], [p[2], p[2] + up[2]],
            color="green", lw=base_lw * 0.8, alpha=alpha)


def plot_3d(d, out):
    P = np.vstack([d["px"], d["py"], d["pz"]]).T
    Q = np.vstack([d["qw"], d["qx"], d["qy"], d["qz"]]).T
    F = np.vstack([d["f0"], d["f1"], d["f2"], d["f3"]]).T
    fig = plt.figure(figsize=(9, 8))
    ax = fig.add_subplot(111, projection="3d")
    ax.plot(P[:, 0], P[:, 1], P[:, 2], color="0.5", lw=1.5, label="flight path")
    ax.scatter(*P[0], color="crimson", s=60, label="start (disturbed)")
    ax.scatter(*P[-1], color="green", s=60, label="end (hover)")
    ax.scatter(0, 0, 0, marker="*", color="black", s=160, label="hover target")

    span = max(np.ptp(P[:, 0]), np.ptp(P[:, 1]), np.ptp(P[:, 2]), 1e-3)
    arm = 0.10 * span + 0.05
    n = len(P)
    idx = list(range(0, n, max(1, n // 12)))
    if n - 1 not in idx:
        idx.append(n - 1)
    for i in idx:
        a = 0.30 + 0.70 * i / (n - 1)
        draw_quad(ax, P[i], quat_to_R(Q[i]), arm, thrusts=F[i], alpha=a)

    ax.set_xlabel("x [m]"); ax.set_ylabel("y [m]"); ax.set_zlabel("z [m]")
    ax.set_title("Quadrotor recovery — rotor disk size ∝ thrust, green = body-up")
    ax.legend(loc="upper left", fontsize=9)
    _set_equal(ax, P)
    fig.tight_layout()
    fig.savefig(out, dpi=130)
    print("wrote", out)


def _set_equal(ax, P):
    mins = P.min(0); maxs = P.max(0)
    c = 0.5 * (mins + maxs)
    r = 0.5 * max((maxs - mins).max(), 0.3)
    ax.set_xlim(c[0] - r, c[0] + r)
    ax.set_ylim(c[1] - r, c[1] + r)
    ax.set_zlim(c[2] - r, c[2] + r)


def animate(d, out):
    try:
        from matplotlib.animation import FuncAnimation, PillowWriter
    except Exception as e:
        print("skip animation:", e)
        return
    P = np.vstack([d["px"], d["py"], d["pz"]]).T
    Q = np.vstack([d["qw"], d["qx"], d["qy"], d["qz"]]).T
    F = np.vstack([d["f0"], d["f1"], d["f2"], d["f3"]]).T
    n = len(P)
    fig = plt.figure(figsize=(8, 7))
    ax = fig.add_subplot(111, projection="3d")
    _set_equal(ax, P)
    ax.set_xlabel("x [m]"); ax.set_ylabel("y [m]"); ax.set_zlabel("z [m]")
    span = max(np.ptp(P[:, 0]), np.ptp(P[:, 1]), np.ptp(P[:, 2]), 1e-3)
    arm = 0.14 * span + 0.06
    ax.scatter(0, 0, 0, marker="*", color="black", s=160)
    (path_line,) = ax.plot([], [], [], color="0.6", lw=1.2)

    def update(i):
        # redraw: remove the previous body (everything but the persistent path + target)
        for art in list(ax.lines):
            if art is not path_line:
                art.remove()
        path_line.set_data(P[:i + 1, 0], P[:i + 1, 1])
        path_line.set_3d_properties(P[:i + 1, 2])
        draw_quad(ax, P[i], quat_to_R(Q[i]), arm, thrusts=F[i])
        ax.set_title(f"MPC recovery  t={d['t'][i]:.2f}s   tilt={d['tilt_deg'][i]:.0f}°")
        return ax.lines

    anim = FuncAnimation(fig, update, frames=n, interval=60, blit=False)
    try:
        anim.save(out, writer=PillowWriter(fps=15))
        print("wrote", out)
    except Exception as e:
        print("skip animation save:", e)
    plt.close(fig)


def main():
    csv_path = sys.argv[1] if len(sys.argv) > 1 else "quadrotor_mpc.csv"
    if not os.path.exists(csv_path):
        sys.exit(f"CSV not found: {csv_path}\nRun the quadrotor_mpc_example first.")
    d, goal = load(csv_path)
    outdir = os.path.dirname(os.path.abspath(csv_path))
    plot_summary(d, os.path.join(outdir, "quadrotor_mpc_summary.png"))
    plot_3d(d, os.path.join(outdir, "quadrotor_mpc_3d.png"))
    animate(d, os.path.join(outdir, "quadrotor_mpc.gif"))


if __name__ == "__main__":
    main()
