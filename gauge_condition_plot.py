import os

import numpy as np
import matplotlib.pyplot as plt

input_path = "Data/global_quantities_gauge.dat"
output_dir = os.path.dirname(input_path)

data = np.loadtxt(input_path)

# ============================================================
# What this diagnostic actually is:
#
# The evolution equations amount to a set of constraint
# equations that the numerically-evolved fields should satisfy
# exactly if the simulation were perfect. Every output step,
# the code measures how badly each constraint equation is
# violated at every grid point, then reports two summaries per
# equation:
#
#   - integrated |violation|: summed over the whole grid
#   - max |violation|: the single worst grid point
#
# This is purely a numerical-health check, not physics you're
# trying to extract. If it stays small/flat (or decays), the
# simulation is respecting its own equations of motion. If it
# grows - especially the max - that's a sign of growing
# numerical error, under-resolution, or instability.
#
# ============================================================
# global_quantities_gauge.dat column layout (one row per output):
#
# 0 = timestep
#
# followed by
# [integrated |violation| for each constraint equation]
# followed by
# [max |violation| for each constraint equation]
#
# The number of constraint equations is set by the physics
# model at runtime (not fixed in the config file), so it's
# inferred here from the number of remaining columns rather
# than being hard-coded.
# ============================================================

timestep = data[:, 0]

num_columns = data.shape[1] - 1
num_equations = num_columns // 2

integrated_violation = data[:, 1 : 1 + num_equations]
max_violation = data[:, 1 + num_equations : 1 + 2 * num_equations]

# Summed across all constraint equations - this is the single
# number you actually want to keep an eye on: is the overall
# constraint violation staying small, or growing over time?
total_integrated_violation = integrated_violation.sum(axis=1)
total_max_violation = max_violation.max(axis=1)

# ============================================================
# Main plot: overall constraint violation vs timestep
# ============================================================

fig, ax1 = plt.subplots(figsize=(8, 6))

ax1.plot(
    timestep,
    total_integrated_violation,
    color="tab:blue",
    label="Total integrated |violation|"
)
ax1.set_xlabel("Timestep")
ax1.set_ylabel("Total integrated |violation|", color="tab:blue")
ax1.set_yscale("log")
ax1.tick_params(axis="y", labelcolor="tab:blue")

ax2 = ax1.twinx()
ax2.plot(
    timestep,
    total_max_violation,
    color="tab:red",
    label="Worst-case max |violation|"
)
ax2.set_ylabel("Worst-case max |violation|", color="tab:red")
ax2.set_yscale("log")
ax2.tick_params(axis="y", labelcolor="tab:red")

plt.title("Overall gauge condition violation vs timestep")
fig.tight_layout()

plt.savefig(
    os.path.join(output_dir, "gauge_condition_overall_vs_time.png"),
    dpi=300
)

plt.close()

# ============================================================
# Secondary plot: per-equation breakdown, in case one equation
# is behaving much worse than the others
# ============================================================

fig, (ax_top, ax_bottom) = plt.subplots(2, 1, figsize=(8, 9), sharex=True)

for eq_iter in range(num_equations):
    ax_top.plot(
        timestep,
        integrated_violation[:, eq_iter],
        label=f"Equation {eq_iter}"
    )

ax_top.set_ylabel("Integrated |violation|")
ax_top.set_yscale("log")
ax_top.set_title("Per-equation gauge condition violation vs timestep")
ax_top.legend()

for eq_iter in range(num_equations):
    ax_bottom.plot(
        timestep,
        max_violation[:, eq_iter],
        label=f"Equation {eq_iter}"
    )

ax_bottom.set_xlabel("Timestep")
ax_bottom.set_ylabel("Max |violation|")
ax_bottom.set_yscale("log")

fig.tight_layout()

plt.savefig(
    os.path.join(output_dir, "gauge_condition_per_equation_vs_time.png"),
    dpi=300
)

plt.close()

print("All gauge-condition-vs-timestep plots generated.")
