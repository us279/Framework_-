import numpy as np
import glob
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# Known parameters (must match the Fortran simulation)
a = 2.0
T_final = 0.5

# Load final snapshot from each precision's CSV outputs
precisions = ["fp16", "fp32", "fp64", "fp128"]
final_data = {}
for prec in precisions:
    files = sorted(glob.glob(f"advection_{prec}_*.csv"))
    if not files:
        continue  # skip if no data for this precision
    last_file = files[-1]                        # last file is final state
    data = np.loadtxt(last_file, delimiter=',')  # load x and u columns
    x_vals = data[:, 0]
    u_vals = data[:, 1]
    final_data[prec] = (x_vals, u_vals)

# Ensure we have at least one precision loaded
if len(final_data) == 0:
    raise FileNotFoundError("No advection output files found. Run the Fortran solver first.")

# Use the first loaded dataset's x-grid as reference (all precisions use same grid)
ref_prec = next(iter(final_data))
x = final_data[ref_prec][0]
# Compute exact solution at final time for error analysis
u_exact = np.sin(2 * np.pi * (x - a * T_final))

# Calculate error norms for each precision
error_norms = {"L1": [], "L2": [], "Linf": []}
labels = []
for prec, (x_vals, u_vals) in final_data.items():
    err = u_vals - u_exact
    error_norms["L1"].append(np.mean(np.abs(err)))           # L1 norm (mean absolute error)
    error_norms["L2"].append(np.sqrt(np.mean(err**2)))       # L2 norm (RMSE)
    error_norms["Linf"].append(np.max(np.abs(err)))          # L∞ norm (max error)
    labels.append(prec)

# Create bar plot comparing error norms across precisions
labels = np.array(labels)
x_pos = np.arange(len(labels))
bar_width = 0.2
fig, ax = plt.subplots(figsize=(6,4))
for i, norm in enumerate(["L1", "L2", "Linf"]):
    ax.bar(x_pos + (i-1)*bar_width, error_norms[norm], bar_width, label=norm)
ax.set_xticks(x_pos, labels=labels)
ax.set_ylabel("Error Norm Value")
ax.set.xlabel("Floating-Point Precision")
ax.set_title("Final Solution Error vs. Exact Solution")
ax.legend(title="Norm")
ax.grid(axis='y', linestyle='--', alpha=0.5)
plt.tight_layout()
plt.savefig("advection_error_norms.png")   # save error comparison plot
plt.close()

# Prepare animation of solution evolution using the highest available precision data
# (Assume the highest precision run is the most accurate representation)
if "fp128" in final_data:
    anim_prec = "fp128"
elif "fp64" in final_data:
    anim_prec = "fp64"
elif "fp32" in final_data:
    anim_prec = "fp32"
else:
    anim_prec = list(final_data.keys())[0]

# Load all snapshot files for the chosen precision
snap_files = sorted(glob.glob(f"advection_{anim_prec}_*.csv"))
solution_frames = [np.loadtxt(f, delimiter=',')[:, 1] for f in snap_files]  # list of u arrays
time_steps = len(solution_frames)
time_values = np.linspace(0, T_final, time_steps)
x_grid = np.loadtxt(snap_files[0], delimiter=',')[:, 0]

# Set up the figure and animation writer
fig, ax = plt.subplots(figsize=(6,4))
line, = ax.plot([], [], 'b-')
time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)
ax.set_xlim(x_grid[0], x_grid[-1])
ax.set_ylim(np.min(solution_frames)*1.1, np.max(solution_frames)*1.1)
ax.set_xlabel("x")
ax.set_ylabel("u(x,t)")

def init():
    """Initialize animation frame."""
    line.set_data(x_grid, solution_frames[0])
    time_text.set_text(f"t = {0:.2f}")
    return line, time_text

def update(frame):
    """Update animation for frame index."""
    line.set_data(x_grid, solution_frames[frame])
    time_text.set_text(f"t = {time_values[frame]:.2f}")
    return line, time_text

ani = animation.FuncAnimation(fig, update, frames=time_steps, init_func=init, blit=True)
ani.save("advection_evolution.gif", writer='pillow', fps=10)
plt.close()
print("Error norms (L1, L2, L_inf):", error_norms)
print(f"Animation saved for precision: {anim_prec}")
