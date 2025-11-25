import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import os

input_file = "Your_file.csv"  # Change to your result file name
base_name = os.path.splitext(input_file)[0]

# ----------------------------------------------------------
# 1. Load data
# ----------------------------------------------------------
# Make sure this path points to your CSV file
df = pd.read_csv(input_file)

# Expected columns in your file (as in the one you uploaded):
# 'Column' = well ID (e.g. 'A3')
# 'Rate_per_second'
# 'Percent_activity'
# 'Percent_inhibition'

# ----------------------------------------------------------
# 2. Add row letter and column number from the well ID
# ----------------------------------------------------------
df["Row"] = df["Column"].str[0]
df["ColNum"] = df["Column"].str[1:].astype(int)

# ----------------------------------------------------------
# 3. Map plate column → concentration (µM)
# ----------------------------------------------------------
col_to_conc = {
    3: 50,
    4: 25,
    5: 12.5,
    6: 6.25,
    7: 3.125,
    8: 1.5625,
    9: 0.78125,
}

# Keep only A–C rows and columns 3–9 (your triplicates)
mask = df["Row"].isin(["A", "B", "C"]) & df["ColNum"].isin(col_to_conc.keys())
replicates = df[mask].copy()

# Add concentration column
replicates["Concentration"] = replicates["ColNum"].map(col_to_conc).astype(float)

# ----------------------------------------------------------
# 4. Calculate mean and SD of % inhibition at each concentration
# ----------------------------------------------------------
summary = (
    replicates
    .groupby("Concentration")["Percent_activity"]
    .agg(["mean", "std"])
    .reset_index()
    .sort_values("Concentration", ascending=False)   # high → low conc
)

print("Mean ± SD of % activity for each concentration:")
print(summary)
summary_out = f"{base_name}_summary.csv"
summary.to_csv(summary_out, index=False)
print(f"Saved summary to: {summary_out}")


# ----------------------------------------------------------
# 5. Define 4-parameter logistic (4PL) function
#    y = D + (A - D) / (1 + (x/C)^B)
#    A = response at low x (top or bottom depending on trend)
#    D = response at high x
#    C = EC50/IC50-like parameter
#    B = Hill slope
# ----------------------------------------------------------
def four_pl(x, A, B, C, D):
    return D + (A - D) / (1.0 + (x / C) ** B)

xdata = summary["Concentration"].values
ydata = summary["mean"].values
yerr  = summary["std"].values

# Initial parameter guesses:
A0 = ydata.min()          # one end of the curve
D0 = ydata.max()          # other end
C0 = np.median(xdata)     # mid concentration
B0 = 1.0                  # slope

p0 = [A0, B0, C0, D0]

# Fit the 4PL curve
params, cov = curve_fit(four_pl, xdata, ydata, p0=p0, maxfev=10000)
A_fit, B_fit, C_fit, D_fit = params
print("\n4PL fit parameters:")
print(f"A (low response)  = {A_fit:.3f}")
print(f"B (Hill slope)    = {B_fit:.3f}")
print(f"C (IC50-like)     = {C_fit:.3f} µM")
print(f"D (high response) = {D_fit:.3f}")

# ----------------------------------------------------------
# 6. Plot data (mean ± SD) and fitted 4PL curve
# ----------------------------------------------------------
# Generate smooth x values for the curve
xfit = np.logspace(np.log10(xdata.min()), np.log10(xdata.max()), 200)
yfit = four_pl(xfit, *params)

# ----------------------------------------------------------
# 7. Calculate concentration at % inhibition = 50 (IC50-like)
# ----------------------------------------------------------

target_inhibition = 50  # 50%

A_fit, B_fit, C_fit, D_fit = params

# Prevent invalid values if the 4PL curve does not cross 50%
if (A_fit - D_fit) * (target_inhibition - D_fit) <= 0:
    ic50 = np.nan
    print("\nThe fitted curve does not cross 50% inhibition.")
else:
    ic50 = C_fit * ((A_fit - D_fit) / (target_inhibition - D_fit) - 1) ** (1 / B_fit)
    print(f"\nEstimated concentration at 50% inhibition (IC50-like): {ic50:.4f} µM")

plt.figure(figsize=(6, 5))
ax = plt.gca()

xdata_log = np.log10(xdata)
xfit_log  = np.log10(xfit)
ic50_log  = np.log10(ic50) if not np.isnan(ic50) else np.nan

ax.errorbar(
    xdata_log, ydata, yerr=yerr,
    fmt="o",
    color="#00DED5",
    ecolor="#00DED5",
    markersize=5,
    capsize=3,
    label="Mean ± SD (% activity)"
)

ax.plot(xfit_log, yfit, color="#09D909", linewidth=1.3, label="4PL fit")

ax.set_xlabel("log(Drug concentration (µM))")
ax.set_ylabel("% Activity")
ax.set_title("Dose-response Curve")

plt.autoscale()

x_min, x_max = ax.get_xlim()
y_min, y_max = ax.get_ylim()

y_min = 0                      
ax.set_xlim(x_min, x_max)
ax.set_ylim(y_min, y_max)

target = 50

if not np.isnan(ic50_log):
    ax.hlines(
        target,              # y = 50
        x_min, ic50_log,
        linestyles="dashed",         
        colors="lightblue",
        linewidth=1
    )

    ax.vlines(
        ic50_log,            
        y_min, target,       
        linestyles="dashed",
        colors="lightblue",
        linewidth=1
    )

ax.text(
    0.95, 0.05,                      
    f"IC50 = {ic50:.2f} µM",
    transform=ax.transAxes,         
    ha="right", va="bottom",
    fontsize=10
)

plt.xlabel("log(Drug concentration (µM))")
plt.ylabel("% Activity")
plt.title("Dose-response Curve")
plt.legend()
plt.tight_layout()
figure_out = f"{base_name}_plot.png"
plt.savefig(figure_out, dpi=300, bbox_inches='tight')
print(f"Saved figure to: {figure_out}")
plt.show()
