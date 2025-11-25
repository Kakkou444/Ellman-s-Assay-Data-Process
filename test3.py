import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# ----------------------------------------------------------
# 1. Load data
# ----------------------------------------------------------
# Make sure this path points to your CSV file
df = pd.read_csv("OMe_no_ester_result.csv")

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

# Fit the 4PL curve（仍然在原始浓度上拟合）
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
# Generate smooth x values for the curve（原始浓度）
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

# === 关键改动：用 log10(concentration) 作为 x 轴 ===
xdata_log = np.log10(xdata)
xfit_log  = np.log10(xfit)
ic50_log  = np.log10(ic50) if not np.isnan(ic50) else np.nan

# 1. 画数据点+误差棒（x 用 log10）
ax.errorbar(
    xdata_log, ydata, yerr=yerr,
    fmt="o",
    color="#00DED5",
    ecolor="#00DED5",
    markersize=5,
    capsize=3,
    label="Mean ± SD (% activity)"
)

# 2. 画 4PL 拟合曲线（x 用 log10）
ax.plot(xfit_log, yfit, color="#09D909", linewidth=1.3, label="4PL fit")

# 3. 设置坐标轴尺度/标题（不再用对数坐标，而是 log(浓度) 作为变量）
ax.set_xlabel("log(Drug concentration (µM))")
ax.set_ylabel("% Activity")
ax.set_title("Dose-response Curve")

# 4. 让 matplotlib 自动决定一个合适的范围
plt.autoscale()

# 5. 如果你想竖线从 0 开始，这里把 y 轴下限固定为 0
x_min, x_max = ax.get_xlim()
y_min, y_max = ax.get_ylim()

y_min = 0                      # 竖线从 0 开始
ax.set_xlim(x_min, x_max)
ax.set_ylim(y_min, y_max)

# 6. 画 IC50 辅助线（用 log10(IC50)）
target = 50    # 50% Activity（你也可以改成别的）

if not np.isnan(ic50_log):
    # 横线：从最左边到 IC50 的 log 值
    ax.hlines(
        target,              # y = 50
        x_min, ic50_log,
        linestyles="dashed",         # 从最左边到 IC50
        colors="lightblue",
        linewidth=1
    )

    # 竖线：从 y=0 画到 y=50，x 为 log10(IC50)
    ax.vlines(
        ic50_log,            # x = log10(IC50)
        y_min, target,       # 从 0 到 50
        linestyles="dashed",
        colors="lightblue",
        linewidth=1
    )

# IC50 黑点 + 文字（保持原来的位置，只是数值显示还是 µM）
ax.text(
    0.95, 0.05,                      # ← 相对坐标：右下角
    f"IC50 = {ic50:.2f} µM",
    transform=ax.transAxes,         # ← 使用坐标轴的相对坐标
    ha="right", va="bottom",
    fontsize=10
)

# 不再使用对数坐标，只改 label
plt.xlabel("log(Drug concentration (µM))")
plt.ylabel("% Activity")
plt.title("OMe-Tacrine Derivative Dose-response Curve")
plt.legend()
plt.tight_layout()
plt.show()
