import os
import pandas as pd
import numpy as np
import datetime


def parse_time_seconds(t):
    if isinstance(t.iloc[0], datetime.time):
        secs = t.apply(
            lambda x: x.hour * 3600 + x.minute * 60 + x.second
            if isinstance(x, datetime.time) else np.nan
        )
        return secs - np.nanmin(secs)

    t_td = pd.to_timedelta(t, errors='coerce')
    if not t_td.isna().all():
        x = t_td.dt.total_seconds().astype(float)
        return x - np.nanmin(x)

    t_dt = pd.to_datetime(t, errors='coerce')
    if not t_dt.isna().all():
        x = (t_dt - t_dt.min()).dt.total_seconds().astype(float)
        return x

    t_num = pd.to_numeric(t, errors='coerce')
    if t_num.notna().any():
        x_sec = t_num * 86400.0
        return x_sec - np.nanmin(x_sec)

    raise ValueError("Time column cannot be parsed as time/timedelta/datetime/float")


def analyse_plate(df, control_col="A1", no_enzyme_col="A2"):
    time_col = df.columns[0]
    data_cols = df.columns[2:]

    t_raw = df[time_col]
    x = parse_time_seconds(t_raw)

    if not np.isfinite(x).any() or np.nanstd(x) == 0:
        x = np.arange(len(df), dtype=float)

    x_mask = np.isfinite(x)
    x_mean = np.nanmean(x[x_mask])
    x_std = np.nanstd(x[x_mask])
    if x_std == 0:
        x_std = 1.0
    x_stdized = (x - x_mean) / x_std

    slopes = {}
    intercepts = {}

    for col in data_cols:
        y = pd.to_numeric(df[col], errors='coerce').astype(float)
        mask = np.isfinite(x_stdized) & np.isfinite(y)
        x2, y2 = x_stdized[mask], y[mask]

        if len(x2) < 2 or np.allclose(y2, y2[0], rtol=0, atol=1e-12):
            continue

        # y = m_std * ((x - mean)/std) + b
        m_std, b = np.polyfit(x2, y2, 1)
        m_per_sec = m_std / x_std
        slopes[col] = m_per_sec
        intercepts[col] = b - m_std * (x_mean / x_std)

    percent_activity = {}
    control_slope = slopes.get(control_col, np.nan)

    if not np.isfinite(control_slope) or control_slope == 0:
        print("[warning] unable to calculate, control well slope = 0 or missing")
    else:
        for col, m in slopes.items():
            if col == no_enzyme_col:
                continue
            if np.isfinite(m):
                percent_activity[col] = 100.0 * (m / control_slope)

    percent_inhibition = {}
    for col, pa in percent_activity.items():
        if np.isfinite(pa):
            percent_inhibition[col] = 100.0 - pa

    cols_out, rates_out, pct_out, inhib_out = [], [], [], []

    for col, m in slopes.items():
        if col == no_enzyme_col:
            continue
        cols_out.append(col)
        rates_out.append(m)
        pa = percent_activity.get(col, np.nan)
        pct_out.append(pa)
        inhib_out.append(100.0 - pa if np.isfinite(pa) else np.nan)

    result_df = pd.DataFrame({
        "Column": cols_out,
        "Rate_per_second": rates_out,
        "Percent_activity": pct_out,
        "Percent_inhibition": inhib_out,
    })

    return result_df


def process_file(filepath, control_col="A1", no_enzyme_col="A2"):

    output_folder = "bulkresults"
    os.makedirs(output_folder, exist_ok=True)

    print(f"\nProcessing: {filepath}")
    ext = os.path.splitext(filepath)[1].lower()

    if ext in [".xlsx", ".xls"]:
        df = pd.read_excel(filepath)
    elif ext == ".csv":
        df = pd.read_csv(filepath)
    else:
        print("skip: file type is not supported")
        return

    result_df = analyse_plate(df, control_col=control_col, no_enzyme_col=no_enzyme_col)

    base, _ = os.path.splitext(os.path.basename(filepath))
    out_path = os.path.join(output_folder, f"{base}_results.csv")

    result_df.to_csv(out_path, index=False, encoding="utf-8-sig")
    print(f"Saved: {out_path}")


if __name__ == "__main__":
    folder = "bulkdata"   # DATA GO IN TO THIS FOLDER
    for name in os.listdir(folder):
        path = os.path.join(folder, name)
        if os.path.isfile(path):
            process_file(path)
