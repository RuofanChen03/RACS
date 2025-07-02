import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import plotly.io as pio
import os
import numpy as np

# Helper function to load Excel data
def load_data(filename, sheet_index=0):
    if not os.path.exists(filename):
        raise FileNotFoundError(f"{filename} not found!")
    print(f"Loading file {filename}")
    return pd.read_excel(filename, sheet_name=sheet_index)

# Overlap function
def overlap(x1, x2, y1, y2):
    if x1 > y2:
        return +1.0
    if x2 < y1:
        return -1.0
    ax = x2 - x1
    ay = y2 - y1
    if x1 <= y1:
        return +ax / ay
    if x1 > y1:
        return -ay / ax
    if x1 == y1 and x2 == y2:
        return 0
    return np.nan

# Full comparison with scaffold filtering
def comparison_all(sample1, sample2):
    results = []
    scflds_smpl1 = sample1['Region'].unique()

    for scfld in scflds_smpl1:
        print(scfld)
        sset1 = sample1[sample1['Region'] == scfld].sort_values('start')
        sset2 = sample2[sample2['Region'] == scfld].sort_values('start')

        for _, row1 in sset1.iterrows():
            for _, row2 in sset2.iterrows():
                ov = overlap(row1['start'], row1['end'], row2['start'], row2['end'])
                if isinstance(ov, float):
                    print(row1.name, row2.name)
                    print(ov)
                    results.append([scfld, row1['start'], row1['end'], row2['start'], row2['end'], ov])

    df = pd.DataFrame(results, columns=['scfld', 'x1', 'x2', 'y1', 'y2', 'overlap'])
    return df[df['overlap'].abs() != 1]

# Comparison with logging/debugging
def comparison(sample1, sample2, DBG=True):
    sample1 = sample1.sort_values('Region')
    sample2 = sample2.sort_values('Region')

    scflds = sorted(set(sample1['Region'].unique()).union(sample2['Region'].unique()))
    results = []

    for scfld in scflds:
        print(scfld)
        sset1 = sample1[sample1['Region'] == scfld].sort_values('start')
        sset2 = sample2[sample2['Region'] == scfld].sort_values('start')

        if sset1.empty:
            if DBG:
                print(f"No data for scaffold {scfld} in sample1")
                print(sset1)
                print(sset2)
            for _, row in sset2.iterrows():
                results.append([scfld, "--", "--", row['start'], row['end'], -1])
        elif sset2.empty:
            if DBG:
                print(f"No data for scaffold {scfld} in sample2")
                print(sset1)
                print(sset2)
            for _, row in sset1.iterrows():
                results.append([scfld, row['start'], row['end'], "--", "--", 1])
        else:
            for _, row1 in sset1.iterrows():
                reg_results = []
                for _, row2 in sset2.iterrows():
                    ov = overlap(row1['start'], row1['end'], row2['start'], row2['end'])
                    if not pd.isna(ov):
                        print(scfld, row1.name, row2.name, "---", row1['start'], row1['end'], row2['start'], row2['end'])
                        print(ov)
                        reg_results.append([scfld, row1['start'], row1['end'], row2['start'], row2['end'], ov])

                reg_df = pd.DataFrame(reg_results, columns=['scafold', 'x1', 'x2', 'y1', 'y2', 'overlap'])
                filter_cond = reg_df['overlap'].abs()
                if (filter_cond == 1).sum() == len(sset2):
                    results.append([scfld, row1['start'], row1['end'], "--", "--", reg_df.iloc[0]['overlap']])
                else:
                    for _, r in reg_df[filter_cond != 1].iterrows():
                        results.append([r['scafold'], r['x1'], r['x2'], r['y1'], r['y2'], r['overlap']])

    df_results = pd.DataFrame(results, columns=['scafold', 'x1', 'x2', 'y1', 'y2', 'overlap'])
    df_results['overlap'] = pd.to_numeric(df_results['overlap'], errors='coerce')
    return df_results[df_results['overlap'].abs() != 1]

# Visualization

def viz_diffs(results, filename="", threshold=1000):
    res = results.dropna(subset=['overlap'])
    res = res[res['overlap'] < threshold]
    yvar = res['overlap']

    if filename:
        pio.write_html(px.bar(res, x='scafold', y='overlap'), file=filename.replace('.pdf', '.html'))
    else:
        fig = go.Figure()
        fig.add_trace(go.Bar(x=res['scafold'], y=res['overlap']))
        fig.show()

    # Static histogram
    import matplotlib.pyplot as plt
    plt.figure(figsize=(10, 6))
    plt.hist(yvar, bins=75, color="lightgray", density=True)
    mean_y = np.mean(yvar)
    std_y = np.std(yvar)
    xmin, xmax = plt.xlim()
    xfit = np.linspace(xmin, xmax, 100)
    yfit = (1 / (std_y * np.sqrt(2 * np.pi))) * np.exp(-0.5 * ((xfit - mean_y) / std_y) ** 2)
    plt.plot(xfit, yfit, color="black", lw=2)
    plt.title("RACS vs MACS")
    plt.xlabel("overlap")
    plt.text(0, max(yfit)*0.9, f"mean={mean_y:.2f}  --  sd={std_y:.2f}")
    plt.grid(True)
    if filename:
        plt.savefig(filename)
    else:
        plt.show()
