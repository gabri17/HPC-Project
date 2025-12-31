import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import matplotlib.ticker as ticker

# --- CONFIGURATION ---
INPUT_FILE = 'benchmark_results_v2.csv'
OUTPUT_DPI = 300
SNS_THEME = "whitegrid"

def load_and_prep_data(filename):
    try:
        df = pd.read_csv(filename)
    except FileNotFoundError:
        print(f"Error: {filename} not found.")
        return None

    # Sort
    df = df.sort_values(by=['Type', 'BufferSize', 'Workers'])

    # --- Pre-calculate Metrics for ALL rows ---
    # We need to find the T(1) baseline for every specific Buffer Size
    
    # 1. Get Baseline Times (Workers=1) for each BufferSize
    # This creates a map: {BufferSize: Time_at_1_Worker}
    strong_baselines = df[(df['Type'] == 'Strong') & (df['Workers'] == 1)].set_index('BufferSize')['Time'].to_dict()
    weak_baselines = df[(df['Type'] == 'Weak') & (df['Workers'] == 1)].set_index('BufferSize')['Time'].to_dict()

    def calculate_metrics(row):
        t1 = 0
        if row['Type'] == 'Strong':
            t1 = strong_baselines.get(row['BufferSize'], row['Time'])
            speedup = t1 / row['Time']
            efficiency = speedup / row['Workers']
        else: # Weak
            t1 = weak_baselines.get(row['BufferSize'], row['Time'])
            # Weak Scaling: Ideal Time is Constant (T1)
            # Efficiency = Ideal_Time / Actual_Time
            efficiency = t1 / row['Time']
            speedup = efficiency * row['Workers'] # Scaled speedup
        return pd.Series([speedup, efficiency])

    df[['Speedup', 'Efficiency']] = df.apply(calculate_metrics, axis=1)
    return df

def get_optimal_buffer(df, scaling_type):
    # Find best buffer at MAX workers
    max_workers = df['Workers'].max()
    subset = df[(df['Type'] == scaling_type) & (df['Workers'] == max_workers)]
    return subset.loc[subset['Time'].idxmin(), 'BufferSize']

# ==========================================
# GRAPH 1: The "Landscape" Heatmap
# ==========================================
def plot_efficiency_heatmap(df):
    subset = df[df['Type'] == 'Strong']
    # Pivot data: Rows=Buffer, Cols=Workers, Values=Efficiency
    matrix = subset.pivot_table(index='BufferSize', columns='Workers', values='Efficiency')
    
    plt.figure(figsize=(12, 8))
    sns.heatmap(matrix, annot=True, fmt=".2f", cmap="RdYlGn", vmin=0, vmax=1.1, 
                cbar_kws={'label': 'Efficiency (1.0 = Perfect)'})
    
    plt.title('Strong Scaling Efficiency Landscape\n(Dark Green is Best)', fontsize=16)
    plt.xlabel('Number of Workers', fontsize=12)
    plt.ylabel('Buffer Size', fontsize=12)
    plt.tight_layout()
    plt.savefig('insight_1_heatmap.png', dpi=OUTPUT_DPI)
    print("Generated insight_1_heatmap.png")

# ==========================================
# GRAPH 2: Buffer Sensitivity (Small vs Opt vs Large)
# ==========================================
def plot_buffer_sensitivity(df, best_buf):
    subset = df[df['Type'] == 'Strong']
    
    # Pick 3 buffers: Smallest, Optimal, Largest
    min_buf = subset['BufferSize'].min()
    max_buf = subset['BufferSize'].max()
    buffers_to_plot = sorted(list(set([min_buf, best_buf, max_buf])))
    
    plt.figure(figsize=(10, 6))
    
    colors = ['red', 'green', 'blue']
    styles = [':', '-', '--']
    
    for i, buf in enumerate(buffers_to_plot):
        data = subset[subset['BufferSize'] == buf]
        label = f'Buffer {buf} (Small)' if buf == min_buf else \
                f'Buffer {buf} (Large)' if buf == max_buf else \
                f'Buffer {buf} (Optimal)'
        
        plt.plot(data['Workers'], data['Efficiency'], 
                 marker='o', linestyle=styles[i], color=colors[i], linewidth=2, label=label)

    plt.axhline(1.0, color='gray', linestyle='-', alpha=0.3)
    plt.axhline(0.7, color='red', linestyle=':', alpha=0.3, label='Threshold 0.7')
    
    plt.title('Does Buffer Size Matter?\nEfficiency Comparison', fontsize=16)
    plt.xlabel('Number of Workers', fontsize=12)
    plt.ylabel('Efficiency', fontsize=12)
    plt.legend()
    plt.ylim(0, 1.1)
    plt.grid(True, alpha=0.5)
    plt.tight_layout()
    plt.savefig('insight_2_buffer_sensitivity.png', dpi=OUTPUT_DPI)
    print("Generated insight_2_buffer_sensitivity.png")

# ==========================================
# GRAPH 3: Classic Log-Log Strong Scaling
# ==========================================
def plot_log_log_scaling(df, best_buf):
    subset = df[(df['Type'] == 'Strong') & (df['BufferSize'] == best_buf)]
    
    plt.figure(figsize=(10, 6))
    
    # Plot Actual Time
    plt.loglog(subset['Workers'], subset['Time'], marker='o', 
               color='royalblue', linewidth=2, markersize=8, label='Actual Time')
    
    # Plot Ideal Time (T1 / N)
    t1 = subset[subset['Workers'] == 1]['Time'].values[0]
    ideal_time = [t1 / w for w in subset['Workers']]
    plt.loglog(subset['Workers'], ideal_time, linestyle='--', color='black', label='Ideal Scaling')
    
    # Formatting
    plt.title(f'Strong Scaling: Time vs Workers (Log-Log)\nBuffer Size {best_buf}', fontsize=16)
    plt.xlabel('Number of Workers (Log Scale)', fontsize=12)
    plt.ylabel('Execution Time (s) (Log Scale)', fontsize=12)
    
    # Fix ticks to show integers 1, 2, 4...
    ax = plt.gca()
    ax.xaxis.set_major_formatter(ticker.ScalarFormatter())
    ax.set_xticks(subset['Workers'])
    
    plt.legend()
    plt.grid(True, which="both", ls="-", alpha=0.2)
    plt.tight_layout()
    plt.savefig('insight_3_log_scaling.png', dpi=OUTPUT_DPI)
    print("Generated insight_3_log_scaling.png")

# ==========================================
# GRAPH 4: Weak Scaling Stability
# ==========================================
def plot_weak_stability(df, best_buf):
    subset = df[(df['Type'] == 'Weak') & (df['BufferSize'] == best_buf)]
    
    plt.figure(figsize=(10, 6))
    
    # Plot Actual Time
    plt.plot(subset['Workers'], subset['Time'], marker='D', color='darkorange', linewidth=2, label='Actual Time')
    
    # Ideal Time (Should be constant T1)
    t1 = subset[subset['Workers'] == 1]['Time'].values[0]
    plt.axhline(t1, color='black', linestyle='--', label=f'Ideal Constant ({t1:.2f}s)')
    
    plt.title(f'Weak Scaling Stability\n(Ideally, this line is flat)', fontsize=16)
    plt.xlabel('Number of Workers (Problem Size Scaled)', fontsize=12)
    plt.ylabel('Execution Time (s)', fontsize=12)
    plt.legend()
    plt.grid(True, alpha=0.5)
    plt.tight_layout()
    plt.savefig('insight_4_weak_stability.png', dpi=OUTPUT_DPI)
    print("Generated insight_4_weak_stability.png")

# ==========================================
# MAIN
# ==========================================
if __name__ == "__main__":
    df = load_and_prep_data(INPUT_FILE)
    
    if df is not None:
        sns.set_theme(style=SNS_THEME)
        
        # 1. Find the global best buffer for Strong scaling to act as our "Champion"
        best_strong_buf = get_optimal_buffer(df, 'Strong')
        print(f"Optimal Buffer Size identified as: {best_strong_buf}")
        
        # 2. Generate Graphs
        plot_efficiency_heatmap(df)
        plot_buffer_sensitivity(df, best_strong_buf)
        plot_log_log_scaling(df, best_strong_buf)
        
        # For Weak scaling, we check if the optimal buffer is different
        best_weak_buf = get_optimal_buffer(df, 'Weak')
        plot_weak_stability(df, best_weak_buf)
        
        print("\nAll insightful plots generated successfully.")