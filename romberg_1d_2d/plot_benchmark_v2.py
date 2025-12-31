import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# --- CONFIGURATION ---
INPUT_FILE = 'benchmark_results_v2.csv'
SNS_THEME = "whitegrid"

def load_data(filename):
    try:
        df = pd.read_csv(filename)
        return df
    except FileNotFoundError:
        print(f"Error: {filename} not found.")
        return None

def get_best_buffer(df, scaling_type):
    """
    Finds the optimal buffer size based on the run with the MOST workers.
    """
    max_workers = df['Workers'].max()
    # Filter for the specific scaling type at max workers
    subset = df[(df['Type'] == scaling_type) & (df['Workers'] == max_workers)]
    
    if subset.empty:
        return df['BufferSize'].iloc[0] # Fallback
        
    # Find row with minimum time
    best_row = subset.loc[subset['Time'].idxmin()]
    return best_row['BufferSize']

def plot_graphs(df):
    # ==========================================
    # 1. Buffer Size Impact (Sweet Spot)
    # ==========================================
    max_workers = df['Workers'].max()
    subset_buf = df[(df['Type'] == 'Strong') & (df['Workers'] == max_workers)].sort_values('BufferSize')
    
    plt.figure(figsize=(10, 6))
    plt.plot(subset_buf['BufferSize'], subset_buf['Time'], marker='o', color='royalblue', linewidth=2)
    
    # Highlight best
    min_time = subset_buf['Time'].min()
    best_buf_val = subset_buf.loc[subset_buf['Time'].idxmin(), 'BufferSize']
    plt.scatter(best_buf_val, min_time, color='red', s=100, zorder=5, label=f'Optimal: {best_buf_val}')

    plt.title(f'Buffer Size Impact (Workers={max_workers})', fontsize=14)
    plt.xlabel('Buffer Size', fontsize=12)
    plt.ylabel('Execution Time (s)', fontsize=12)
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.tight_layout()
    plt.savefig('opt2_1_buffer_impact.png')
    print(f"Generated opt2_1_buffer_impact.png (Optimal Buffer found: {best_buf_val})")

    # ==========================================
    # Prepare Data for Scaling Plots
    # ==========================================
    # We use the best buffer found above to draw the scaling lines
    best_buf = best_buf_val
    
    # Filter for Strong Scaling
    strong_df = df[(df['Type'] == 'Strong') & (df['BufferSize'] == best_buf)].sort_values('Workers')
    
    # Calculate Metrics
    # Speedup = T(1 Worker) / T(N Workers)
    t1_strong = strong_df[strong_df['Workers'] == 1]['Time'].values[0]
    strong_df['Speedup'] = t1_strong / strong_df['Time']
    strong_df['Efficiency'] = strong_df['Speedup'] / strong_df['Workers']

    # Filter for Weak Scaling
    weak_df = df[(df['Type'] == 'Weak') & (df['BufferSize'] == best_buf)].sort_values('Workers')
    
    # Weak Efficiency = T(1 Worker) / T(N Workers) 
    # (Ideally time stays constant as we add workers + work)
    t1_weak = weak_df[weak_df['Workers'] == 1]['Time'].values[0]
    weak_df['Efficiency'] = t1_weak / weak_df['Time']

    # ==========================================
    # 2. Strong Scaling - Speedup
    # ==========================================
    plt.figure(figsize=(10, 6))
    plt.plot(strong_df['Workers'], strong_df['Speedup'], marker='o', color='green', linewidth=2, label='Actual Speedup')
    # Ideal Line (y=x)
    plt.plot(strong_df['Workers'], strong_df['Workers'], 'k--', label='Ideal Linear', alpha=0.5)
    
    plt.title(f'Strong Scaling Speedup (Buffer {best_buf})', fontsize=14)
    plt.xlabel('Number of Workers', fontsize=12)
    plt.ylabel('Speedup Factor', fontsize=12)
    plt.xticks(strong_df['Workers'])
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.tight_layout()
    plt.savefig('opt2_2_strong_speedup.png')
    print("Generated opt2_2_strong_speedup.png")

    # ==========================================
    # 3. Strong Scaling - Efficiency
    # ==========================================
    plt.figure(figsize=(10, 6))
    plt.plot(strong_df['Workers'], strong_df['Efficiency'], marker='s', color='purple', linewidth=2, label='Actual Efficiency')
    plt.axhline(1.0, color='k', linestyle='--', alpha=0.5, label='Ideal')
    plt.axhline(0.7, color='r', linestyle=':', alpha=0.5, label='Threshold 0.7')
    
    plt.title(f'Strong Scaling Efficiency (Buffer {best_buf})', fontsize=14)
    plt.xlabel('Number of Workers', fontsize=12)
    plt.ylabel('Efficiency (0.0 - 1.0)', fontsize=12)
    plt.ylim(0, 1.1)
    plt.xticks(strong_df['Workers'])
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.tight_layout()
    plt.savefig('opt2_3_strong_efficiency.png')
    print("Generated opt2_3_strong_efficiency.png")

    # ==========================================
    # 4. Weak Scaling - Efficiency
    # ==========================================
    plt.figure(figsize=(10, 6))
    plt.plot(weak_df['Workers'], weak_df['Efficiency'], marker='^', color='orange', linewidth=2, label='Actual Efficiency')
    plt.axhline(1.0, color='k', linestyle='--', alpha=0.5, label='Ideal')
    
    plt.title(f'Weak Scaling Efficiency (Buffer {best_buf})', fontsize=14)
    plt.xlabel('Number of Workers', fontsize=12)
    plt.ylabel('Efficiency (Normalized to 1 Worker)', fontsize=12)
    plt.ylim(0, 1.2)
    plt.xticks(weak_df['Workers'])
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.tight_layout()
    plt.savefig('opt2_4_weak_efficiency.png')
    print("Generated opt2_4_weak_efficiency.png")

if __name__ == "__main__":
    df = load_data(INPUT_FILE)
    if df is not None:
        sns.set(style=SNS_THEME)
        plot_graphs(df)