import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# --- CONFIGURATION ---
INPUT_FILE = 'benchmark_results.csv'  
SNS_THEME = "whitegrid"

def load_data(filename):
    try:
        df = pd.read_csv(filename)
        return df
    except FileNotFoundError:
        print(f"Error: '{filename}' not found. Please rename your file or update INPUT_FILE in the script.")
        return None

def get_best_buffer(df, scaling_type):
    """
    Finds the 'Optimal Buffer Size' by looking at the run with the 
    Maximum Processors and finding which Buffer Size gave the minimum time.
    """
    max_procs = df['Procs'].max()
    subset = df[(df['Type'] == scaling_type) & (df['Procs'] == max_procs)]
    
    best_row = subset.loc[subset['Time'].idxmin()]
    return best_row['BufferSize']

def plot_buffer_impact(df):
    """
    Plot 1: Execution Time vs Buffer Size (at Max Procs).
    Shows the 'Sweet Spot'.
    """
    max_procs = df['Procs'].max()
    subset = df[(df['Type'] == 'Strong') & (df['Procs'] == max_procs)].sort_values('BufferSize')
    
    plt.figure(figsize=(10, 6))
    plt.plot(subset['BufferSize'], subset['Time'], marker='o', linewidth=2, color='royalblue')
    
    min_time = subset['Time'].min()
    best_buf = subset.loc[subset['Time'].idxmin(), 'BufferSize']
    
    plt.scatter(best_buf, min_time, color='red', s=100, zorder=5, label=f'Optimal: {best_buf}')
    
    plt.title(f'Impact of Buffer Size on Execution Time\n(Strong Scaling, P={max_procs})', fontsize=14)
    plt.xlabel('Buffer Size', fontsize=12)
    plt.ylabel('Time (seconds) - Lower is Better', fontsize=12)
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.7)
    
    plt.tight_layout()
    plt.savefig('simple_1_buffer_impact.png')
    print(f"Generated 'simple_1_buffer_impact.png' (Optimal Buffer: {best_buf})")

def plot_strong_scaling(df, best_buffer):
    """
    Plot 2 & 3: Speedup and Efficiency vs Processors 
    (Using ONLY the data from the Optimal Buffer Size).
    """
    subset = df[(df['Type'] == 'Strong') & (df['BufferSize'] == best_buffer)].sort_values('Procs')
    
    t1 = subset[subset['Procs'] == 1]['Time'].values[0]
    subset['Speedup'] = t1 / subset['Time']
    subset['Efficiency'] = subset['Speedup'] / subset['Procs']

    plt.figure(figsize=(10, 6))
    plt.plot(subset['Procs'], subset['Speedup'], marker='o', linewidth=2, color='green', label=f'Actual (Buffer {best_buffer})')
    plt.plot(subset['Procs'], subset['Procs'], 'k--', label='Ideal Linear', alpha=0.5)
    
    plt.title(f'Strong Scaling Speedup (Buffer Size {best_buffer})', fontsize=14)
    plt.xlabel('Number of Processors', fontsize=12)
    plt.ylabel('Speedup', fontsize=12)
    plt.xticks(subset['Procs']) 
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.tight_layout()
    plt.savefig('simple_2_strong_speedup.png')
    print("Generated 'simple_2_strong_speedup.png'")

    # --- Graph 3: Efficiency ---
    plt.figure(figsize=(10, 6))
    plt.plot(subset['Procs'], subset['Efficiency'], marker='s', linewidth=2, color='purple', label=f'Actual (Buffer {best_buffer})')
    plt.axhline(y=1.0, color='k', linestyle='--', alpha=0.5, label='Ideal')
    plt.axhline(y=0.7, color='r', linestyle=':', alpha=0.5, label='Good Threshold (0.7)')
    
    plt.title(f'Strong Scaling Efficiency (Buffer Size {best_buffer})', fontsize=14)
    plt.xlabel('Number of Processors', fontsize=12)
    plt.ylabel('Efficiency (0.0 - 1.0)', fontsize=12)
    plt.ylim(0, 1.1)
    plt.xticks(subset['Procs'])
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.tight_layout()
    plt.savefig('simple_3_strong_efficiency.png')
    print("Generated 'simple_3_strong_efficiency.png'")

def plot_weak_scaling(df, best_buffer):
    """
    Plot 4: Weak Scaling Efficiency vs Processors
    """
    subset = df[(df['Type'] == 'Weak') & (df['BufferSize'] == best_buffer)].sort_values('Procs')
    
    
    t1 = subset[subset['Procs'] == 1]['Time'].values[0]
    subset['Efficiency'] = t1 / subset['Time']

    plt.figure(figsize=(10, 6))
    plt.plot(subset['Procs'], subset['Efficiency'], marker='^', linewidth=2, color='orange', label=f'Actual (Buffer {best_buffer})')
    plt.axhline(y=1.0, color='k', linestyle='--', alpha=0.5, label='Ideal')
    
    plt.title(f'Weak Scaling Efficiency (Buffer Size {best_buffer})', fontsize=14)
    plt.xlabel('Number of Processors', fontsize=12)
    plt.ylabel('Efficiency (Normalized to P=1)', fontsize=12)
    plt.ylim(0, 1.2)
    plt.xticks(subset['Procs'])
    plt.legend()
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.tight_layout()
    plt.savefig('simple_4_weak_efficiency.png')
    print("Generated 'simple_4_weak_efficiency.png'")

if __name__ == "__main__":
    df = load_data(INPUT_FILE)
    
    if df is not None:
        sns.set(style=SNS_THEME)
        
        plot_buffer_impact(df)
        
        best_strong_buffer = get_best_buffer(df, 'Strong')
        plot_strong_scaling(df, best_strong_buffer)
        
        best_weak_buffer = get_best_buffer(df, 'Weak')
        plot_weak_scaling(df, best_weak_buffer)
        
        print("\nAll plots generated successfully!")