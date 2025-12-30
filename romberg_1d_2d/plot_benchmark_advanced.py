import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns
import numpy as np

# Configuration
INPUT_FILE = 'benchmark_results.csv'
SNS_THEME = "whitegrid"
CMAP = "viridis"

def load_and_prep_data(filename):
    try:
        df = pd.read_csv(filename)
    except FileNotFoundError:
        print(f"Error: {filename} not found. Please run the benchmark first.")
        return None

    # Calculate Speedup and Efficiency
    # We need to find the base time T(1) for each BufferSize to be accurate
    
    # 1. Create a dictionary of T(1) for each Buffer Size in Strong Scaling
    base_times_strong = df[(df['Type'] == 'Strong') & (df['Procs'] == 1)].set_index('BufferSize')['Time'].to_dict()
    
    # 2. Create a dictionary of T(1) for each Buffer Size in Weak Scaling
    base_times_weak = df[(df['Type'] == 'Weak') & (df['Procs'] == 1)].set_index('BufferSize')['Time'].to_dict()

    def calc_metrics(row):
        t1 = 0
        if row['Type'] == 'Strong':
            t1 = base_times_strong.get(row['BufferSize'], row['Time'])
            speedup = t1 / row['Time']
            efficiency = speedup / row['Procs']
        else: # Weak
            t1 = base_times_weak.get(row['BufferSize'], row['Time'])
            # Weak Speedup = (Time(1) * Procs) / Time(N) ? 
            # Standard Weak Efficiency = Time(1) / Time(N)
            efficiency = t1 / row['Time'] 
            speedup = efficiency * row['Procs'] # scaled speedup
            
        return pd.Series([speedup, efficiency])

    df[['Speedup', 'Efficiency']] = df.apply(calc_metrics, axis=1)
    return df

def plot_buffer_sweet_spot(df):
    """
    Plots Time vs Buffer Size for the Maximum Processor count available.
    Helps identify the optimal buffer size.
    """
    max_procs = df['Procs'].max()
    subset = df[(df['Type'] == 'Strong') & (df['Procs'] == max_procs)].sort_values('BufferSize')
    
    if subset.empty:
        print("No data found for Strong scaling at max procs.")
        return

    plt.figure(figsize=(10, 6))
    plt.plot(subset['BufferSize'], subset['Time'], marker='o', linewidth=2, color='b')
    
    plt.title(f'Impact of Buffer Size on Execution Time\n(Strong Scaling, P={max_procs})')
    plt.xlabel('Buffer Size (Elements)')
    plt.ylabel('Execution Time (s)')
    plt.grid(True, which="both", ls="-", alpha=0.5)
    
    # Annotate the minimum point
    min_time = subset['Time'].min()
    best_buf = subset.loc[subset['Time'].idxmin(), 'BufferSize']
    plt.annotate(f'Optimal: {best_buf}', xy=(best_buf, min_time), xytext=(best_buf, min_time*1.1),
                 arrowprops=dict(facecolor='black', shrink=0.05))

    plt.tight_layout()
    plt.savefig('plot_1_buffer_sweet_spot.png')
    print("Generated plot_1_buffer_sweet_spot.png")

def plot_3d_surface(df, scaling_type):
    """
    Creates a 3D surface plot: X=Procs, Y=Buffer, Z=Efficiency
    """
    subset = df[df['Type'] == scaling_type]
    
    if subset.empty:
        return

    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(111, projection='3d')

    # Data for 3D plot must be grid format
    pivot_table = subset.pivot_table(index='BufferSize', columns='Procs', values='Efficiency')
    
    X, Y = np.meshgrid(pivot_table.columns, pivot_table.index)
    Z = pivot_table.values

    surf = ax.plot_surface(X, Y, Z, cmap=CMAP, edgecolor='none', alpha=0.9)
    
    ax.set_title(f'3D Analysis: {scaling_type} Scaling Efficiency\n(Procs vs Buffer Size)')
    ax.set_xlabel('Number of Processors')
    ax.set_ylabel('Buffer Size')
    ax.set_zlabel('Efficiency (0-1)')
    ax.view_init(30, 120) # Adjust angle for better view

    fig.colorbar(surf, shrink=0.5, aspect=5, label='Efficiency')
    
    plt.tight_layout()
    plt.savefig(f'plot_2_3d_{scaling_type.lower()}.png')
    print(f"Generated plot_2_3d_{scaling_type.lower()}.png")

def plot_heatmap(df, scaling_type):
    """
    Creates a Heatmap: X=Procs, Y=Buffer, Color=Efficiency
    Often easier to read than 3D.
    """
    subset = df[df['Type'] == scaling_type]
    if subset.empty: return

    pivot_table = subset.pivot_table(index='BufferSize', columns='Procs', values='Efficiency')

    plt.figure(figsize=(10, 8))
    sns.heatmap(pivot_table, annot=True, fmt=".2f", cmap=CMAP, vmin=0, vmax=1.1)
    
    plt.title(f'{scaling_type} Scaling Efficiency Heatmap')
    plt.xlabel('Number of Processors')
    plt.ylabel('Buffer Size')
    
    plt.tight_layout()
    plt.savefig(f'plot_3_heatmap_{scaling_type.lower()}.png')
    print(f"Generated plot_3_heatmap_{scaling_type.lower()}.png")

if __name__ == "__main__":
    df = load_and_prep_data(INPUT_FILE)
    
    if df is not None:
        sns.set(style=SNS_THEME)
        
        # 1. Sweet Spot Analysis (Line Plot)
        plot_buffer_sweet_spot(df)
        
        # 2. 3D All-in-One Graphs
        plot_3d_surface(df, 'Strong')
        plot_3d_surface(df, 'Weak')
        
        # 3. Heatmaps (Alternative to 3D)
        plot_heatmap(df, 'Strong')
        plot_heatmap(df, 'Weak')
        
        print("\nDone! Check the generated PNG files.")