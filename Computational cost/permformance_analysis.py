#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

# Set plot style
plt.style.use('ggplot')
sns.set_context("talk")

# Load the data
df = pd.read_csv('timing.csv')

# Ensure data types are correct
df['n_elements'] = pd.to_numeric(df['n_elements'])
df['seconds'] = pd.to_numeric(df['seconds'])

# Create a figure with multiple subplots
fig, axes = plt.subplots(2, 2, figsize=(16, 12))
fig.suptitle('Burgers Equation Solver Performance Analysis', fontsize=20)

# 1. Runtime vs Problem Size (separate by precision)
ax1 = axes[0, 0]
for precision in df['precision'].unique():
    for version in df['version'].unique():
        subset = df[(df['precision'] == precision) & (df['version'] == version)]
        label = f"{precision}, {version}"
        ax1.plot(subset['n_elements'], subset['seconds'], 
                marker='o', linewidth=2, label=label)
        
        # Add power law fit
        if len(subset) > 2:
            x = subset['n_elements']
            y = subset['seconds']
            log_x = np.log(x)
            log_y = np.log(y)
            slope, intercept = np.polyfit(log_x, log_y, 1)
            fit_y = np.exp(intercept) * x**slope
            ax1.plot(x, fit_y, '--', alpha=0.7)

ax1.set_xscale('log', base=2)
ax1.set_yscale('log')
ax1.set_xlabel('Number of Elements')
ax1.set_ylabel('Runtime (seconds)')
ax1.set_title('Runtime vs Problem Size')
ax1.grid(True, which="both", ls="-", alpha=0.2)
ax1.legend(title='Config', fontsize=8)

# 2. Computational Efficiency (n_elements / runtime)
ax2 = axes[0, 1]
df['efficiency'] = df['n_elements'] / df['seconds']

for precision in df['precision'].unique():
    for version in df['version'].unique():
        subset = df[(df['precision'] == precision) & (df['version'] == version)]
        label = f"{precision}, {version}"
        ax2.plot(subset['n_elements'], subset['efficiency'], 
                marker='o', linewidth=2, label=label)

ax2.set_xscale('log', base=2)
ax2.set_yscale('log')
ax2.set_xlabel('Number of Elements')
ax2.set_ylabel('Elements/Second')
ax2.set_title('Computational Efficiency')
ax2.grid(True, which="both", ls="-", alpha=0.2)
ax2.legend(title='Config', fontsize=8)

# 3. Scaling Analysis - table of exponents
scaling_data = []
for precision in df['precision'].unique():
    for version in df['version'].unique():
        subset = df[(df['precision'] == precision) & (df['version'] == version)]
        if len(subset) > 2:
            log_x = np.log(subset['n_elements'])
            log_y = np.log(subset['seconds'])
            slope, intercept = np.polyfit(log_x, log_y, 1)
            scaling_data.append({
                'precision': precision,
                'version': version,
                'scaling_exponent': slope,
                'coefficient': np.exp(intercept)
            })

scaling_df = pd.DataFrame(scaling_data)
print("\nScaling Analysis (Runtime ~ N^exponent):")
print(scaling_df)

# Create bar plot for scaling exponents
ax3 = axes[1, 0]
sns.barplot(x='precision', y='scaling_exponent', hue='version', data=scaling_df, ax=ax3)
ax3.set_title('Scaling Exponents')
ax3.set_ylabel('Exponent (Runtime ~ N^exponent)')
ax3.axhline(y=1, color='r', linestyle='--', alpha=0.7)  # Linear scaling reference
ax3.axhline(y=2, color='g', linestyle='--', alpha=0.7)  # Quadratic scaling reference
ax3.grid(axis='y', alpha=0.3)

# 4. Performance improvement of KS over plain
improvement_data = []
for precision in df['precision'].unique():
    for n_elem in df['n_elements'].unique():
        plain_data = df[(df['precision'] == precision) & 
                      (df['version'] == 'plain') & 
                      (df['n_elements'] == n_elem)]
        
        ks_data = df[(df['precision'] == precision) & 
                   (df['version'] == 'ks') & 
                   (df['n_elements'] == n_elem)]
        
        if not plain_data.empty and not ks_data.empty:
            plain_time = plain_data.iloc[0]['seconds']
            ks_time = ks_data.iloc[0]['seconds']
            
            improvement = (plain_time - ks_time) / plain_time * 100  # Percentage
            
            improvement_data.append({
                'precision': precision,
                'n_elements': n_elem,
                'improvement': improvement
            })

if improvement_data:
    imp_df = pd.DataFrame(improvement_data)
    
    ax4 = axes[1, 1]
    sns.barplot(x='precision', y='improvement', hue='n_elements', data=imp_df, ax=ax4)
    ax4.set_title('KS vs Plain Improvement')
    ax4.set_ylabel('Improvement (%)')
    ax4.axhline(y=0, color='r', linestyle='--', alpha=0.7)
    ax4.grid(axis='y', alpha=0.3)
    
    # Print average improvement by precision
    avg_imp = imp_df.groupby('precision')['improvement'].mean()
    print("\nAverage Improvement (KS vs Plain):")
    print(avg_imp)

plt.tight_layout()
plt.savefig('burgers_performance_analysis.png', dpi=300)
print("Analysis complete. Figure saved as 'burgers_performance_analysis.png'")

# Create a more detailed view by precision
plt.figure(figsize=(20, 15))

for i, precision in enumerate(df['precision'].unique()):
    plt.subplot(2, 2, i+1)
    
    for version in df['version'].unique():
        subset = df[(df['precision'] == precision) & (df['version'] == version)]
        
        plt.plot(subset['n_elements'], subset['seconds'], 
                marker='o', linewidth=2, markersize=10, label=f"{version}")
        
        # Add power law fit
        if len(subset) > 2:
            x = subset['n_elements']
            y = subset['seconds']
            log_x = np.log(x)
            log_y = np.log(y)
            slope, intercept = np.polyfit(log_x, log_y, 1)
            fit_y = np.exp(intercept) * x**slope
            plt.plot(x, fit_y, '--', alpha=0.7)
            
            # Add equation
            plt.text(x.iloc[-1], y.iloc[-1]*1.1, 
                    f"T ~ N^{slope:.2f}", 
                    ha='right', fontsize=12)
    
    plt.xscale('log', base=2)
    plt.yscale('log')
    plt.xlabel('Number of Elements')
    plt.ylabel('Runtime (seconds)')
    plt.title(f'Runtime vs Problem Size - {precision.upper()}')
    plt.grid(True, which="both", ls="-", alpha=0.2)
    plt.legend(title='Version')

plt.tight_layout()
plt.savefig('precision_comparison.png', dpi=300)
print("Precision comparison saved as 'precision_comparison.png'")

# Calculate theoretical memory usage
memory_data = []
for precision in df['precision'].unique():
    for n_elem in df['n_elements'].unique():
        # Calculate bytes per value based on precision
        if precision == 'fp16':
            bytes_per_val = 2
        elif precision == 'fp32':
            bytes_per_val = 4
        elif precision == 'fp64':
            bytes_per_val = 8
        else:  # fp128
            bytes_per_val = 16
        
        # Calculate memory usage
        n_nodes = n_elem + 1
        # Main arrays: x, u, u_old
        vector_mem = 3 * n_nodes * bytes_per_val
        # Matrices: M, K, C, A
        matrix_mem = 4 * n_nodes * n_nodes * bytes_per_val
        total_mem = (vector_mem + matrix_mem) / (1024 * 1024)  # Convert to MB
        
        for version in df['version'].unique():
            subset = df[(df['precision'] == precision) & 
                      (df['version'] == version) & 
                      (df['n_elements'] == n_elem)]
            
            if not subset.empty:
                runtime = subset.iloc[0]['seconds']
                
                memory_data.append({
                    'precision': precision,
                    'version': version,
                    'n_elements': n_elem,
                    'memory_mb': total_mem,
                    'runtime': runtime
                })

memory_df = pd.DataFrame(memory_data)

# Plot memory usage vs runtime
plt.figure(figsize=(12, 8))

for precision in memory_df['precision'].unique():
    for version in memory_df['version'].unique():
        subset = memory_df[(memory_df['precision'] == precision) & 
                         (memory_df['version'] == version)]
        
        plt.plot(subset['memory_mb'], subset['runtime'], 
                marker='o', linewidth=2, markersize=8, 
                label=f"{precision}, {version}")

plt.xscale('log')
plt.yscale('log')
plt.xlabel('Estimated Memory Usage (MB)')
plt.ylabel('Runtime (seconds)')
plt.title('Runtime vs Memory Usage')
plt.grid(True, which="both", ls="-", alpha=0.2)
plt.legend(title='Configuration')
plt.tight_layout()
plt.savefig('memory_vs_runtime.png', dpi=300)
print("Memory analysis saved as 'memory_vs_runtime.png'")

# Create summary table
summary = pd.DataFrame({
    'Configuration': [
        f"{df.loc[df['seconds'].idxmin(), 'precision']} {df.loc[df['seconds'].idxmin(), 'version']} (n={df.loc[df['seconds'].idxmin(), 'n_elements']})",
        f"{df.loc[df['seconds'].idxmax(), 'precision']} {df.loc[df['seconds'].idxmax(), 'version']} (n={df.loc[df['seconds'].idxmax(), 'n_elements']})",
        f"{scaling_df.loc[scaling_df['scaling_exponent'].idxmin(), 'precision']} {scaling_df.loc[scaling_df['scaling_exponent'].idxmin(), 'version']}",
        f"{scaling_df.loc[scaling_df['scaling_exponent'].idxmax(), 'precision']} {scaling_df.loc[scaling_df['scaling_exponent'].idxmax(), 'version']}"
    ],
    'Metric': [
        'Fastest configuration',
        'Slowest configuration',
        'Best scaling (lowest exponent)',
        'Worst scaling (highest exponent)'
    ],
    'Value': [
        f"{df['seconds'].min():.6f} s",
        f"{df['seconds'].max():.6f} s",
        f"O(n^{scaling_df['scaling_exponent'].min():.2f})",
        f"O(n^{scaling_df['scaling_exponent'].max():.2f})"
    ]
})

print("\nPerformance Summary:")
print(summary)

# Save summary to CSV
summary.to_csv('performance_summary.csv', index=False)
print("Summary saved to 'performance_summary.csv'")