#!/usr/bin/env python3

"""
Plot KMC k-mer frequency histogram
"""
import argparse
import matplotlib.pyplot as plt
import pandas as pd
import sys


def plot_kmc_histogram(histogram_file, output_file=None, max_freq=None, min_freq=None, log_scale=False):
    """
    Plot k-mer frequency histogram from KMC output
    
    Args:
        histogram_file: Path to histogram file from kmc_tools
        output_file: Path to save plot (if None, display interactively)
        max_freq: Maximum frequency to plot (None = plot all)
        min_freq: Minimum frequency to plot (None = plot all)
        log_scale: Use log scale for y-axis
    """
    # Read histogram data
    df = pd.read_csv(histogram_file, sep='\t', header=None, names=['frequency', 'count'])
    
    # Filter by max frequency if specified
    if max_freq:
        df = df[df['frequency'] <= max_freq]

    if min_freq:
        df = df[df['frequency'] >= min_freq]

    # Create plot
    fig, ax = plt.subplots(figsize=(12, 6))
    
    ax.plot(df['frequency'], df['count'], linewidth=1.5)
    ax.fill_between(df['frequency'], df['count'], alpha=0.3)
    
    ax.set_xlabel('K-mer Frequency', fontsize=12)
    ax.set_ylabel('Count', fontsize=12)
    ax.set_title('K-mer Frequency Distribution', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3)
    
    if log_scale:
        ax.set_yscale('log')
    
    # Add statistics
    total_kmers = df['count'].sum()
    peak_freq = df.loc[df['count'].idxmax(), 'frequency']
    
    stats_text = f'Total k-mers: {total_kmers:,}\nPeak frequency: {peak_freq}'
    ax.text(0.98, 0.98, stats_text, 
            transform=ax.transAxes, 
            fontsize=10,
            verticalalignment='top',
            horizontalalignment='right',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    plt.tight_layout()
    
    # Save or show
    if output_file:
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Plot saved to {output_file}")
    else:
        plt.show()
    
    # Print statistics
    print(f"\nK-mer Statistics:")
    print(f"  Total unique k-mers: {len(df):,}")
    print(f"  Total k-mer count: {total_kmers:,}")
    print(f"  Peak frequency: {peak_freq}")
    print(f"  Mean frequency: {(df['frequency'] * df['count']).sum() / total_kmers:.2f}")


def main():
    parser = argparse.ArgumentParser(
        description='Plot KMC k-mer frequency histogram',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Generate histogram from KMC database first:
  kmc_tools transform kmc_test histogram kmc_histogram.txt -cx1000
  
  # Plot the histogram:
  python plot_kmc_histogram.py kmc_histogram.txt -o kmer_plot.png
  
  # Plot with log scale and max frequency:
  python plot_kmc_histogram.py kmc_histogram.txt -o kmer_plot.png --log --max-freq 200
        """
    )
    
    parser.add_argument('histogram', help='KMC histogram file from kmc_tools')
    parser.add_argument('-o', '--output', help='Output plot file (PNG, PDF, SVG)')
    parser.add_argument('--max-freq', type=int, help='Maximum frequency to plot')
    parser.add_argument('--min-freq', type=int, help='Minimum frequency to plot')
    parser.add_argument('--log', action='store_true', help='Use log scale for y-axis')
    
    args = parser.parse_args()
    
    try:
        plot_kmc_histogram(args.histogram, args.output, args.max_freq, args.min_freq, args.log)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
