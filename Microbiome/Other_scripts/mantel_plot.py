import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

def plot_distance_matrices(dm_geodesic_path, dm_other_path, output_plot_path=None, 
                           x_label='Geodesic Distance', y_label='Diversity Distance'):
    """
    Generate a scatter plot between two distance matrices.

    Parameters:
    dm_geodesic_path (str): Path to the geodesic distance matrix (TSV file).
    dm_other_path (str): Path to the other distance matrix (TSV file).
    output_plot_path (str, optional): Path to save the resulting plot. If None, the plot is not saved.
    x_label (str): Label for the x-axis of the scatter plot.
    y_label (str): Label for the y-axis of the scatter plot.
    """
    # Load distance matrices
    dm1 = pd.read_csv(dm_geodesic_path, sep='\t', index_col=0)
    dm2 = pd.read_csv(dm_other_path, sep='\t', index_col=0)

    # Align matrices by common indices
    common_ids = dm1.index.intersection(dm2.index)
    if len(common_ids) == 0:
        raise ValueError("No common IDs found between the two matrices.")

    dm1 = dm1.loc[common_ids, common_ids]
    dm2 = dm2.loc[common_ids, common_ids]

    # Extract upper triangular values (excluding diagonal)
    dm1_flat = dm1.where(np.triu(np.ones(dm1.shape), k=1).astype(bool)).stack()
    dm2_flat = dm2.where(np.triu(np.ones(dm2.shape), k=1).astype(bool)).stack()

    # Create scatter plot with regression line
    plt.figure(figsize=(6, 4))
    sns.regplot(x=dm1_flat, y=dm2_flat, scatter_kws={'s': 10, 'color': '#56C667FF'}, line_kws={'color': '#3F4788FF'})
    plt.xlabel(x_label)
    plt.ylabel(y_label)

    # Show the plot with grid
    plt.grid(True)
    plt.tight_layout()

    # Save plot if output path is specified
    if output_plot_path:
        plt.savefig(output_plot_path, dpi=300)
    plt.show()

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Generate scatter plot between two distance matrices.")
    parser.add_argument("dm_geodesic_path", help="Path to the geodesic distance matrix (TSV file).")
    parser.add_argument("dm_other_path", help="Path to the other distance matrix (TSV file).")
    parser.add_argument("--output_plot_path", help="Path to save the resulting plot.", default=None)
    parser.add_argument("--x_label", help="Label for the x-axis of the scatter plot.", default="Geodesic Distance")
    parser.add_argument("--y_label", help="Label for the y-axis of the scatter plot.", default="Diversity Distance")

    args = parser.parse_args()

    plot_distance_matrices(
        dm_geodesic_path=args.dm_geodesic_path,
        dm_other_path=args.dm_other_path,
        output_plot_path=args.output_plot_path,
        x_label=args.x_label,
        y_label=args.y_label
    )
