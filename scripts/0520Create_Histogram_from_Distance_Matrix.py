#!/opt/conda/bin/python

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from typing import Tuple


def read_distance_matrix(file_path: str) -> np.ndarray:
    """
    Read distance matrix from CSV file

    Args:
        file_path: Path to input CSV file

    Returns:
        np.ndarray: Array of distances against reference
    """
    try:
        # Read first line to get sample names
        with open(file_path, 'r') as f:
            sample_names = f.readline().strip().split(',')

        # Read the second line (distances against reference)
        distances = pd.read_csv(file_path, skiprows=0, nrows=1)

        # Convert to numpy array and skip the first value (self-distance)
        return distances.iloc[0, 1:].values

    except Exception as e:
        print(f"Error reading file: {str(e)}")
        raise


def calculate_fd_bins(data: np.ndarray) -> Tuple[int, float]:
    """
    Calculate number of bins using Freedman-Diaconis rule

    Args:
        data: Input data array

    Returns:
        Tuple[int, float]: Number of bins and bin width
    """
    # Calculate IQR
    q75, q25 = np.percentile(data, [75, 25])
    iqr = q75 - q25

    # Calculate bin width using F-D rule
    n = len(data)
    bin_width = 2 * iqr * (n ** (-1/3))

    # Calculate number of bins
    data_range = np.max(data) - np.min(data)
    n_bins = int(np.ceil(data_range / bin_width))

    return n_bins, bin_width


def create_histogram(
        distances: np.ndarray,
        output_dir: str
) -> None:
    """
    Create histograms using both automatic and F-D binning

    Args:
        distances: Array of distance values
        output_dir: Output directory path
    """
    # Calculate F-D bins
    n_bins, bin_width = calculate_fd_bins(distances)

    # Create figure with two subplots
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))

    # Plot using automatic binning
    ax1.hist(distances, bins='auto', edgecolor='black', alpha=0.7)
    ax1.set_title('Automatic Binning')
    ax1.set_xlabel('Distance from Reference')
    ax1.set_ylabel('Frequency')
    ax1.grid(True, alpha=0.3)

    # Plot using F-D rule
    ax2.hist(distances, bins=n_bins, edgecolor='black', alpha=0.7)
    ax2.set_title(f'Freedman-Diaconis Rule\n(n_bins={n_bins}, width={bin_width:.2f})')
    ax2.set_xlabel('Distance from Reference')
    ax2.set_ylabel('Frequency')
    ax2.grid(True, alpha=0.3)

    # Adjust layout
    plt.tight_layout()

    # Print binning information
    print(f"Freedman-Diaconis rule:")
    print(f"Number of bins: {n_bins}")
    print(f"Bin width: {bin_width:.2f}")

    # Save plots
    plt.savefig(
        os.path.join(output_dir, "distance_histogram_comparison.png"),
        dpi=300,
        bbox_inches='tight'
    )
    plt.close()


def main():
    # Set input/output paths
    input_dir = "./0510_distance_matrix"
    output_dir = input_dir
    input_file = os.path.join(input_dir, "haplotype_distance_matrix.csv")

    try:
        # Read distance matrix
        distances = read_distance_matrix(input_file)

        # Create histograms
        create_histogram(distances, output_dir)

        print(f"Successfully created histogram comparison")

    except Exception as e:
        print(f"An error occurred: {str(e)}")
        raise


if __name__ == "__main__":
    main()