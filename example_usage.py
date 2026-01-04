"""
Example usage of the professional scRNA perturbation analysis pipeline.

This script demonstrates how to customize and run the analysis pipeline.
"""

from scrna_perturb_professional import (
    Config, load_data, calculate_qc_metrics, filter_data,
    preprocess_data, find_hvg_and_reduce_dimensions,
    prepare_perturbation_labels, plot_perturbation_umap,
    perform_differential_expression, plot_de_results,
    generate_summary_report, setup_logging
)
from pathlib import Path
import logging

# Example 1: Run with default configuration
def example_default():
    """Run analysis with default settings."""
    logger = setup_logging(Config.OUTPUT_DIR)
    
    # Load and process data
    adata = load_data(Config.DATA_FILE, Config.N_CELLS_TO_SUBSAMPLE, logger)
    adata = calculate_qc_metrics(adata, logger)
    adata = filter_data(adata, logger)
    adata = preprocess_data(adata, logger)
    adata = find_hvg_and_reduce_dimensions(adata, logger)
    adata = prepare_perturbation_labels(adata, logger)
    
    # Generate outputs
    plot_perturbation_umap(adata, Config.OUTPUT_DIR, logger)
    de_results = perform_differential_expression(adata, Config.OUTPUT_DIR, logger)
    plot_de_results(adata, de_results, Config.OUTPUT_DIR, logger)
    generate_summary_report(adata, de_results, Config.OUTPUT_DIR, logger)
    
    # Save processed data
    output_file = Path(Config.OUTPUT_DIR) / 'processed_data.h5ad'
    adata.write(output_file)
    logger.info(f"Analysis complete! Results saved to {Config.OUTPUT_DIR}")


# Example 2: Customize configuration for different perturbation
def example_custom_perturbation():
    """Run analysis with custom perturbation target."""
    # Modify config
    Config.TARGET_PERTURBATION = 'BRCA1'  # Change target gene
    Config.OUTPUT_DIR = 'perturbation_analysis_BRCA1'
    
    logger = setup_logging(Config.OUTPUT_DIR)
    logger.info(f"Analyzing {Config.TARGET_PERTURBATION} perturbation")
    
    # Run analysis (same as above)
    adata = load_data(Config.DATA_FILE, Config.N_CELLS_TO_SUBSAMPLE, logger)
    adata = calculate_qc_metrics(adata, logger)
    adata = filter_data(adata, logger)
    adata = preprocess_data(adata, logger)
    adata = find_hvg_and_reduce_dimensions(adata, logger)
    adata = prepare_perturbation_labels(adata, logger)
    
    plot_perturbation_umap(adata, Config.OUTPUT_DIR, logger)
    de_results = perform_differential_expression(adata, Config.OUTPUT_DIR, logger)
    plot_de_results(adata, de_results, Config.OUTPUT_DIR, logger)
    generate_summary_report(adata, de_results, Config.OUTPUT_DIR, logger)
    
    output_file = Path(Config.OUTPUT_DIR) / 'processed_data.h5ad'
    adata.write(output_file)


# Example 3: Analyze multiple perturbations
def example_multiple_perturbations():
    """Run analysis for multiple perturbation targets."""
    perturbations = ['TP53', 'BRCA1', 'RB1']  # List of genes to analyze
    
    for pert in perturbations:
        logger = logging.getLogger(__name__)
        logger.info(f"\n{'='*60}")
        logger.info(f"Analyzing {pert} perturbation")
        logger.info(f"{'='*60}")
        
        # Update config for this perturbation
        Config.TARGET_PERTURBATION = pert
        Config.OUTPUT_DIR = f'perturbation_analysis_{pert}'
        
        # Run analysis
        adata = load_data(Config.DATA_FILE, Config.N_CELLS_TO_SUBSAMPLE, logger)
        adata = calculate_qc_metrics(adata, logger)
        adata = filter_data(adata, logger)
        adata = preprocess_data(adata, logger)
        adata = find_hvg_and_reduce_dimensions(adata, logger)
        adata = prepare_perturbation_labels(adata, logger)
        
        plot_perturbation_umap(adata, Config.OUTPUT_DIR, logger)
        de_results = perform_differential_expression(adata, Config.OUTPUT_DIR, logger)
        plot_de_results(adata, de_results, Config.OUTPUT_DIR, logger)
        generate_summary_report(adata, de_results, Config.OUTPUT_DIR, logger)
        
        output_file = Path(Config.OUTPUT_DIR) / 'processed_data.h5ad'
        adata.write(output_file)


if __name__ == '__main__':
    # Run default example
    print("Running example with default configuration...")
    example_default()
    
    # Uncomment to run other examples:
    # example_custom_perturbation()
    # example_multiple_perturbations()

