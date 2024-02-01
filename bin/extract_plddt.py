import argparse
import os
import pandas as pd

from Bio import PDB
import statistics
import numpy as np
import logging

from concurrent.futures import ProcessPoolExecutor

# Dictionary to map 3-letter amino acid codes to 1-letter codes
AMINO_ACID_DICT = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D',
    'CYS': 'C', 'GLN': 'Q', 'GLU': 'E', 'GLY': 'G',
    'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
    'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
    'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
}

def setup_logging():
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def compute_statistics(values):
    logging.info('Computing statistics...')
    stats = {
        'sequence_length': len(values),
        'plddt_mean': statistics.mean(values),
        'plddt_median': statistics.median(values),
        'plddt_stdev': statistics.stdev(values),
        'plddt_variance': statistics.variance(values),
        'plddt_pstdev': statistics.pstdev(values),
        'plddt_pvariance': statistics.pvariance(values),
        'plddt_min': min(values),
        'plddt_max': max(values),
        'plddt_q1': np.percentile(values, 25),
        'plddt_q2': np.percentile(values, 50),
        'plddt_q3': np.percentile(values, 75),
    }

    try:
        stats['plddt_mode'] = statistics.mode(values)
    except statistics.StatisticsError:
        stats['plddt_mode'] = float('na')

    logging.info('Statistics computed.')
    return stats

def extract_b_factors(pdb_file):
    logging.info(f'Processing {pdb_file}...')
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('protein', pdb_file)
    data = []
    index = 0  # Initialize index
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    if atom.get_name() == 'CA':  # Only considering alpha-carbon atoms
                        index += 1  # Increment index
                        amino_acid = residue.get_resname()
                        b_factor = atom.get_bfactor()
                        data.append((index, amino_acid, b_factor))
    logging.info(f'Finished processing {pdb_file}.')
    return pdb_file, data  # Return pdb_file as well to identify which file the data belongs to

def write_csv(pdb_file, data, output_dir):
    df = pd.DataFrame(data, columns=['index', 'amino_acid', 'plddt'])  # Adjusted columns
    df['residue'] = df['amino_acid'].map(AMINO_ACID_DICT)
    df = df.reindex(columns=['index', 'residue', 'amino_acid', 'plddt'])

    filename_with_suffix = os.path.basename(pdb_file)
    filename_without_suffix = filename_with_suffix.replace('_model_1_ptm_relaxed.pdb', '')
    csv_file_name = f'{output_dir}/prot_{filename_without_suffix}.csv'
    df.to_csv(csv_file_name, index=False, header=True, sep=";", decimal=".")
    logging.info(f'Data CSV file created: {csv_file_name}')

def generate_csvs(directory, statistics_name, output_dir):
    logging.info(f'Scanning directory {directory} for PDB files...')
    pdb_files = [f for f in os.listdir(directory) if f.endswith('.pdb')]
    logging.info(f'PDB files found: {len(pdb_files)}')

    all_data = []
    stats_data = []

    with ProcessPoolExecutor() as executor:
        results = list(executor.map(extract_b_factors, [os.path.join(directory, pdb_file) for pdb_file in pdb_files]))

    for result in results:
        pdb_file, data = result
        all_data.extend(data)
        write_csv(pdb_file, data, output_dir)

        values = [item[2] for item in data]  # Extracting plddt values from data
        stats = compute_statistics(values)

        pdb_file_protein = os.path.basename(pdb_file).replace('.pdb', '')
        stats_data.append((pdb_file_protein, stats['sequence_length'], stats['plddt_mean'], stats['plddt_median'], stats['plddt_stdev'], stats['plddt_variance'],
                           stats['plddt_pstdev'], stats['plddt_pvariance'], stats['plddt_min'], stats['plddt_max'],
                           stats['plddt_q1'], stats['plddt_q2'], stats['plddt_q3'], stats['plddt_mode']))

    stats_df = pd.DataFrame(stats_data, columns=['prot_name', 'sequence_length', 'plddt_mean', 'plddt_median', 'plddt_stdev', 'plddt_variance',
                                                 'plddt_pstdev', 'plddt_pvariance', 'plddt_min', 'plddt_max',
                                                 'plddt_q1', 'plddt_q2', 'plddt_q3', 'plddt_mode'])

    stats_df.to_csv(statistics_name, index=False, header=True, sep=";", decimal=".")
    logging.info('Statistics CSV file created: statistics.csv')

def main():
    parser = argparse.ArgumentParser(description='Process PDB files and extract B-Factor values.')
    parser.add_argument('directory', type=str, help='Path to the directory containing PDB files.')
    parser.add_argument('statistics_name', type=str, help='')
    parser.add_argument('output_dir', type=str, help='')
    args = parser.parse_args()

    setup_logging()
    logging.info('Script initiated.')
    generate_csvs(args.directory, args.statistics_name, args.output_dir)
    logging.info('Script completed successfully.')

if __name__ == '__main__':
    main()