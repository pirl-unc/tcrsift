# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from os import PathLike
from os.path import join 
import sys

import argh
import scanpy as sc
import pandas as pd


def load_sample(
        sample_name : str,
        gene_expression__dir : str | PathLike,
        vdj_dir : str | PathLike):
    gene_expression_matrix_dir = join(gene_expression__dir, 'filtered_feature_bc_matrix')
    vdj_clonotypes_csv_path = join(vdj_dir, 'clonotypes.csv')
    vdj_annotations_csv_path = join(vdj_dir, 'all_annotations.csv')
    gene_expression_data = sc.read_10x_mtx(
        gene_expression_matrix_dir, 
        var_names='gene_ids')
    df_clonotypes = pd.read_csv(vdj_clonotypes_csv_path)
    df_annotations = pd.read_csv(vdj_annotations_csv_path)
    print("Loaded data")

def run_from_sample_sheet_for_count_and_vdj_outputs(sample_sheet_csv_path : str | PathLike):
    """
    Load CellRanger 5'GEX+3'VDJ data for each sample based on paths in a sample sheet CSV file.
    Assumes you have run 'cellranger vdj' and 'cellranger count' on the samples.
    
    Parameters
    ----------
    sample_sheet_csv_path
        CSV file with columns 'sample', 'gene_expression_dir', 'vdj_dir'.
    """
    print("Running from sample sheet")
    sample_sheet = pd.read_csv(sample_sheet_csv_path)
    for _, row in sample_sheet.iterrows():
        print(row)
        load_sample(
            sample_name=row.sample,
            row['gene_expression_matrix_dir'], 
            row['vdj_clonotypes_csv_path'], 
            row['vdj_annotations_csv_path'])

def main(args=None):
    if args is None:
        args = sys.argv[1:]
    parser = argh.ArghParser()
    
    argh.set_default_command(parser, run_from_sample_sheet_for_count_and_vdj_outputs)
    argh.dispatch(parser, argv=args)