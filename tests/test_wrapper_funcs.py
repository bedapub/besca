import pathlib
import pytest
from typing import List
from os.path import join
import filecmp
from random import sample
import anndata
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix

from scvelo import AnnData

import besca as bc
import scanpy as sc
from besca.st._wrapper_funcs import additional_labeling, highly_variable_genes
from besca.Import._read import all_symbols_ok

pytest.raw_data_subset_size = 1000


def assert_gct_equal(reference_file: str, output_file: str) -> None:
    """Compare two .gct files allowing for rounding differences.

    Values are rounded when written (%.3f, p-values %.3e), so floating point
    differences between platforms can flip the last digit: fixed-point
    values may differ by one unit in the last digit, p-values by 0.15%.
    Ranked files list genes by score, and genes with tied scores can be
    ordered differently, so rows are matched by gene.
    """
    with open(reference_file) as ref, open(output_file) as out:
        header = [next(ref), next(ref)]
        assert header == [next(out), next(out)], f"{output_file}: header differs"

    if reference_file.endswith(".pvalues.gct"):
        tolerance = {"rtol": 1.5e-3, "atol": 0}
    else:
        tolerance = {"rtol": 0, "atol": 1.5e-3}
    ref_df = pd.read_csv(reference_file, sep="\t", skiprows=2, index_col="Description")
    out_df = pd.read_csv(output_file, sep="\t", skiprows=2, index_col="Description")
    pd.testing.assert_frame_equal(
        ref_df.sort_index(),
        out_df.sort_index(),
        check_exact=False,
        obj=output_file,
        **tolerance,
    )


@pytest.fixture(scope="session")
def load_kotliarov2020_processed_data() -> AnnData:
    return bc.datasets.Kotliarov2020_processed()


@pytest.fixture(scope="session")
def load_kotliarov2020_raw_data() -> AnnData:
    return bc.datasets.Kotliarov2020_raw()


@pytest.fixture(scope="session")
def subset_kotliarov2020_raw_data(load_kotliarov2020_raw_data: AnnData) -> AnnData:
    subset = sample(
        [i for i in range(load_kotliarov2020_raw_data.shape[0])],
        pytest.raw_data_subset_size,
    )
    return load_kotliarov2020_raw_data[subset, :]


@pytest.fixture(scope="session")
def create_random_anndata_object() -> AnnData:

    counts = csr_matrix(np.random.poisson(1, size=(1000, 2000)))
    return anndata.AnnData(counts)


@pytest.fixture
def reference_files() -> List[str]:

    file_list: List[str] = [
        "WilxRank.gct",
        "WilxRank.logFC.gct",
        "WilxRank.pvalues.gct",
        "average.gct",
        "cell2labels.tsv",
        "fract_pos.gct",
        "labelinfo.tsv",
        "celltype_labelinfo.tsv",
    ]
    return file_list


@pytest.fixture
def reference_folder() -> str:

    root_path: str = (
        "./tests/data/st/wrapper_funcs/labeling_test_CTL_OLD/labelings/celltype/"
    )

    return root_path


def test_additional_labeling(
    load_kotliarov2020_processed_data: AnnData,
    reference_files: List[str],
    reference_folder: str,
    tmp_path_factory: pathlib.Path,
):

    tmp_path = tmp_path_factory.mktemp("additional_labeling")

    additional_labeling(
        adata=load_kotliarov2020_processed_data,
        labeling_author="MK",
        results_folder=join(tmp_path),
        labeling_to_use="celltype1",
        labeling_name="celltype",
        labeling_description="ctl_new",
        is_celltype_labeling=False,
    )

    additional_labeling(
        adata=load_kotliarov2020_processed_data,
        labeling_author="MK",
        results_folder=join(tmp_path),
        labeling_to_use="celltype1",
        labeling_name="celltype",
        labeling_description="manual celltype annotation",
        is_celltype_labeling=True,
        filename="celltype_labelinfo.tsv",
    )

    path_to_output_files = join(tmp_path, "labelings", "celltype")

    for file in reference_files:
        reference_file: str = join(reference_folder, file)
        output_file: str = join(path_to_output_files, file)

        if file.endswith(".gct"):
            assert_gct_equal(reference_file, output_file)
        else:
            identical: bool = filecmp.cmp(reference_file, output_file, shallow=False)
            assert identical is True, f"File {file} is not as expected!"


def test_highly_variable_genes(create_random_anndata_object: AnnData):

    adata = highly_variable_genes(
        adata=create_random_anndata_object, batch_key=None, n_shared=2
    )
    # ensure that we get all hvgs and not just one
    assert adata.shape[0] == pytest.raw_data_subset_size
    assert adata.shape[1] == 608


def test_all_symbols_ok():

    good_symbols_list: List[str] = [
        "AC004.2",
        "AC43.3",
        "A1",
        "AC004.2",
        "ACDF3",
        "A..4A",
        "YYYY",
        "782387324A",
        "Nature",
        "GoodgenedNA",
        "NAGene",
        "abcNAcde",
        "nature",
        "HLA-DOB",
        "Mafa-DQA",
        "-ABC",
    ]
    good: bool = all_symbols_ok(symbols_list=good_symbols_list)
    assert good == True

    bad_symbols_list: List[str] = [
        "434443",
        "NA",
        ".",
        "...",
        ".1",
        ".NA",
        "",
        "na",
        "Na",
        "N/A",
        "CORRECTGeneName",
    ]
    bad: bool = all_symbols_ok(symbols_list=bad_symbols_list)
    assert bad == False
