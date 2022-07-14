import pathlib
import pytest
from typing import List
from os import listdir
from os.path import isfile, join
import filecmp

from scvelo import AnnData

import besca as bc
from besca.st._wrapper_funcs import additional_labeling_refactored


@pytest.fixture(scope="session")
def load_kotliarov2020_processed_data() -> AnnData:
    return bc.datasets.Kotliarov2020_processed()


@pytest.fixture
def reference_files() -> List[str]:

    file_list: List[str] = [
        # "WilxRank.gct",
        # "WilxRank.logFC.gct",
        # "WilxRank.pvalues.gct",
        "average.gct",
        "cell2labels.tsv",
        "fract_pos.gct",
        # "labelinfo.tsv",
    ]
    return file_list


@pytest.fixture
def reference_folder() -> str:

    root_path: str = (
        "./tests/data/st/wrapper_funcs/labeling_test_CTL_OLD/labelings/celltype/"
    )

    return root_path


def test_additional_labeling_refactored(
    load_kotliarov2020_processed_data: AnnData,
    reference_files: List[str],
    reference_folder: str,
    tmp_path_factory: pathlib.Path,
):

    tmp_path = tmp_path_factory.mktemp("additional_labeling")

    additional_labeling_refactored(
        adata=load_kotliarov2020_processed_data,
        labeling_author="MK",
        results_folder=join(tmp_path),
        labeling_to_use="celltype1",
        labeling_name="celltype",
        labeling_description="ctl_new",
        is_celltype_labeling=True,
    )

    path_to_output_files = join(tmp_path, "labelings", "celltype")

    for file in reference_files:
        reference_file: str = join(reference_folder, file)
        output_file: str = join(path_to_output_files, file)

        print(f"Comparing {reference_file} and {output_file})")
        identical: bool = filecmp.cmp(reference_file, output_file)

        # TODO: Print more info on failure
        assert identical is True
