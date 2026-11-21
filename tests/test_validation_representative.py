"""pick_representative_cell must return one cell even when barcodes repeat."""

import pandas as pd

from tcrsift.validation import pick_representative_cell


def test_repeated_barcode_still_gives_one_row():
    # The same barcode in two samples of one patient, as load_samples keeps it.
    cells = pd.DataFrame(
        {"TRA_1_umis": [9, 2, 1], "TRB_1_umis": [9, 2, 1], "TRA_1_v_gene": ["TRAV12-2", "TRAV12-2", "TRAV1"]},
        index=["AAAC-1", "AAAC-1", "GGGT-1"],
    )
    rep = pick_representative_cell(cells, umi_cols=("TRA_1_umis", "TRB_1_umis"))
    assert isinstance(rep, pd.Series)
    assert rep["TRA_1_v_gene"] == "TRAV12-2"
    assert rep["TRA_1_umis"] == 9
