import pandas as pd
import pytest

from scripts.add_radicality import annotate_radicality


def test_adds_published_grantham_and_miyata_distances() -> None:
    data = pd.DataFrame({"anc_aa": ["A", "D"], "der_aa": ["R", "E"]})

    result = annotate_radicality(data)

    assert result["gr_rad_anc"].tolist() == [112.0, 45.0]
    assert result["miyata_rad"].tolist() == [2.92, 0.9]


def test_rejects_nonstandard_amino_acids() -> None:
    data = pd.DataFrame({"anc_aa": ["X"], "der_aa": ["A"]})

    with pytest.raises(ValueError, match="Non-standard amino acid"):
        annotate_radicality(data)
