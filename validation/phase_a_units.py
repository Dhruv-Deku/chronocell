"""
The Phase A units (validation/phase_a.py), fixed before any of them was built. Practice units come
only from practice datasets, test units only from test datasets (validation/datasets.py roles); a
test unit is built only when a pre-registered gate reads it.

Hi-C pairings: each imaged region with Hi-C of the same cell line. Bintu practice regions get Rao 2014
(K562 only) and ENCODE in situ / intact Hi-C; depth variants by binomial thinning (1, 1/4, 1/16) let a
calibration learn how size depends on depth. Bintu test regions get Rao 2014 (IMR-90) and, as
secondary units, ENCODE IMR-90 intact and dilution Hi-C and ENCODE A549 in situ Hi-C.
"""

from __future__ import annotations

from phase_a import UnitSpec, genome_units

THIN = (1.0, 0.25, 0.0625)

PRACTICE_HIC: list[UnitSpec] = (
    [UnitSpec("bintu_k562_28_30", "hic", s, t) for s in ("rao2014", "encode_k562_intact") for t in THIN]
    + [UnitSpec(k, "hic", s, t) for k in ("bintu_hct116_28_30", "bintu_hct116_34_37")
       for s in ("encode_hct116_insitu", "encode_hct116_intact") for t in THIN]
    + [UnitSpec("su_chr2", "hic", "rao2014", t) for t in THIN]
    + [UnitSpec("su_chr2_parm_rep", "hic", "rao2014", 1.0)]
)
PRACTICE_IMAGING: list[UnitSpec] = [UnitSpec(k, "imaging") for k in (
    "bintu_k562_28_30", "bintu_hct116_28_30", "bintu_hct116_28_30_auxin", "bintu_hct116_34_37", "su_chr2",
    "su_chr2_parm_rep")]


def practice_genome() -> list[UnitSpec]:
    return (genome_units("su_genome_tx", "hic") + genome_units("su_genome_tx", "hic", thin=0.25)
            + genome_units("su_genome_tx", "imaging"))


TEST_HIC_MAIN: list[UnitSpec] = [UnitSpec(k, "hic", "rao2014") for k in (
    "bintu_imr90_28_30", "bintu_imr90_18_20", "su_chr21", "su_chr21_rep")]
TEST_HIC_SECONDARY: list[UnitSpec] = (
    [UnitSpec(k, "hic", s) for k in ("bintu_imr90_28_30", "bintu_imr90_18_20")
     for s in ("encode_imr90_intact", "encode_imr90_dilution")]
    + [UnitSpec("bintu_a549_28_30", "hic", "encode_a549_insitu")]
)
TEST_IMAGING: list[UnitSpec] = [UnitSpec(k, "imaging") for k in (
    "bintu_imr90_28_30", "bintu_imr90_18_20", "bintu_a549_28_30", "bintu_hct116_34_37_auxin", "su_chr21",
    "su_chr21_rep")]


def test_genome(input_kind: str) -> list[UnitSpec]:
    return genome_units("su_genome", input_kind) + genome_units("su_genome_amanitin", input_kind)


def __getattr__(name):              # PRACTICE / TEST are built lazily: listing the genome units loads the traces
    if name == "PRACTICE":
        return PRACTICE_HIC + PRACTICE_IMAGING + practice_genome()
    if name == "TEST":
        return TEST_HIC_MAIN + TEST_HIC_SECONDARY + TEST_IMAGING + test_genome("hic") + test_genome("imaging")
    raise AttributeError(name)
