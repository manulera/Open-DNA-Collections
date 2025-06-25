import glob
from pydna.parsers import parse
from Bio.Restriction import BsaI
from pydna.dseqrecord import Dseqrecord
from pydna.dseq import Dseq
from Bio import BiopythonParserWarning

# Supress BiopythonParserWarning:
import warnings

warnings.filterwarnings("ignore", category=BiopythonParserWarning)


edges = [
    "GGAG",
    "TGAC",
    "TCCC",
    "TACT",
    "CCAT",
    "AATG",
    "AGCC",
    "TTCG",
    "GCAG",
    "GCTT",
    "GGTA",
    "CGCT",
]


gb_files = glob.glob("./**/*.gb", recursive=True)

# Exclude the genbank dir
gb_files = [f for f in gb_files if not f.startswith("./genbank/")]
enzyme = BsaI
ovhg = abs(enzyme.ovhg)

# gb_files = ["./Ecoli Nanobody Toolkit/genbank_seq/BC_RJ_SD8.gb"]

for gb_file in gb_files:
    seq = parse(gb_file)[0]
    if not seq.circular:
        seq = seq.looped()

    out = seq.cut(enzyme)
    if len(out) > 1:
        found_edges = []
        for f in out:
            for rc in [False, True]:
                if rc:
                    f = f.reverse_complement()
                left_ovhg = str(f.seq[:ovhg]).upper()
                right_ovhg = str(f.seq[-ovhg:]).upper()
                left_edge = edges.index(left_ovhg) if left_ovhg in edges else None
                right_edge = edges.index(right_ovhg) if right_ovhg in edges else None
                if (
                    left_edge is not None
                    and right_edge is not None
                    and left_edge < right_edge
                ):
                    found_edges.append((left_edge, right_edge, rc))
        if len(found_edges):
            print(gb_file, *found_edges)
            # print(">>", gb_file, found_edges)
        else:
            print(" ", gb_file, "--")
