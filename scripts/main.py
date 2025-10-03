import glob
from pydna.parsers import parse
from Bio.Restriction import BsaI
from pydna.dseqrecord import Dseqrecord
from pydna.dseq import Dseq
from Bio import BiopythonParserWarning
import yaml
import pandas as pd
import networkx as nx

# Supress BiopythonParserWarning:
import warnings

warnings.filterwarnings("ignore", category=BiopythonParserWarning)

odc_plasmids = pd.read_csv("../odc_plasmids.csv")

with open("syntax.yaml", "r") as restriction_fragment:
    syntax = yaml.safe_load(restriction_fragment)

overhangs = [edge["overhang"] for edge in syntax["edges"]]


gb_files = glob.glob("../**/*.gb", recursive=True)

# Exclude the genbank dir
gb_files = [f for f in gb_files if not f.startswith("../genbank/")]
enzyme = BsaI
ovhg = abs(enzyme.ovhg)

# gb_files = ["./Ecoli Nanobody Toolkit/genbank_seq/BC_RJ_SD8.gb"]
table = list()
for gb_file in gb_files:
    collection = gb_file.split("/")[-3]
    plasmid_id = gb_file.split("/")[-1].split(".")[0]
    left_ovhg = None
    right_ovhg = None
    longest_feature_type = None
    info = None
    seq = parse(gb_file)[0]
    if not seq.circular:
        seq = seq.looped()

    out = seq.cut(enzyme)
    found_edges = []
    if len(out) > 1:
        for restriction_fragment in out:
            for rc in [False, True]:
                if rc:
                    restriction_fragment = restriction_fragment.reverse_complement()
                _left_ovhg = str(restriction_fragment.seq[:ovhg]).upper()
                _right_ovhg = str(restriction_fragment.seq[-ovhg:]).upper()
                left_edge = (
                    overhangs.index(_left_ovhg) if _left_ovhg in overhangs else None
                )
                right_edge = (
                    overhangs.index(_right_ovhg) if _right_ovhg in overhangs else None
                )
                if (
                    left_edge is not None
                    and right_edge is not None
                    and left_edge < right_edge
                ):
                    found_edges.append((_left_ovhg, _right_ovhg, restriction_fragment))
    if len(found_edges) > 1:
        info = "Multiple edges found"
    elif len(found_edges) == 0:
        info = "No edges found"
    else:
        fragment = found_edges[0][2]
        left_ovhg = found_edges[0][0]
        right_ovhg = found_edges[0][1]
        sorted_features = sorted(fragment.features, key=lambda x: len(x.location))
        if len(sorted_features) > 0:
            longest_feature_type = sorted_features[0].type
        else:
            info = "No features found"

    table.append(
        {
            "collection": collection,
            "id": plasmid_id,
            "left_overhang": left_ovhg,
            "right_overhang": right_ovhg,
            "longest_feature_type": longest_feature_type,
            "info": info,
            "path": gb_file,
        }
    )
df = pd.DataFrame(table)
df["plasmid_name"] = df["id"].map(odc_plasmids.set_index("ODC ID")["Name"])
# Sort columns
df = df[
    [
        "collection",
        "id",
        "plasmid_name",
        "left_overhang",
        "right_overhang",
        "longest_feature_type",
        "info",
        "path",
    ]
]
df.to_csv("summary.csv", index=False)
df.to_json("index.json", indent=4, orient="records")

df = df[df["info"].isna()]
df.to_json("index_overhangs.json", indent=4, orient="records")

# Get all possible pairs of overhangs
pairs = list(set(df["left_overhang"] + "|" + df["right_overhang"]))

G = nx.DiGraph()
G.add_nodes_from(overhangs)
for pair in pairs:
    left_overhang, right_overhang = pair.split("|")
    G.add_edge(left_overhang, right_overhang)

# Get all linear paths from the first to the last node
paths = list(nx.all_simple_paths(G, source=overhangs[0], target=overhangs[-1]))
