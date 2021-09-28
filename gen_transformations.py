from itertools import combinations

FRAMES = ["ECEF", "TEME", "PEF", "TOD", "MOD", "J2000"]
IMPLED_TRANS = [("ECEF", "PEF"), ("PEF", "TOD"), ("TOD", "MOD"), ("MOD", "J2000"), ("TEME", "PEF")]

OUT_FILE = "test_out.cpp"

TRANSFORMATION_TEMPLATE = """
template <>
Transformation<Frame::{}, Frame::{}>::Transformation(const double jd, const EOPData eop) {{
{}
}}
"""

INV_TRANS_TEMPLATE = """\
    auto const inverse = Transformation<Frame::{}, Frame::{}>(jd, eop);
    *this = inverse.inverse();\
"""

out = ["\n// *** AUTO-GENERATED START ***\n"]

frame_combs = [("MOD", "J2000"), ("TOD", "MOD"), ("PEF", "TOD"), ("ECEF", "PEF"), ("TEME", "PEF")]
frame_combs += [
    ("TOD", "J2000"), ("PEF", "J2000"), ("ECEF", "J2000"), ("TEME", "J2000"),
    ("PEF", "MOD"), ("ECEF", "MOD"), ("TEME", "MOD"),
    ("ECEF", "TOD"), ("TEME", "TOD"),
    ("TEME", "ECEF"),
]

# frame_combs = list(IMPLED_TRANS)
# for from_to in combinations(FRAMES, 2):
#     if from_to in IMPLED_TRANS or from_to[::-1] in IMPLED_TRANS:
#         continue
#     frame_combs.append(from_to)
#     # Special case
#     if "ECEF" in from_to and "TEME" in from_to:
#         func_txt = "    // FIXME: Impl special case"
#     else:
#     out.append(TRANSFORMATION_TEMPLATE.format(*from_to, func_txt))

# Inverse transformations
out.append("\n// Inverse Transformations\n")
for from_to in frame_combs:
    inv_code = INV_TRANS_TEMPLATE.format(*from_to)
    out.append(TRANSFORMATION_TEMPLATE.format(*from_to[::-1], inv_code))

out.append("\n// *** AUTO-GENERATED END ***\n")
out_str = "".join(out)
with open(OUT_FILE, "w") as f:
    f.write(out_str)
