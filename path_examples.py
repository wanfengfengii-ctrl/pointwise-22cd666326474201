"""Fixed input/output examples; no finite-path or limit solver is provided.

Default: verify the supplied V1 and the hand-calculated inputs.
--v2: exercise the new public APIs after implementation.
"""
from fractions import Fraction as F
import sys

from spectral.archive.catalog import Catalog
from spectral.archive.codec import encode_payload
from spectral.archive.reading import read
from spectral.dictionary.revisions import save_dictionary
from spectral.measurements.intake import save_batch
from spectral.measurements.cohorts import freeze
from spectral.unmixing.jobs import unmix
from spectral.unmixing.solver import solve_canonical


def segment(lo, hi, constant, slope):
    return {"lo": F(lo), "hi": F(hi),
            "constant": [F(x) for x in constant],
            "slope": [F(x) for x in slope]}


TRIANGLE = [
    segment(0, F(1, 2), [1, 0], [-2, 1]),
    segment(F(1, 2), 1, [0, F(1, 2)], [0, 0]),
]

CASES = [
    {"name": "finite_boundary", "matrix": [[1, 1], [0, 1]],
     "values": [[1, 0], [0, 1]], "direction": [-1, 1], "lo": 0, "hi": 1,
     "paths": [TRIANGLE, [segment(0, 1, [0, F(1, 2)], [0, 0])]],
     "limits": [
         ({}, [[[F(0), F(1)]], [[F(0), F(1)]]]),
         ({"k1": (None, F(3, 4))}, [[[F(0), F(1)]], [[F(0), F(1)]]]),
         ({"k1": (F(1, 4), F(1, 4))}, [[[F(1, 4), F(1, 4)]], []]),
         ({"k1": (F(1, 2), F(1, 2))}, [[[F(1, 2), F(1)]], [[F(0), F(1)]]]),
         ({"k1": (F(3, 4), None)}, [[], []]),
         ({"k0": (0, F(1, 2)), "k1": (F(1, 4), None)},
          [[[F(1, 4), F(1)]], [[F(0), F(1)]]]),
     ]},
    {"name": "duplicate_simultaneous", "matrix": [[1, 1]],
     "values": [[1]], "direction": [-2], "lo": 0, "hi": 1,
     "paths": [[segment(0, F(1, 2), [F(1, 2), F(1, 2)], [-1, -1]),
                segment(F(1, 2), 1, [0, 0], [0, 0])]],
     "limits": [({"k0": (None, 0)}, [[[F(1, 2), F(1)]]])]},
    {"name": "proportional", "matrix": [[1, 2]],
     "values": [[1]], "direction": [-2], "lo": 0, "hi": 1,
     "paths": [[segment(0, F(1, 2), [F(1, 5), F(2, 5)], [-F(2, 5), -F(4, 5)]),
                segment(F(1, 2), 1, [0, 0], [0, 0])]],
     "limits": [({}, [[[F(0), F(1)]]])]},
    {"name": "zero_matrix", "matrix": [[0, 0]],
     "values": [[2]], "direction": [1], "lo": 0, "hi": 2,
     "paths": [[segment(0, 2, [0, 0], [0, 0])]],
     "limits": [({"k1": (0, 0)}, [[[F(0), F(2)]]])]},
    {"name": "singleton", "matrix": [[1, 1], [0, 1]],
     "values": [[1, 0]], "direction": [-1, 1], "lo": F(3, 4), "hi": F(3, 4),
     "paths": [[segment(F(3, 4), F(3, 4), [0, F(1, 2)], [0, 0])]],
     "limits": [({"k1": (F(1, 2), F(1, 2))}, [[[F(3, 4), F(3, 4)]]])]},
    {"name": "nonzero_left", "matrix": [[1, 1], [0, 1]],
     "values": [[1, 0]], "direction": [-1, 1], "lo": F(1, 4), "hi": F(3, 4),
     "paths": [[segment(F(1, 4), F(1, 2), [1, 0], [-2, 1]),
                segment(F(1, 2), F(3, 4), [0, F(1, 2)], [0, 0])]],
     "limits": [({}, [[[F(1, 4), F(3, 4)]]])]},
]


def make_parent(c, case):
    m, n = len(case["matrix"]), len(case["matrix"][0])
    channels = ["c" + str(i) for i in range(m)]
    components = ["k" + str(i) for i in range(n)]
    dictionary = save_dictionary(c, channels, components,
                                 ["signal"] * m, ["mg/L"] * n, case["matrix"])
    batch = save_batch(c, channels, ["signal"] * m,
                       [{"sample_id": "s" + str(i), "values": y}
                        for i, y in enumerate(case["values"])])
    revisions = c.read_record(batch, "batch")["data"]["revision_ids"]
    return unmix(c, freeze(c, dictionary, revisions)), revisions


def verify_v1_inputs():
    checked = 0
    c = Catalog(":memory:")
    try:
        for case in CASES:
            parent, revisions = make_parent(c, case)
            rows = read(c, parent)["data"]["rows"]
            assert [r["measurement_revision"] for r in rows] == revisions
            for y, pieces in zip(case["values"], case["paths"]):
                for piece in pieces:
                    for t in (piece["lo"], (piece["lo"] + piece["hi"]) / 2, piece["hi"]):
                        measured = [F(a) + t * F(b) for a, b in zip(y, case["direction"])]
                        actual = solve_canonical(case["matrix"], measured).coefficients
                        expected = tuple(a + t * b for a, b in
                                         zip(piece["constant"], piece["slope"]))
                        assert actual == expected, (case["name"], t, expected, actual)
                        checked += 1
        print(encode_payload({"mode": "V1_input_check", "point_checks": checked,
                              "target_implemented": False,
                              "note": "Hand examples agree with V1; not a complete-path proof."}))
    finally:
        c.close()


def verify_v2_examples():
    from spectral.response.requests import trace
    from spectral.response.limits import acceptable
    c = Catalog(":memory:")
    try:
        for case in CASES:
            parent, revisions = make_parent(c, case)
            path_id = trace(c, parent, case["direction"], case["lo"], case["hi"])
            saved = read(c, path_id)
            assert saved["kind"] == "path"
            assert saved["parents"] == [{"type": "artifact", "id": parent}]
            data = saved["data"]
            assert data["parent_id"] == parent
            assert data["direction"] == [F(x) for x in case["direction"]]
            assert (data["lo"], data["hi"]) == (F(case["lo"]), F(case["hi"]))
            assert data["components"] == ["k" + str(i) for i in range(len(case["matrix"][0]))]
            assert [r["measurement_revision"] for r in data["rows"]] == revisions
            assert [r["sample_id"] for r in data["rows"]] == ["s" + str(i) for i in range(len(revisions))]
            assert [r["segments"] for r in data["rows"]] == case["paths"], case["name"]
            assert trace(c, parent, case["direction"], case["lo"], case["hi"]) == path_id
            reports = []
            for bounds, expected in case["limits"]:
                result_id = acceptable(c, path_id, bounds)
                result = read(c, result_id)
                assert result["kind"] == "acceptable"
                assert result["parents"] == [{"type": "artifact", "id": path_id}]
                assert result["data"]["path_id"] == path_id
                assert result["data"]["bounds"] == {
                    k: [None if v is None else F(v) for v in pair]
                    for k, pair in bounds.items()}
                rows = result["data"]["rows"]
                assert [r["measurement_revision"] for r in rows] == revisions
                assert [r["sample_id"] for r in rows] == ["s" + str(i) for i in range(len(revisions))]
                assert [r["intervals"] for r in rows] == expected, case["name"]
                assert acceptable(c, path_id, bounds) == result_id
                reports.append(result)
            print(encode_payload({"case": case["name"], "path": saved, "reports": reports}))
    finally:
        c.close()


if __name__ == "__main__":
    if sys.argv[1:] == []:
        verify_v1_inputs()
    elif sys.argv[1:] == ["--v2"]:
        verify_v2_examples()
    else:
        raise SystemExit("usage: python3 -B path_examples.py [--v2]")
