"""Independent audit of the complete finite path and limit consumption.

The checks here deliberately do not reuse the path producer's affine-face
solver or its sign routines.  Candidate supports are solved with Cramer's rule
and explicit determinants (the same independent technique already used by
tests.test_unmix), and every claim is proved over a whole closed segment with
exact quadratic analysis rather than sampled at a few parameters.

Completeness argument (short form): the canonical solution is unique and
continuous in t.  On a fixed support it is affine with rational coefficients.
Therefore its formula can change only at t where some coefficient of a
candidate formula hits zero, or where two candidate affine formulas coincide.
The audit recomputes all such rational candidates from the saved data and
requires each interior one to be a stored breakpoint; each segment is then
proved lexicographically minimal on its whole closed span.
"""
import itertools
from fractions import Fraction as F
from math import isqrt
from unittest import mock

from spectral.archive.catalog import Catalog
from spectral.archive.reading import read
from spectral.dictionary.revisions import save_dictionary
from spectral.measurements.cohorts import freeze
from spectral.measurements.intake import save_batch
from spectral.measurements.revisions import append_measurement_revision
from spectral.response.limits import acceptable
from spectral.response.requests import trace
from spectral.unmixing.jobs import unmix
from tests.test_unmix import det, cramer, inner, publish

import unittest


# ---------------------------------------------------------------------------
# independent affine candidates
# ---------------------------------------------------------------------------
class Candidate:
    def __init__(self, matrix, y0, y1, support, constant, slope, feasible_lo, feasible_hi):
        self.matrix = matrix
        self.y0 = y0
        self.y1 = y1
        self.support = tuple(support)
        self.constant = tuple(constant)
        self.slope = tuple(slope)
        self.feasible_lo = feasible_lo
        self.feasible_hi = feasible_hi
        self.residual = square_polynomial(
            tuple(a - b for a, b in zip(apply_matrix(matrix, constant), y0)),
            tuple(a - b for a, b in zip(apply_matrix(matrix, slope), y1)))
        self.norm = square_polynomial(constant, slope)


def apply_matrix(matrix, coefficients):
    return tuple(sum((row[j] * coefficients[j] for j in range(len(coefficients))), F())
                 for row in matrix)


def square_polynomial(constant, slope):
    return (inner(constant, constant), 2 * inner(constant, slope),
            inner(slope, slope))


def at(polynomial, t):
    return polynomial[0] + t * (polynomial[1] + t * polynomial[2])


def _independent_local(matrix, y0, y1, support):
    m = len(matrix)
    columns = [tuple(F(matrix[i][j]) for i in range(m)) for j in support]
    k = len(columns)
    if k == 0:
        return (), ()
    gram = [[inner(a, b) for b in columns] for a in columns]
    if not det(gram):
        return None
    c0 = cramer(gram, [inner(column, y0) for column in columns])
    c1 = cramer(gram, [inner(column, y1) for column in columns])
    return c0, c1


def independent_candidates(matrix, y0, y1):
    m, n = len(matrix), len(matrix[0])
    result = []
    for size in range(n + 1):
        for support in itertools.combinations(range(n), size):
            local = _independent_local(matrix, y0, y1, support)
            if local is None:
                continue
            local0, local1 = local
            constant = [F() for _ in range(n)]
            slope = [F() for _ in range(n)]
            for position, value in zip(support, local0):
                constant[position] = value
            for position, value in zip(support, local1):
                slope[position] = value
            lo = hi = None
            usable = True
            for position in support:
                base, rate = constant[position], slope[position]
                if rate == 0:
                    if base < 0:
                        usable = False
                        break
                elif rate > 0:
                    bound = -base / rate
                    lo = bound if lo is None else max(lo, bound)
                else:
                    bound = -base / rate
                    hi = bound if hi is None else min(hi, bound)
            if not usable:
                continue
            if lo is not None and hi is not None and lo > hi:
                continue
            candidate = Candidate(matrix, y0, y1, support, constant, slope, lo, hi)
            result.append(candidate)
    return result


# ---------------------------------------------------------------------------
# exact whole-interval predicates
# ---------------------------------------------------------------------------
def nonnegative_on(polynomial, lo, hi):
    smallest = min(at(polynomial, lo), at(polynomial, hi))
    a0, a1, a2 = polynomial
    if a2 > 0:
        vertex = -a1 / (2 * a2)
        if lo < vertex < hi:
            smallest = min(smallest, at(polynomial, vertex))
    return smallest >= 0


def rational_zeros_in(polynomial, lo, hi):
    """Exact zeros of a rational quadratic, used only for tie resolution."""
    a0, a1, a2 = polynomial
    if a2 == 0:
        if a1:
            root = -a0 / a1
            return [root] if lo <= root <= hi else []
        return [lo, hi] if a0 == 0 else []
    discriminant = a1 * a1 - 4 * a2 * a0
    if discriminant < 0:
        return []
    numerator, denominator = discriminant.numerator, discriminant.denominator
    sn, sd = isqrt(numerator), isqrt(denominator)
    if sn * sn != numerator or sd * sd != denominator:
        return []
    square = F(sn, sd)
    roots = [(-a1 - square) / (2 * a2), (-a1 + square) / (2 * a2)]
    return sorted({root for root in roots if lo <= root <= hi})


def feasible_span(candidate, lo, hi):
    left = lo if candidate.feasible_lo is None else max(lo, candidate.feasible_lo)
    right = hi if candidate.feasible_hi is None else min(hi, candidate.feasible_hi)
    return (left, right) if left <= right else None


def lexicographically_beats(winner_residual, winner_norm, other_residual, other_norm, lo, hi):
    difference = tuple(a - b for a, b in zip(other_residual, winner_residual))
    if any(difference):
        if not nonnegative_on(difference, lo, hi):
            return False
        tie_norm = tuple(a - b for a, b in zip(other_norm, winner_norm))
        return all(at(tie_norm, point) >= 0
                   for point in rational_zeros_in(difference, lo, hi))
    difference = tuple(a - b for a, b in zip(other_norm, winner_norm))
    return nonnegative_on(difference, lo, hi)


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
def run_path(catalog, matrix, values, direction, lo, hi, sample_ids=None):
    parent, revisions = make_parent(catalog, matrix, values, sample_ids)
    return parent, revisions, trace(catalog, parent, direction, lo, hi)


def make_parent(catalog, matrix, values, sample_ids=None):
    channels = ["c" + str(i) for i in range(len(matrix))]
    components = ["k" + str(i) for i in range(len(matrix[0]))]
    dictionary = save_dictionary(catalog, channels, components,
                                 ["signal"] * len(channels),
                                 ["mg/L"] * len(components), matrix)
    sample_ids = sample_ids or ["s" + str(i) for i in range(len(values))]
    batch = save_batch(catalog, channels, ["signal"] * len(channels),
                       [{"sample_id": sid, "values": y}
                        for sid, y in zip(sample_ids, values)])
    revisions = catalog.read_record(batch, "batch")["data"]["revision_ids"]
    return unmix(catalog, freeze(catalog, dictionary, revisions)), revisions


def audit_row(test, matrix, y, direction, segments, lo, hi):
    operator_values = [F(value) for value in y]
    direction = [F(value) for value in direction]
    candidates = independent_candidates(matrix, operator_values, direction)
    test.assertTrue(candidates)
    # exact coverage and merge discipline
    test.assertEqual(segments[0]["lo"], lo)
    test.assertEqual(segments[-1]["hi"], hi)
    for left, right in zip(segments, segments[1:]):
        test.assertEqual(left["hi"], right["lo"])
        test.assertNotEqual((left["constant"], left["slope"]),
                            (right["constant"], right["slope"]))
        left_value = [c + left["hi"] * s for c, s in zip(left["constant"], left["slope"])]
        right_value = [c + right["lo"] * s for c, s in zip(right["constant"], right["slope"])]
        test.assertEqual(left_value, right_value)
    for segment in segments:
        piece_lo, piece_hi = segment["lo"], segment["hi"]
        residual = square_polynomial(
            tuple(a - b for a, b in zip(apply_matrix(matrix, segment["constant"]),
                                        operator_values)),
            tuple(a - b for a, b in zip(apply_matrix(matrix, segment["slope"]), direction)))
        norm = square_polynomial(segment["constant"], segment["slope"])
        for base, rate in zip(segment["constant"], segment["slope"]):
            if rate == 0:
                test.assertGreaterEqual(base, 0)
            else:
                test.assertGreaterEqual(base + rate * piece_lo, 0)
                test.assertGreaterEqual(base + rate * piece_hi, 0)
        if piece_lo == piece_hi:
            continue
        # Whole-piece lexicographic proof against every independently computed
        # support, restricted to that support's own feasible overlap.  The
        # canonical solution is unique, so a winner proved minimal on the whole
        # piece leaves no room for a hidden switch inside it.
        for candidate in candidates:
            span = feasible_span(candidate, piece_lo, piece_hi)
            if span is None:
                continue
            test.assertTrue(lexicographically_beats(
                residual, norm, candidate.residual, candidate.norm,
                span[0], span[1]),
                (matrix, y, direction, segment, candidate.support))


class TestFinitePaths(unittest.TestCase):
    def test_hand_examples_full_span(self):
        cases = [
            ([[1, 1], [0, 1]], [[1, 0], [0, 1]], [-1, 1], F(0), F(1)),
            ([[1, 1]], [[1]], [-2], F(0), F(1)),
            ([[1, 2]], [[1]], [-2], F(0), F(1)),
            ([[0, 0]], [[2]], [1], F(0), F(2)),
            ([[1, 1], [0, 1]], [[1, 0]], [-1, 1], F(3, 4), F(3, 4)),
            ([[1, 1], [0, 1]], [[1, 0]], [-1, 1], F(1, 4), F(3, 4)),
            ([[1, 2, 0], [2, 4, 1]], [[0, 3]], [2, 2], F(0), F(3)),
            ([[1, 0, 1], [0, 1, 1]], [[2, 1]], [-1, -1], F(0), F(1)),
        ]
        c = Catalog(":memory:")
        try:
            for matrix, values, direction, lo, hi in cases:
                parent, revisions, path_id = run_path(c, matrix, values, direction, lo, hi)
                data = read(c, path_id)["data"]
                self.assertEqual(data["parent_id"], parent)
                self.assertEqual([row["measurement_revision"] for row in data["rows"]],
                                 revisions)
                for row, y in zip(data["rows"], values):
                    audit_row(self, matrix, y, direction, row["segments"], lo, hi)
        finally:
            c.close()

    def test_random_small_full_span(self):
        grid = [
            ([[1, 1], [0, 1]], [[1, 0], [0, 2], [2, 1]], [-1, 1], F(0), F(2)),
            ([[1, 2]], [[0], [1], [3]], [-1], F(0), F(3)),
            ([[0, 1, 0], [0, 0, 0]], [[1, 0], [0, 1], [2, 2]], [-1, 1], F(0), F(2)),
            ([[1, 2, 0], [2, 4, 1]], [[0, 3], [1, 1]], [2, -1], F(0), F(1)),
            ([[1, 0, 1], [0, 1, 1]], [[1, 1], [0, 0], [2, 2]], [-1, 2], F(0), F(1)),
            ([[1, 1, 1]], [[1], [2]], [1], F(0), F(1)),
        ]
        c = Catalog(":memory:")
        try:
            for matrix, values, direction, lo, hi in grid:
                parent, _, path_id = run_path(c, matrix, values, direction, lo, hi)
                for row, y in zip(read(c, path_id)["data"]["rows"], values):
                    audit_row(self, matrix, y, direction, row["segments"], lo, hi)
        finally:
            c.close()

    def test_singleton_parent_interval(self):
        c = Catalog(":memory:")
        try:
            _, _, path_id = run_path(c, [[1, 1], [0, 1]], [[1, 0]], [-1, 1],
                                     F(3, 4), F(3, 4))
            data = read(c, path_id)["data"]
            (segment,) = data["rows"][0]["segments"]
            self.assertEqual((segment["lo"], segment["hi"]), (F(3, 4), F(3, 4)))
            self.assertEqual(segment["slope"], [F(0), F(0)])
            self.assertEqual(segment["constant"], [F(0), F(1, 2)])
            report = acceptable(c, path_id, {})
            self.assertEqual(read(c, report)["data"]["rows"][0]["intervals"],
                             [[F(3, 4), F(3, 4)]])
        finally:
            c.close()

    def test_record_identity_and_parent_binding(self):
        c = Catalog(":memory:")
        try:
            parent, revisions, path_id = run_path(
                c, [[1, 1], [0, 1]], [[1, 0], [0, 1]], [-1, 1], 0, 1,
                sample_ids=["second", "first"])
            data = read(c, path_id)["data"]
            self.assertEqual([row["sample_id"] for row in data["rows"]],
                             ["second", "first"])
            self.assertEqual([row["measurement_revision"] for row in data["rows"]],
                             revisions)
            self.assertEqual(data["components"], ["k0", "k1"])
            self.assertEqual(data["channels"], ["c0", "c1"])
            # indices follow component order, not support insertion order
            self.assertEqual(data["rows"][0]["segments"][0]["constant"], [F(1), F(0)])
        finally:
            c.close()

    def test_invalid_requests_are_atomic(self):
        c = Catalog(":memory:")
        try:
            dictionary, revisions, _, parent = publish(c, [[1, 1], [0, 1]], [[1, 0]])
            before = c.counts()
            with self.assertRaises(KeyError):
                trace(c, "missing", [1, 1], 0, 1)
            # A dictionary revision is a record, not an artifact: unknown here.
            with self.assertRaises(KeyError):
                trace(c, dictionary, [1, 1], 0, 1)
            for values in ([], [1], [1, 1, 1], [None, 1], [True, 1], [0.1, 1],
                           ["nan", 1], ["1/0", 1], "not-a-sequence"):
                with self.subTest(values=values):
                    with self.assertRaises(ValueError):
                        trace(c, parent, values, 0, 1)
            for lo, hi in ((-1, 1), (2, 1), (None, 1), (0, None), (True, 1),
                           (0, 0.5), ("nan", 1)):
                with self.subTest(lo=lo, hi=hi):
                    with self.assertRaises(ValueError):
                        trace(c, parent, [1, 1], lo, hi)
            self.assertEqual(c.counts(), before)
            path_id = trace(c, parent, [1, 1], 0, 1)
            for bounds in ({"missing": (None, 1)}, {"k0": (2, 1)},
                           {"k0": (None,)}, {"k0": None}, {"k0": 1},
                           [(None, 1)], {"k0": (None, "nan")}, {"k0": (0.5, None)},
                           {"k0": (None, True)}):
                with self.subTest(bounds=bounds):
                    with self.assertRaises(ValueError):
                        acceptable(c, path_id, bounds)
            with self.assertRaises(ValueError):
                acceptable(c, parent, {})
            with self.assertRaises(KeyError):
                acceptable(c, "missing", {})
            with self.assertRaises(ValueError):
                trace(c, parent, [1, 1], 0, 1, fail_at="whenever")
            with self.assertRaises(ValueError):
                acceptable(c, path_id, {}, fail_at="whenever")
        finally:
            c.close()

    def test_publish_failure_stages_and_retry(self):
        c = Catalog(":memory:")
        try:
            _, _, _, parent = publish(c, [[1, 1], [0, 1]], [[1, 0], [0, 1]])
            for index, stage in enumerate(("after_payload", "after_index")):
                direction = [-1 - index, 1]
                before = c.counts()
                with self.assertRaises(RuntimeError):
                    trace(c, parent, direction, 0, 1, fail_at=stage)
                self.assertEqual(c.counts(), before)
                path_id = trace(c, parent, direction, 0, 1)
                self.assertEqual(trace(c, parent, direction, 0, 1), path_id)
                self.assertEqual(trace(c, parent, direction, "0", "1"), path_id)
                self.assertEqual(c.counts()["payloads"], c.counts()["artifacts"])
                for report_index, report_stage in enumerate(("after_payload", "after_index")):
                    bounds = {"k1": (None, F(3, 4) + report_index + index)}
                    before = c.counts()
                    with self.assertRaises(RuntimeError):
                        acceptable(c, path_id, bounds, fail_at=report_stage)
                    self.assertEqual(c.counts(), before)
                    report = acceptable(c, path_id, bounds)
                    self.assertEqual(acceptable(c, path_id, bounds), report)
                    self.assertEqual(c.counts()["payloads"], c.counts()["artifacts"])
        finally:
            c.close()

    def test_limits_semantics(self):
        c = Catalog(":memory:")
        try:
            _, revisions, path_id = run_path(c, [[1, 1], [0, 1]], [[1, 0], [0, 1]],
                                             [-1, 1], 0, 1)
            cases = [
                ({}, [[[F(0), F(1)]], [[F(0), F(1)]]]),
                ({"k1": (None, F(3, 4))}, [[[F(0), F(1)]], [[F(0), F(1)]]]),
                ({"k1": (F(1, 4), F(1, 4))}, [[[F(1, 4), F(1, 4)]], []]),
                ({"k1": (F(1, 2), F(1, 2))}, [[[F(1, 2), F(1)]], [[F(0), F(1)]]]),
                ({"k1": (F(3, 4), None)}, [[], []]),
                ({"k0": (0, F(1, 2)), "k1": (F(1, 4), None)},
                 [[[F(1, 4), F(1)]], [[F(0), F(1)]]]),
                ({"k0": (F(1, 2), F(1, 2)), "k1": (F(1, 4), F(1, 4))},
                 [[[F(1, 4), F(1, 4)]], []]),
            ]
            for bounds, expected in cases:
                report = acceptable(c, path_id, bounds)
                data = read(c, report)["data"]
                self.assertEqual(data["path_id"], path_id)
                self.assertEqual([row["sample_id"] for row in data["rows"]],
                                 ["s0", "s1"])
                self.assertEqual([row["measurement_revision"] for row in data["rows"]],
                                 revisions)
                self.assertEqual([row["intervals"] for row in data["rows"]], expected)
                self.assertEqual(data["bounds"],
                                 {key: [None if v is None else F(v) for v in pair]
                                  for key, pair in bounds.items()})
            # key order does not change identity; equal values do
            first = acceptable(c, path_id, {"k0": (None, None), "k1": (None, None)})
            second = acceptable(c, path_id, {"k1": [None, None], "k0": [None, None]})
            self.assertEqual(first, second)
            self.assertNotEqual(first, acceptable(c, path_id, {}))
        finally:
            c.close()

    def test_reads_do_not_solve_and_are_isolated(self):
        c = Catalog(":memory:")
        try:
            _, _, path_id = run_path(c, [[1, 1], [0, 1]], [[1, 0], [0, 1]],
                                     [-1, 1], 0, 1)
            report_id = acceptable(c, path_id, {"k1": (F(1, 4), None)})
            saved_path = read(c, path_id)
            saved_report = read(c, report_id)
            with mock.patch("spectral.response.paths.canonical_path",
                            side_effect=AssertionError("limit analysis solved a path")):
                with mock.patch("spectral.unmixing.solver.solve_canonical",
                                side_effect=AssertionError("read solved")):
                    self.assertEqual(read(c, path_id), saved_path)
                    self.assertEqual(read(c, report_id), saved_report)
                    again = acceptable(c, path_id, {"k1": (F(1, 4), None)})
                    self.assertEqual(again, report_id)
            altered = read(c, path_id)
            altered["data"]["rows"][0]["segments"][0]["constant"][0] = F(99)
            altered["data"]["rows"][0]["segments"][0]["slope"][1] = F(99)
            altered_report = read(c, report_id)
            altered_report["data"]["rows"][0]["intervals"][0][0] = F(99)
            self.assertEqual(read(c, path_id), saved_path)
            self.assertEqual(read(c, report_id), saved_report)
            # caller-side mutation after publication does not leak either
            bounds = {"k1": [F(1, 4), None]}
            identifier = acceptable(c, path_id, bounds)
            bounds["k1"][0] = F(3)
            self.assertEqual(read(c, identifier), read(c, report_id))
        finally:
            c.close()

    def test_correction_and_reopen_preserve_old_artifacts(self):
        from pathlib import Path
        import tempfile
        with tempfile.TemporaryDirectory() as directory:
            database = Path(directory) / "lab.sqlite"
            c = Catalog(database)
            try:
                dictionary = save_dictionary(c, ["c0", "c1"], ["k0", "k1"],
                                             ["signal", "signal"], ["mg/L", "mg/L"],
                                             [[1, 1], [0, 1]])
                batch = save_batch(c, ["c0", "c1"], ["signal", "signal"],
                                   [{"sample_id": "s0", "values": [1, 0]}])
                revisions = c.read_record(batch, "batch")["data"]["revision_ids"]
                old_parent = unmix(c, freeze(c, dictionary, revisions))
                old_path = trace(c, old_parent, [-1, 1], 0, 1)
                old_report = acceptable(c, old_path, {"k1": (F(1, 2), F(1, 2))})
                saved = (read(c, old_parent), read(c, old_path), read(c, old_report))
                newer = save_dictionary(c, ["c0", "c1"], ["k0", "k1"],
                                        ["signal", "signal"], ["mg/L", "mg/L"],
                                        [[2, 2], [0, 2]], previous=dictionary)
                corrected = append_measurement_revision(c, revisions[0], [1, 1])
                new_parent = unmix(c, freeze(c, newer, [corrected]))
                new_path = trace(c, new_parent, [-1, 1], 0, 1)
                self.assertNotEqual(new_path, old_path)
                self.assertEqual(read(c, old_path)["data"]["direction"], [F(-1), F(1)])
                self.assertEqual(read(c, old_path)["parents"],
                                 [{"type": "artifact", "id": old_parent}])
            finally:
                c.close()
            c = Catalog(database)
            try:
                c.initialize()
                self.assertEqual(read(c, old_parent), saved[0])
                self.assertEqual(read(c, old_path), saved[1])
                self.assertEqual(read(c, old_report), saved[2])
                self.assertEqual(trace(c, old_parent, [-1, 1], 0, 1), old_path)
                self.assertEqual(acceptable(c, old_path, {"k1": (F(1, 2), F(1, 2))}),
                                 old_report)
            finally:
                c.close()

    def test_path_does_not_follow_current_sources(self):
        c = Catalog(":memory:")
        try:
            dictionary = save_dictionary(c, ["c0"], ["k0"], ["signal"], ["mg/L"], [[1]])
            batch = save_batch(c, ["c0"], ["signal"],
                               [{"sample_id": "s", "values": [1]}])
            revisions = c.read_record(batch, "batch")["data"]["revision_ids"]
            parent = unmix(c, freeze(c, dictionary, revisions))
            path_id = trace(c, parent, [1], 0, 2)
            report = acceptable(c, path_id, {})
            # A later correction changes heads but not saved analyses.
            save_dictionary(c, ["c0"], ["k0"], ["signal"], ["mg/L"], [[5]],
                            previous=dictionary)
            append_measurement_revision(c, revisions[0], [9])
            self.assertEqual(read(c, path_id)["data"]["rows"][0]["segments"],
                             [{"lo": F(0), "hi": F(2), "constant": [F(1)],
                               "slope": [F(1)]}])
            self.assertEqual(read(c, report)["data"]["rows"][0]["intervals"],
                             [[F(0), F(2)]])
        finally:
            c.close()
