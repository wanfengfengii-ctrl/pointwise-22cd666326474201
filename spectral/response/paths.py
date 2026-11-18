"""Complete piecewise-affine canonical paths over a closed rational interval.

Only the frozen snapshot stored in the parent unmix artifact is used.  Every
support is solved once as an affine family in the global parameter ``t``:
x(t) = constant + t * slope.  Candidate breakpoints are the exact rational
values where a support enters or leaves the nonnegative cone, or where two
affine expressions coincide.  Because the canonical solution is unique and
continuous in t, a change of formula can only happen at such a rational point;
each resulting piece is nevertheless proved analytically on the whole closed
piece against every competing support.  Sampling is never used as a proof.
"""
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations
from math import isqrt

from spectral.archive.codec import encode_payload
from spectral.archive.publication import publish_artifact
from spectral.dictionary.operators import DictionaryOperator
from spectral.unmixing.linear import dot, minimum_norm_solution
from spectral.unmixing.solver import solve_canonical


@dataclass(frozen=True)
class AffineFace:
    support: tuple
    constant: tuple
    slope: tuple
    feasible_lo: object
    feasible_hi: object
    residual_polynomial: tuple
    norm_polynomial: tuple


def _face(operator, y0, y1, support):
    """Solve x(t) = constant + t * slope on one support, minimum norm."""
    n = operator.n
    columns = [tuple(row[j] for row in operator.matrix) for j in support]
    count = len(columns)
    gram = [[dot(left, right) for right in columns] for left in columns]
    target0 = [dot(column, y0) for column in columns]
    target1 = [dot(column, y1) for column in columns]
    local0 = minimum_norm_solution(gram, target0, count) if count else ()
    local1 = minimum_norm_solution(gram, target1, count) if count else ()
    constant = [Fraction() for _ in range(n)]
    slope = [Fraction() for _ in range(n)]
    for position, value in zip(support, local0):
        constant[position] = value
    for position, value in zip(support, local1):
        slope[position] = value
    feasible_lo = None
    feasible_hi = None
    for position in support:
        base, rate = constant[position], slope[position]
        if rate == 0:
            if base < 0:
                return None
        else:
            bound = -base / rate
            if rate > 0:
                feasible_lo = bound if feasible_lo is None else max(feasible_lo, bound)
            else:
                feasible_hi = bound if feasible_hi is None else min(feasible_hi, bound)
    if (feasible_lo is not None and feasible_hi is not None
            and feasible_lo > feasible_hi):
        return None
    residual0 = tuple(a - b for a, b in zip(operator.apply(constant), y0))
    residual1 = tuple(a - b for a, b in zip(operator.apply(slope), y1))
    return AffineFace(
        tuple(support), tuple(constant), tuple(slope), feasible_lo, feasible_hi,
        (dot(residual0, residual0), 2 * dot(residual0, residual1),
         dot(residual1, residual1)),
        (dot(constant, constant), 2 * dot(constant, slope), dot(slope, slope)))


def enumerate_affine_faces(matrix, values, direction):
    operator = DictionaryOperator(matrix)
    y0 = tuple(Fraction(value) for value in values)
    y1 = tuple(Fraction(value) for value in direction)
    for size in range(operator.n + 1):
        for support in combinations(range(operator.n), size):
            face = _face(operator, y0, y1, support)
            if face is not None:
                yield face


def _value(polynomial, t):
    return polynomial[0] + t * (polynomial[1] + t * polynomial[2])


def _nonnegative_on(polynomial, lo, hi):
    """Exact sign proof for a rational quadratic on a closed interval."""
    a0, a1, a2 = polynomial
    smallest = min(_value(polynomial, lo), _value(polynomial, hi))
    if a2 > 0:
        vertex = -a1 / (2 * a2)
        if lo < vertex < hi:
            smallest = min(smallest, _value(polynomial, vertex))
    return smallest >= 0


def _touches(face, lo, hi):
    return ((face.feasible_lo is None or face.feasible_lo <= hi)
            and (face.feasible_hi is None or face.feasible_hi >= lo))


def _overlaps_open(face, lo, hi):
    return ((face.feasible_lo is None or face.feasible_lo < hi)
            and (face.feasible_hi is None or face.feasible_hi > lo))


def _covers(face, lo, hi):
    return ((face.feasible_lo is None or face.feasible_lo <= lo)
            and (face.feasible_hi is None or face.feasible_hi >= hi))


def _feasible_at(face, t):
    return ((face.feasible_lo is None or face.feasible_lo <= t)
            and (face.feasible_hi is None or face.feasible_hi >= t))


def _coincidence(left, right):
    """Return the unique t with left(t) == right(t), else None.

    None also covers identical affine pieces and parallel pieces that never
    meet; neither requires a breakpoint.
    """
    point = None
    for a0, a1, b0, b1 in zip(left.constant, left.slope,
                              right.constant, right.slope):
        delta0, delta1 = a0 - b0, a1 - b1
        if delta1:
            candidate = -delta0 / delta1
            if point is None:
                point = candidate
            elif point != candidate:
                return None
        elif delta0:
            return None
    return point


def _breakpoints(faces, lo, hi):
    points = {lo, hi}
    for face in faces:
        for position in face.support:
            rate = face.slope[position]
            if rate:
                root = -face.constant[position] / rate
                if lo < root < hi:
                    points.add(root)
    for index, left in enumerate(faces):
        for right in faces[index + 1:]:
            if not (_overlaps_open(left, lo, hi) and _overlaps_open(right, lo, hi)):
                continue
            point = _coincidence(left, right)
            if point is not None and lo < point < hi:
                points.add(point)
    return sorted(points)


def _zeros_in(polynomial, lo, hi):
    """Exact zeros of a rational quadratic on a closed interval."""
    a0, a1, a2 = polynomial
    if a2 == 0:
        if a1 == 0:
            return [lo, hi] if a0 == 0 else []
        root = -a0 / a1
        return [root] if lo <= root <= hi else []
    discriminant = a1 * a1 - 4 * a2 * a0
    if discriminant < 0:
        return []
    numerator, denominator = discriminant.numerator, discriminant.denominator
    sn, sd = isqrt(numerator), isqrt(denominator)
    if sn * sn != numerator or sd * sd != denominator:
        # A nonnegative quadratic cannot have two irrational roots and an
        # irrational double root would force a rational vertex; no zero exists.
        return []
    square = Fraction(sn, sd)
    roots = [(-a1 - square) / (2 * a2), (-a1 + square) / (2 * a2)]
    return sorted({root for root in roots if lo <= root <= hi})


def _beats(winner, other, lo, hi):
    """Prove lexicographic residual-then-norm ordering on one closed span.

    The span passed here is always the competitor's own feasible overlap with
    the piece, so an infeasible competitor is never asked to be optimal.
    """
    difference = tuple(a - b for a, b in
                       zip(other.residual_polynomial, winner.residual_polynomial))
    if any(difference):
        if not _nonnegative_on(difference, lo, hi):
            return False
        # Residual ties only matter at exact zeros; the norm must win there.
        norm_difference = tuple(
            a - b for a, b in
            zip(other.norm_polynomial, winner.norm_polynomial))
        return all(_value(norm_difference, t) >= 0
                   for t in _zeros_in(difference, lo, hi))
    difference = tuple(a - b for a, b in
                       zip(other.norm_polynomial, winner.norm_polynomial))
    return _nonnegative_on(difference, lo, hi)


def _feasible_span(face, lo, hi):
    left = lo if face.feasible_lo is None else max(lo, face.feasible_lo)
    right = hi if face.feasible_hi is None else min(hi, face.feasible_hi)
    return (left, right) if left <= right else None


def _segment(lo, hi, constant, slope):
    return {"lo": Fraction(lo), "hi": Fraction(hi),
            "constant": [Fraction(value) for value in constant],
            "slope": [Fraction(value) for value in slope]}


def canonical_path(matrix, values, direction, lo, hi):
    operator = DictionaryOperator(matrix)
    y0 = tuple(Fraction(value) for value in values)
    y1 = tuple(Fraction(value) for value in direction)
    lo, hi = Fraction(lo), Fraction(hi)
    if lo > hi:
        raise ValueError("inverted path interval")
    if lo == hi:
        solution = solve_canonical(matrix, [a + lo * b for a, b in zip(y0, y1)])
        return [_segment(lo, hi, solution.coefficients,
                         [Fraction() for _ in range(operator.n)])]
    faces = [face for face in enumerate_affine_faces(matrix, y0, y1)
             if _touches(face, lo, hi)]
    if not faces:
        raise ValueError("no feasible affine face on the path interval")
    points = _breakpoints(faces, lo, hi)
    pieces = []
    for left, right in zip(points, points[1:]):
        middle = (left + right) / 2
        feasible = [face for face in faces if _feasible_at(face, middle)]
        winner = min(feasible, key=lambda face: (
            _value(face.residual_polynomial, middle),
            _value(face.norm_polynomial, middle)))
        if not _covers(winner, left, right):
            raise ValueError("path piece leaves its support cone")
        for other in faces:
            if other is winner:
                continue
            span = _feasible_span(other, left, right)
            if span is None:
                continue
            if not _beats(winner, other, span[0], span[1]):
                raise ValueError("path piece is not lexicographically minimal")
        pieces.append(_segment(left, right, winner.constant, winner.slope))
    merged = []
    for piece in pieces:
        if (merged and merged[-1]["constant"] == piece["constant"]
                and merged[-1]["slope"] == piece["slope"]):
            merged[-1]["hi"] = piece["hi"]
        else:
            merged.append(piece)
    _check_continuity(merged, lo, hi)
    return merged


def _check_continuity(segments, lo, hi):
    if not segments or segments[0]["lo"] != lo or segments[-1]["hi"] != hi:
        raise ValueError("path does not cover the requested interval")
    for left, right in zip(segments, segments[1:]):
        if left["hi"] != right["lo"]:
            raise ValueError("path has a gap or overlap")
        if (left["constant"], left["slope"]) == (right["constant"], right["slope"]):
            raise ValueError("adjacent identical path segments were not merged")
        left_value = [a + left["hi"] * b
                      for a, b in zip(left["constant"], left["slope"])]
        right_value = [a + right["lo"] * b
                       for a, b in zip(right["constant"], right["slope"])]
        if left_value != right_value:
            raise ValueError("path is discontinuous at a shared breakpoint")


def analyze_path(parent_payload, direction, lo, hi):
    snapshot = parent_payload["input"]
    dictionary = snapshot["dictionary"]
    matrix = dictionary["matrix"]
    operator = DictionaryOperator(matrix)
    if len(direction) != operator.m:
        raise ValueError("direction dimension mismatch")
    if len(parent_payload["rows"]) != len(snapshot["measurements"]):
        raise ValueError("incomplete parent rows")
    rows = []
    for measurement, saved in zip(snapshot["measurements"],
                                  parent_payload["rows"]):
        if (saved["sample_id"] != measurement["sample_id"]
                or saved["measurement_revision"] != measurement["revision_id"]):
            raise ValueError("parent record identity mismatch")
        segments = canonical_path(matrix, measurement["values"], direction, lo, hi)
        rows.append({"sample_id": measurement["sample_id"],
                     "measurement_revision": measurement["revision_id"],
                     "segments": segments})
    return rows


def build_path_payload(parent, direction, lo, hi, rows):
    data = parent["data"]
    dictionary = data["input"]["dictionary"]
    if len(rows) != len(data["rows"]):
        raise ValueError("incomplete path results")
    for saved, result in zip(data["rows"], rows):
        if (saved["sample_id"] != result["sample_id"]
                or saved["measurement_revision"] != result["measurement_revision"]):
            raise ValueError("path record identity mismatch")
    return {"parent_id": parent["id"], "cohort_id": data["cohort_id"],
            "channels": dictionary["channels"],
            "channel_units": dictionary["channel_units"],
            "components": dictionary["components"],
            "component_units": dictionary["component_units"],
            "direction": tuple(direction), "lo": lo, "hi": hi, "rows": rows}


def publish_path_result(catalog, parent, direction, lo, hi, rows, fail_at=None):
    payload = build_path_payload(parent, direction, lo, hi, rows)
    key = "path:" + parent["id"] + ":" + encode_payload(
        {"direction": tuple(direction), "lo": lo, "hi": hi})
    return publish_artifact(catalog, "path", key, payload,
                            [{"type": "artifact", "id": parent["id"]}], fail_at)
