"""Component bounds applied to an already saved complete path.

Nothing here solves or unmixes.  The saved segments are treated as immutable
input; only their affine coefficients and the component coordinates decide the
acceptable parameter set.  Every endpoint is an exact Fraction and the result
is a list of maximal closed intervals in ascending order.
"""
from spectral.archive.codec import encode_payload, parse_rational
from spectral.archive.publication import publish_artifact
from spectral.archive.reading import read


def _parse_endpoint(value, label):
    if value is None:
        return None
    try:
        return parse_rational(value)
    except ValueError as exc:
        raise ValueError("invalid " + label + " endpoint") from exc


def parse_bounds(components, bounds):
    if not isinstance(bounds, dict):
        raise ValueError("bounds must be a mapping")
    known = set(components)
    for key in bounds:
        if not isinstance(key, str) or key not in known:
            raise ValueError("unknown component bound")
    parsed = {}
    for key, pair in bounds.items():
        if not isinstance(pair, (list, tuple)) or len(pair) != 2:
            raise ValueError("bound must be a lower/upper pair")
        lower = _parse_endpoint(pair[0], "lower")
        upper = _parse_endpoint(pair[1], "upper")
        if lower is not None and upper is not None and lower > upper:
            raise ValueError("inverted bound")
        parsed[key] = (lower, upper)
    return parsed


def _clip_affine(lo, hi, base, rate, lower, upper):
    """Closed t-subinterval of [lo,hi] with lower <= base+rate*t <= upper."""
    if rate == 0:
        if (lower is not None and base < lower) or (
                upper is not None and base > upper):
            return None
        return (lo, hi)
    if lower is not None:
        bound = (lower - base) / rate
        if rate > 0:
            lo = max(lo, bound)
        else:
            hi = min(hi, bound)
    if upper is not None:
        bound = (upper - base) / rate
        if rate > 0:
            hi = min(hi, bound)
        else:
            lo = max(lo, bound)
    if lo > hi:
        return None
    return (lo, hi)


def _row_intervals(components, segments, bounds):
    positions = {name: index for index, name in enumerate(components)}
    pieces = []
    for segment in segments:
        lo, hi = segment["lo"], segment["hi"]
        for key, (lower, upper) in bounds.items():
            index = positions[key]
            clipped = _clip_affine(lo, hi, segment["constant"][index],
                                   segment["slope"][index], lower, upper)
            if clipped is None:
                break
            lo, hi = clipped
        else:
            pieces.append((lo, hi))
    pieces.sort()
    merged = []
    for interval in pieces:
        if merged and interval[0] <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], interval[1]))
        else:
            merged.append(interval)
    return [[lo, hi] for lo, hi in merged]


def build_acceptable_payload(path_id, path_data, bounds):
    components = path_data["components"]
    rows = []
    for path_row in path_data["rows"]:
        rows.append({"sample_id": path_row["sample_id"],
                     "measurement_revision": path_row["measurement_revision"],
                     "intervals": _row_intervals(components,
                                                 path_row["segments"], bounds)})
    ordered = {key: [lower, upper] for key, (lower, upper) in sorted(bounds.items())}
    return {"path_id": path_id, "bounds": ordered, "rows": rows}


def publish_acceptable_result(catalog, path, bounds, fail_at=None):
    payload = build_acceptable_payload(path["id"], path["data"], bounds)
    ordered = {key: list(pair) for key, pair in payload["bounds"].items()}
    key = "acceptable:" + path["id"] + ":" + encode_payload(
        [[name, pair] for name, pair in ordered.items()])
    return publish_artifact(catalog, "acceptable", key, payload,
                            [{"type": "artifact", "id": path["id"]}], fail_at)


def acceptable(catalog, path_id, bounds, fail_at=None):
    path = read(catalog, path_id)
    if path["kind"] != "path":
        raise ValueError("acceptable requires a path parent")
    parsed = parse_bounds(path["data"]["components"], bounds)
    return publish_acceptable_result(catalog, path, parsed, fail_at)
