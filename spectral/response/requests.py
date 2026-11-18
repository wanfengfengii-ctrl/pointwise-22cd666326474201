from dataclasses import dataclass

from spectral.archive.codec import parse_rational, sequence
from spectral.archive.reading import read
from spectral.response.analysis import analyze_direction
from spectral.response.limits import acceptable
from spectral.response.paths import analyze_path, publish_path_result
from spectral.response.results import publish_direction_result


@dataclass(frozen=True)
class DirectionRequest:
    parent_id: str
    values: tuple


def direction(catalog, parent_id, values, fail_at=None):
    parent = read(catalog, parent_id)
    if parent["kind"] != "unmix":
        raise ValueError("direction requires an unmix parent")
    request = DirectionRequest(parent_id, tuple(parse_rational(value) for value in sequence(values, "direction")))
    if len(request.values) != len(parent["data"]["input"]["dictionary"]["channels"]):
        raise ValueError("direction dimension mismatch")
    responses = analyze_direction(parent["data"], request.values)
    return publish_direction_result(catalog, parent, request.values, responses, fail_at)


@dataclass(frozen=True)
class PathRequest:
    parent_id: str
    values: tuple
    lo: object
    hi: object


def _endpoint(value, label):
    try:
        return parse_rational(value)
    except ValueError as exc:
        raise ValueError("invalid path " + label) from exc


def trace(catalog, parent_id, d, lo, hi, fail_at=None):
    parent = read(catalog, parent_id)
    if parent["kind"] != "unmix":
        raise ValueError("trace requires an unmix parent")
    values = tuple(parse_rational(value) for value in sequence(d, "direction"))
    if len(values) != len(parent["data"]["input"]["dictionary"]["channels"]):
        raise ValueError("direction dimension mismatch")
    request = PathRequest(parent_id, values, _endpoint(lo, "lower"), _endpoint(hi, "upper"))
    if request.lo < 0:
        raise ValueError("path lower endpoint must be nonnegative")
    if request.lo > request.hi:
        raise ValueError("inverted path interval")
    rows = analyze_path(parent["data"], request.values, request.lo, request.hi)
    return publish_path_result(catalog, parent, request.values, request.lo,
                               request.hi, rows, fail_at)
