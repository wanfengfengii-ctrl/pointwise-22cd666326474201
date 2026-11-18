from dataclasses import dataclass, replace
from fractions import Fraction

from spectral.dictionary.operators import DictionaryOperator
from spectral.unmixing.faces import enumerate_face_candidates, solve_face
from spectral.unmixing.linear import dot


@dataclass(frozen=True)
class ResponseBranch:
    constant: tuple
    slope: tuple
    residual_polynomial: tuple
    norm_polynomial: tuple
    right_limit: Fraction


def _square_polynomial(constant, slope):
    return (dot(constant, constant), 2 * dot(constant, slope), dot(slope, slope))


def enumerate_right_branches(matrix, values, direction):
    operator = DictionaryOperator(matrix)
    if len(direction) != operator.m:
        raise ValueError("direction dimension mismatch")
    branches = []
    exclusion_limit = Fraction(1)
    for candidate in enumerate_face_candidates(matrix, values):
        u = candidate.coefficients
        v = solve_face(matrix, direction, candidate.support).coefficients
        if any(a < 0 or (a == 0 and b < 0) for a, b in zip(u, v)):
            permanent = any((a <= 0 and b < 0) or (a < 0 and b == 0)
                            for a, b in zip(u, v))
            if not permanent:
                entrance = max(-a / b for a, b in zip(u, v) if a < 0)
                exclusion_limit = min(exclusion_limit, entrance)
            continue
        limit = min([Fraction(1)] + [-a / b for a, b in zip(u, v) if a > 0 and b < 0])
        residual = tuple(a - b for a, b in zip(operator.apply(u), values))
        residual_slope = tuple(a - b for a, b in zip(operator.apply(v), direction))
        branches.append(ResponseBranch(u, v, _square_polynomial(residual, residual_slope),
                                       _square_polynomial(u, v), limit))
    for branch in branches:
        yield replace(branch, right_limit=min(branch.right_limit, exclusion_limit))


def _sign_interval(polynomial):
    for i, coefficient in enumerate(polynomial):
        if coefficient:
            tail = sum(abs(x) for x in polynomial[i + 1:])
            return min(Fraction(1), abs(coefficient) / (2 * tail + 1))
    return Fraction(1)


def choose_right_branch(branches):
    branches = tuple(branches)
    if not branches:
        raise ValueError("no feasible right branch")
    winner = min(branches, key=lambda b: (b.residual_polynomial, b.norm_polynomial))
    limit = min(branch.right_limit for branch in branches)
    for branch in branches:
        difference = tuple(a - b for a, b in
                           zip(branch.residual_polynomial, winner.residual_polynomial))
        if not any(difference):
            difference = tuple(a - b for a, b in
                               zip(branch.norm_polynomial, winner.norm_polynomial))
        limit = min(limit, _sign_interval(difference))
    return replace(winner, right_limit=limit)
