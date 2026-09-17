"""The AR table's span test, against a direct linear-algebra reference.

`Table` keeps an incremental echelon basis and memoises results, so it is worth checking
against a straightforward rank computation on the same rows.
"""

import random

import numpy as np
import pytest

from ar import Table
from relations import RelationNode


def reference_is_spanned(rows, query, tol=1e-10):
    """Ground truth: is `query` in the row space of `rows`?"""
    if not any(query):
        return True  # the zero vector is in every span, including the empty one
    if not rows:
        return False
    matrix = np.array([list(r) + [0] * (len(query) - len(r)) for r in rows], dtype=float)
    vector = np.array(query, dtype=float)
    rank = np.linalg.matrix_rank(matrix, tol=tol)
    rank_aug = np.linalg.matrix_rank(np.vstack([matrix, vector]), tol=tol)
    return rank == rank_aug


def test_matches_reference_on_random_tables():
    rng = random.Random(7)
    mismatches = 0
    spanned_seen = 0
    for _ in range(200):
        width = rng.randint(3, 8)
        table = Table([f"c{i}" for i in range(width)])
        stored = []
        for step in range(rng.randint(2, 14)):
            query = [rng.choice([0, 0, 1, -1]) for _ in range(width)]
            got, _ = table.is_spanned(query)
            want = reference_is_spanned(stored, query)
            spanned_seen += int(want)
            mismatches += int(got != want)

            row = [rng.choice([0, 0, 0, 1, -1, 2]) for _ in range(width)]
            if stored and rng.random() < 0.3:
                # sometimes feed a row that is a combination of existing ones
                row = list(np.sum([rng.choice([-1, 1, 2]) * np.array(b)
                                   for b in rng.sample(stored, min(2, len(stored)))],
                                  axis=0))
            table.add_row(list(row), RelationNode(f"cong", step, ("k", step)))
            if not reference_is_spanned(stored, row):
                stored.append(list(row))
    assert mismatches == 0
    assert spanned_seen > 0, "test never exercised the spanned branch"


def test_dependency_certificate_reproduces_the_query():
    table = Table(["a", "b", "c"])
    r1 = RelationNode("cong", 1, ("k", 1))
    r2 = RelationNode("cong", 2, ("k", 2))
    table.add_row([1, -1, 0], r1)
    table.add_row([0, 1, -1], r2)

    spanned, used = table.is_spanned([1, 0, -1])
    assert spanned
    assert used == {r1, r2}, "both rows are needed to reach a - c"

    spanned, used = table.is_spanned([1, -1, 0])
    assert spanned and used == {r1}


def test_only_independent_rows_are_stored():
    table = Table(["a", "b"])
    r = RelationNode("cong", 1, ("k", 1))
    table.add_row([1, -1], r)
    table.add_row([2, -2], r)   # a multiple of the first: adds nothing
    table.add_row([0, 0], r)    # vacuous
    assert table._n == 1


def test_adding_a_column_preserves_earlier_answers():
    table = Table(["a", "b"])
    r = RelationNode("cong", 1, ("k", 1))
    table.add_row([1, -1], r)
    assert table.is_spanned([1, -1])[0]
    table.add_col("c")
    assert table.is_spanned([1, -1, 0])[0]
    assert not table.is_spanned([0, 1, -1])[0]


def test_memo_is_invalidated_when_a_row_arrives():
    table = Table(["a", "b", "c"])
    r1 = RelationNode("cong", 1, ("k", 1))
    r2 = RelationNode("cong", 2, ("k", 2))
    table.add_row([1, -1, 0], r1)
    assert not table.is_spanned([0, 1, -1])[0]   # cached as a miss
    table.add_row([0, 1, -1], r2)
    assert table.is_spanned([0, 1, -1])[0], "stale miss was served from the cache"
