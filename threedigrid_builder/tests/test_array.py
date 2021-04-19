from enum import IntEnum
from numpy.testing import assert_equal
from threedigrid_builder.base import array_of
from typing import Tuple

import itertools
import numpy as np
import pytest


class Animal(IntEnum):
    ANT = 1
    BEE = 2
    CAT = 3
    DOG = 4


class Record:
    id: int
    number: float
    yesno: bool
    name: str
    animal: Animal
    xyz: Tuple[float, float, float]


@array_of(Record)
class Records:
    """Records docstring."""

    def __init__(self):
        self.length_for_init_test = len(self.id)

    def does_it_work(self):
        return "yes"


def test_arrays_from_lists():
    records = Records(
        id=[0, 1],
        number=[5, 7],
        yesno=[True, False],
        name=["a", "b"],
        animal=[Animal.ANT, -9999],
        xyz=[[1.0, 1.0, 1.0], [2.0, 2.0, 2.0]],
    )

    assert_equal(records.id, np.array([0, 1], dtype=np.int32))
    assert_equal(records.number, np.array([5.0, 7.0], dtype=np.float64))
    assert_equal(records.yesno, np.array([True, False], dtype=bool))
    assert_equal(records.name, np.array(["a", "b"], dtype=object))
    assert_equal(records.animal, np.array([1, -9999], dtype=np.int32))
    assert_equal(
        records.xyz,
        np.array([[1.0, 1.0, 1.0], [2.0, 2.0, 2.0]], dtype=np.float64, order="F"),
    )


def test_arrays_from_arrays():
    records = Records(
        id=np.array([0, 1], dtype=np.uint16),
        number=np.array([5, 7], dtype=int),
        yesno=np.array([1, 0], dtype=int),
    )

    assert_equal(records.id, np.array([0, 1], dtype=np.int32))
    assert_equal(records.number, np.array([5.0, 7.0], dtype=np.float64))
    assert_equal(records.yesno, np.array([True, False], dtype=bool))


def test_empty_values():
    records = Records(id=[0, 1])

    assert_equal(records.number, np.array([np.nan, np.nan], dtype=np.float64))
    assert_equal(records.yesno, np.array([False, False], dtype=bool))
    assert_equal(records.name, np.array([None, None], dtype=object))
    assert_equal(records.animal, np.array([-9999, -9999], dtype=np.int32))
    assert_equal(records.xyz, np.full((2, 3), np.nan, dtype=np.float64))


def test_broadcast_scalar_values():
    records = Records(id=[0, 1], number=5.2)

    assert_equal(records.number, np.array([5.2, 5.2], dtype=np.float64))


def test_id_from_generator():
    records = Records(id=(x for x in [5, 7]))
    assert_equal(records.id, [5, 7])


def test_id_from_islice():
    records = Records(id=itertools.islice([5, 7, 8], 2))
    assert_equal(records.id, [5, 7])


def test_id_required():
    with pytest.raises(TypeError):
        Records(number=[5, 7])


def test_id_dimensionality():
    with pytest.raises(ValueError):
        Records(id=7)


def test_id_sorted():
    with pytest.raises(ValueError):
        Records(id=[7, 5])


def test_unqual_lengths():
    with pytest.raises(ValueError):
        Records(id=[1], number=[5, 7])


@pytest.mark.parametrize(
    "xyz", [[[1, 1, 1], [2, 2, 2]], np.array([[1.0, 1.0, 1.0], [2.0, 2.0, 2.0]]), None]
)
def test_fortran_order(xyz):
    records = Records(id=[0, 1], xyz=xyz)

    assert records.xyz.flags["F_CONTIGUOUS"]
    assert records.xyz.shape == (2, 3)


def test_enum_invalid():
    with pytest.raises(ValueError):
        Records(id=range(5), animal=[1, 2, 3, 4, 5])


def test_python_attrs():
    assert Records.__name__ == Records.__qualname__ == "Records"
    assert Records.__doc__ == "Records docstring."


def test_init():
    assert Records(id=[1, 2]).length_for_init_test == 2


def test_len():
    assert len(Records(id=[1, 2])) == 2


def test_repr():
    assert repr(Records(id=[2, 3])) == "<Records object, Record array of length 2>"


@pytest.mark.parametrize(
    "id,expected", [(1, 0), ([5], [2]), ([5, 3], [2, 1]), ([1, 1, 1], [0, 0, 0])]
)
def test_id_to_index(id, expected):
    records = Records(id=[1, 3, 5])

    assert_equal(records.id_to_index(id), expected)


@pytest.mark.parametrize(
    "id,expected", [(1, 3), ([2], [5]), ([1, 0], [3, 1]), ([1, 1, 1], [3, 3, 3])]
)
def test_index_to_id(id, expected):
    records = Records(id=[1, 3, 5])

    assert_equal(records.index_to_id(id), expected)


def test_methods():
    assert Records(id=[1]).does_it_work() == "yes"


def test_concatenate():
    a = Records(id=[1], xyz=[[0, 0, 0]])
    b = Records(id=[2], xyz=[[1, 1, 1]])
    actual = a + b
    assert len(actual) == 2
    assert_equal(actual.id, [1, 2])
    assert_equal(actual.xyz, [[0, 0, 0], [1, 1, 1]])


def test_concatenate_inplace():
    a = Records(id=[1], xyz=[[0, 0, 0]])
    b = Records(id=[2], xyz=[[1, 1, 1]])
    a += b
    assert len(a) == 2
    assert_equal(a.id, [1, 2])
    assert_equal(a.xyz, [[0, 0, 0], [1, 1, 1]])


def test_reorder_by():
    records = Records(id=[5, 7], number=[2, 3])
    records.reorder_by([1, 0])

    assert_equal(records.number, [3, 2])  # reordered
    assert_equal(records.id, [5, 7])  # kept the same
