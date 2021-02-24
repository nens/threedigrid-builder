from threedigrid_builder.base import array_of
from numpy.testing import assert_equal
import numpy as np
import pytest
from typing import Tuple
from enum import IntEnum


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
    def does_it_work(self):
        return "yes"


def test_arrays_from_lists():
    records = Records(
        id=[0, 1],
        number=[5, 7],
        yesno=[True, False],
        name=["a", "b"],
        animal=[Animal.ANT, Animal.DOG],
        xyz=[[1.0, 1.0, 1.0], [2.0, 2.0, 2.0]],
    )

    assert_equal(records.id, np.array([0, 1], dtype=np.int32))
    assert_equal(records.number, np.array([5.0, 7.0], dtype=np.float64))
    assert_equal(records.yesno, np.array([True, False], dtype=bool))
    assert_equal(records.name, np.array(["a", "b"], dtype=object))
    assert_equal(records.animal, np.array([1, 4], dtype=np.int32))
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


def test_none():
    records = Records(id=[0, 1])

    assert records.number is None
    assert records.yesno is None
    assert records.name is None


def test_id_required():
    with pytest.raises(TypeError):
        Records(number=[5, 7])


def test_unqual_lengths():
    with pytest.raises(ValueError):
        Records(id=[1], number=[5, 7])


def test_fortran_order():
    records = Records(id=[0, 1], xyz=np.array([[1.0, 1.0, 1.0], [2.0, 2.0, 2.0]]))

    assert records.xyz.flags["F_CONTIGUOUS"]
    assert records.xyz.shape == (2, 3)


def test_python_attrs():
    assert Records.__name__ == Records.__qualname__ == "Records"
    assert Records.__doc__ == "Records docstring."
    assert Records.__module__ == "threedigrid_builder.base.tests.test_array"


def test_len():
    assert len(Records(id=[2, 1])) == 2


def test_repr():
    assert repr(Records(id=[2, 3])) == "<Records object, Record array of length 2>"


def test_methods():
    assert Records(id=[1]).does_it_work() == "yes"
