import itertools
from enum import IntEnum
from typing import Tuple

import numpy as np
import pytest
from numpy.testing import assert_equal

from threedigrid_builder.base import Array, replace, search, TooManyExist


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


class Records(Array[Record]):
    """Records docstring."""

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
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

    assert_equal(records.id_to_index(id, check_exists=True), expected)


def test_id_to_index_check_exists():
    records = Records(id=[1, 3, 5])

    with pytest.raises(KeyError) as e:
        records.id_to_index([1, 2, 3], check_exists=True)
        assert_equal(e.indices, [1])
        assert_equal(e.values, [2])

    with pytest.raises(KeyError) as e:
        records.id_to_index([1, 1, 2, 2], check_exists=True)
        assert_equal(e.indices, [2, 3])
        assert_equal(e.values, [2, 2])

    records.id_to_index([1, 1, 2, 2])  # no check, no error


@pytest.mark.parametrize(
    "id,expected",
    [
        (1, 3),
        ([2], [5]),
        ([1, 0], [3, 1]),
        ([1, 1, 1], [3, 3, 3]),
        (-9999, -9999),
        ([-9999], [-9999]),
        ([1, -9999], [3, -9999]),
    ],
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


def test_reorder():
    records = Records(id=[5, 7], number=[3, 2])
    records.reorder([1, 0])

    assert_equal(records.number, [2, 3])  # reordered
    assert_equal(records.id, [5, 7])  # kept the same


def test_reorder_by():
    records = Records(id=[5, 7], number=[3, 2])
    records.reorder_by("number")

    assert_equal(records.number, [2, 3])  # reordered
    assert_equal(records.id, [5, 7])  # kept the same


@pytest.mark.parametrize(
    "mask,expected",
    [
        ([True, True], [2, 3]),
        ([False, True], [3]),
        ([True, False], [2]),
        ([False, False], []),
        ([0, 1], [2, 3]),
        ([0], [2]),
        ([1], [3]),
        ([], []),
    ],
)
def test_indexing(mask, expected):
    records = Records(id=[5, 7], number=[2, 3])
    assert_equal(records[mask].number, records.number[mask])


def test_replace():
    actual = replace(
        np.array([0, 9, -150, 9999999]),
        {
            0: 1,
            9: 9999999,
            -150: -200,
            9999999: 9,
            1: 6,
        },
    )
    expected = np.array([1, 9999999, -200, 9])

    assert_equal(actual, expected)


def test_search():
    actual = search(
        [1, 4, 5, 9, 2],
        [4, 1, 2],
        assume_ordered=False,
        check_exists=True,
    )

    assert_equal(actual, [1, 0, 4])


def test_search_ordered():
    actual = search(
        [1, 2, 4, 5, 9],
        [4, 1, 2],
        assume_ordered=True,
        check_exists=True,
    )

    assert_equal(actual, [2, 0, 1])


def test_search_missing():
    actual = search(
        [1, 4, 5, 9, 2],
        [4, 1, 10],
        assume_ordered=False,
        check_exists=False,
    )

    assert_equal(actual, [1, 0, -9999])


def test_search_bool_mask():
    actual = search(
        [1, 4, 2, 9, 2],
        [4, 1, 2],
        mask=[True, True, False, False, True],
        assume_ordered=False,
        check_exists=True,
    )

    assert_equal(actual, [1, 0, 4])


def test_search_int_mask():
    actual = search(
        [1, 4, 2, 9, 2],
        [4, 1, 2],
        mask=[0, 1, 4],
        assume_ordered=False,
        check_exists=True,
    )

    assert_equal(actual, [1, 0, 4])


def test_search_int_mask_ordered():
    actual = search(
        [1, 4, 2, 9, 2],
        [4, 1, 2],
        mask=[0, 4, 1],  # this orders a
        assume_ordered=True,
        check_exists=True,
    )

    assert_equal(actual, [1, 0, 4])


def test_search_int_mask_missing():
    actual = search(
        [1, 4, 2, 9, 2],
        [4, 1, 9],
        mask=[0, 1, 4],
        assume_ordered=False,
        check_exists=False,
    )

    assert_equal(actual, [1, 0, -9999])


def test_search_check_exists():
    with pytest.raises(KeyError) as e:
        search(
            [1, 2, 4, 5, 9],
            [4, 6, 1, 2, 3],
            assume_ordered=False,
            check_exists=True,
        )
        assert_equal(e.indices, [1, 4])
        assert_equal(e.values, [6, 3])


def test_search_check_exists_mask():
    with pytest.raises(KeyError) as e:
        search(
            [1, 2, 4, 5, 9],
            [4, 1],
            mask=[0, 1, 3],
            assume_ordered=False,
            check_exists=True,
        )
        assert_equal(e.indices, [0])
        assert_equal(e.values, [4])


def test_search_empty():
    actual = search(
        [],
        [4, 1, 2],
        assume_ordered=False,
        check_exists=False,
    )

    assert_equal(actual, [-9999, -9999, -9999])


def test_search_empty_check_exists():
    with pytest.raises(KeyError):
        search(
            [],
            [4, 1, 2],
            assume_ordered=False,
            check_exists=True,
        )


def test_search_empty_mask():
    actual = search(
        [0],
        [4, 1, 2],
        mask=[False],
        assume_ordered=False,
        check_exists=False,
    )

    assert_equal(actual, [-9999, -9999, -9999])


def test_search_empty_mask_check_exists():
    with pytest.raises(KeyError):
        search(
            [0],
            [4, 1, 2],
            mask=[False],
            assume_ordered=False,
            check_exists=True,
        )


def test_split_in_two():
    records = Records(id=[0, 1, 2, 3, 4, 5])

    actual_1, actual_2 = records.split_in_two([0, 1, 1, 2, 3, 2])
    assert_equal(actual_1, [0, 1, 3, 4])
    assert_equal(actual_2, [2, 5])


def test_split_in_two_empty():
    records = Records(id=[])

    actual_1, actual_2 = records.split_in_two([])
    assert len(actual_1) == 0
    assert len(actual_2) == 0


def test_split_in_two_unique():
    records = Records(id=[0, 1, 2])

    actual_1, actual_2 = records.split_in_two([1, 2, 7])
    assert_equal(actual_1, [0, 1, 2])
    assert_equal(actual_2, [])


def test_split_in_two_err():
    records = Records(id=[0, 1, 2, 3])

    with pytest.raises(TooManyExist) as e:
        records.split_in_two([5, 2, 5, 5])

    assert_equal(e.value.values, [5])


class RecordsWithScalarArg(Array[Record]):
    """Records docstring."""

    scalars = ("foo",)


def test_scalar_arg_init():
    records = RecordsWithScalarArg(id=[], foo="bar")
    assert records.foo == "bar"


def test_scalar_arg_init_typeerror():
    with pytest.raises(TypeError):
        RecordsWithScalarArg(id=[])


def test_scalar_arg_getitem():
    records = RecordsWithScalarArg(id=[1, 2], foo="bar")
    actual = records[:1]
    assert len(actual) == 1
    assert actual.foo == "bar"


def test_scalar_arg_add():
    a = RecordsWithScalarArg(id=[1, 2], foo="bar")
    b = RecordsWithScalarArg(id=[3, 4], foo="baz")
    actual = a + b
    assert len(actual) == 4
    assert actual.foo == "bar"


@pytest.mark.parametrize(
    "name,index,expected",
    [
        ("number", [1, 2, -9999], [2.3, 3.4, np.nan]),
        ("name", [1, 2, -9999], ["b", "c", None]),
        ("animal", [1, 2, -9999], [2, 3, -9999]),
        ("number", 2, 3.4),
    ],
)
def test_take(name, index, expected):
    records = Records(
        id=range(3), name=["a", "b", "c"], animal=[1, 2, 3], number=[1.2, 2.3, 3.4]
    )

    assert_equal(records.take(name, index), expected)


def test_empty_column():
    records = Records(id=[1, 2, 3])
    assert "number" not in records.__dict__
    arr = records.number
    assert records.__dict__["number"] is arr
    assert_equal(arr, np.nan)
