import typing
from enum import IntEnum
from typing import _GenericAlias, Generic, TypeVar

import numpy as np

__all__ = [
    "is_tuple_type",
    "is_int_enum",
    "unpack_optional_type",
    "replace",
    "search",
    "Array",
    "TooManyExist",
]


def is_tuple_type(_type):
    """Test if the type is a indexed typing.Tuple type.

    Examples:
        is_tuple_type(Tuple[int, int]) == True
    """
    return isinstance(_type, _GenericAlias) and _type.__origin__ is tuple


def is_union_type(_type):
    """Test if the type is an typing.Union type.

    Examples:
        is_union_type(Union[int, str]) == True
    """
    return isinstance(_type, _GenericAlias) and _type.__origin__ is typing.Union


def is_int_enum(_type):
    """Test if the type is a subclass of IntEnum."""
    try:
        return issubclass(_type, IntEnum)
    except TypeError:  # happens when _type is something like List[int]
        return False


def unpack_optional_type(_type):
    """Unpack Optional[<type>] to <type> if necessary

    Examples:
        >>> unpack_optional(Optional[int])
        int
        >>> unpack_optional(int)
        int
    """
    # Optional[x] is actually Union[x, None]
    if not is_union_type(_type):
        return _type
    # Check if it has 2 args and the second one is a NoneType
    if len(_type.__args__) == 2 and isinstance(None, _type.__args__[1]):
        return _type.__args__[0]


def _get_dtype_and_null_value(elem_type):
    if elem_type is int or is_int_enum(elem_type):
        dtype = np.int32
        null_value = -9999
    elif elem_type is float:
        dtype = np.float64
        null_value = np.nan
    elif elem_type is bool:
        dtype = bool
        null_value = False
    else:
        dtype = object
        null_value = None
    return dtype, null_value


def _to_ndarray(value, elem_type, expected_length):
    """Cast value to numpy array, with some type mappings.

    Args:
        value (iterable or None): iterable to convert to ndarray
        elem_type (type)
        expected_length (int or None): the expected number of rows
    """
    # Tuple[...] becomes 2D Array
    if is_tuple_type(elem_type):
        sub_types = elem_type.__args__
        expected_shape = (expected_length, len(sub_types))
        elem_type = sub_types[0]
    else:
        expected_shape = (expected_length,)

    # cast python to numpy dtypes
    dtype, null_value = _get_dtype_and_null_value(elem_type)

    # return empty array if no value or None is supplied
    if value is None:
        return np.full(expected_shape, null_value, dtype=dtype, order="F")
    elif np.isscalar(value) or isinstance(value, IntEnum):
        arr = np.full(expected_shape, value, dtype=dtype, order="F")
    else:
        arr = np.asarray(value, dtype=dtype, order="F")

    if is_int_enum(elem_type):
        allowed = list(elem_type) + [-9999]
        if not np.isin(arr, allowed).all():
            raise ValueError(
                f"{set(value) - set(allowed)} is not a valid {elem_type.__name__}"
            )
    if arr.shape != expected_shape:
        raise ValueError(
            f"Expected an array with shape {expected_shape}, got {arr.shape}."
        )
    return arr


T = TypeVar("T")


class Array(Generic[T]):
    """A dataclass with fields ("columns") that are numpy arrays of equal size."""

    scalars = tuple()

    def __init_subclass__(cls) -> None:
        # Trick to get the record dataclass from the generic typing
        base = cls.__orig_bases__[0]  # type: ignore
        (data_class,) = base.__args__
        cls.validate_data_class(data_class)
        super().__init_subclass__()
        cls.data_class = data_class

    @staticmethod
    def validate_data_class(data_class):
        fields = typing.get_type_hints(data_class)

        # id must be present
        if "id" not in fields:
            raise TypeError("The id field is required")
        if fields["id"] is not int:
            raise TypeError("The id field should be of integer type")

        # check the subtypes of Tuples
        for name, elem_type in fields.items():
            if is_tuple_type(elem_type):
                sub_types = elem_type.__args__
                n_columns = len(sub_types)
                if n_columns == 0:
                    raise TypeError("Tuple cannot contain 0 elements")
                elem_type = sub_types[0]
                if not all(x is elem_type for x in sub_types):
                    raise TypeError("Tuple cannot contain mixed dtypes")

    def __init__(self, **kwargs):
        """Initialize the array dataclass from values.

        None values are converted to arrays with NULL values that are appropriate
        for the type (nan for float, -9999 for int, None for object,
        False for bool).
        """
        # id is required; get its length to use it for validation
        id = kwargs.pop("id", None)
        if id is None:
            raise TypeError("missing required keyword argument 'id'")
        # handle generators
        if hasattr(id, "__iter__") and not hasattr(id, "__len__"):
            id = np.fromiter(id, dtype=np.int32)
        # cast to integer array
        id = np.asarray(id, dtype=np.int32, order="F")
        if id.ndim != 1:
            raise ValueError(f"Expected one-dimensional 'id' array, got {id.ndim}")
        if id.size > 1 and not np.all(id[1:] > id[:-1]):
            raise ValueError("Values in 'id' must be unique and sorted")
        self.id = id

        # fetch the scalar arguments
        for name in self.scalars:
            value = kwargs.pop(name, None)
            if value is None:
                raise TypeError(f"missing required keyword argument '{name}'")
            setattr(self, name, value)

        # convert each field to an array and set it to self
        fields = typing.get_type_hints(self.data_class)
        for name, elem_type in fields.items():
            if name == "id":
                continue
            value = kwargs.pop(name, None)
            if value is not None:
                try:
                    arr = _to_ndarray(value, elem_type, id.size)
                except ValueError as e:
                    raise ValueError(f"Error parsing {name} to array: {e}")
                setattr(self, name, arr)

        # call the init of the wrapped class
        if len(kwargs) > 0:
            raise TypeError(
                f"{self.__class__.__name__} got unexpected arguments {list(kwargs)}"
            )
        super().__init__(**kwargs)

    def __repr__(self):
        return (
            f"<{self.__class__.__name__} object, "
            f"{self.data_class.__name__} array of length {len(self)}>"
        )

    def __len__(self):
        return len(self.id)

    def __add__(self, other: T) -> T:
        """Concatenate two array dataclasses of equal type."""
        if not issubclass(other.__class__, self.__class__):
            raise TypeError(
                f"Cannot concatenate {self} with {other} as they are not of "
                f"equal types."
            )
        new_fields = {
            name: np.concatenate((getattr(self, name), getattr(other, name)))
            for name in typing.get_type_hints(self.data_class).keys()
        }

        for scalar_arg in self.scalars:
            new_fields[scalar_arg] = getattr(self, scalar_arg)
        return self.__class__(**new_fields)

    def id_to_index(self, id, check_exists: bool = False):
        """Find the index of records with given id.

        Args:
            id (int or array_like): The id(s) to find in self.id
            check_exists (bool): Whether to check if the id is actually present. Raises
              a KeyError with added "values" and "indices" attributes if an id is not
              present.

        Returns:
            int or array_like: the indexes into self.id
        """
        # some timings:
        # len(id)    len(self.id)   timing (microseconds)
        # 1000       2000           44
        # 1000       10000          62
        return search(self.id, id, assume_ordered=True, check_exists=check_exists)

    def index_to_id(self, index):
        """Find the id of records with given index.

        Note:
            Index should be 0 <= index < len(self). Values of -9999 will persist.

        Args:
            index (int or array_like): The indices to return from self.id

        Returns:
            int or array_like: the ids from self.id
        """
        index = np.asarray(index)
        return self.take("id", index)

    def __getattr__(self, name):
        try:
            elem_type = typing.get_type_hints(self.data_class)[name]
        except KeyError:
            raise AttributeError
        arr = _to_ndarray(None, elem_type, len(self))
        setattr(self, name, arr)
        return arr

    def take(self, name, index, fill_value=None):
        """Take values from a column, returing 'fill_value' where index == -9999.

        Args:
            name (str): the column name to take values from
            index (int or array_like): The indices to return from self.id
            fill_value: what to return for -9999 indices (defaults to -9999/NaN/None)

        Returns:
            scalar or array_like: the values from column 'name'
        """
        values = getattr(self, name)
        if fill_value is None:
            _, fill_value = _get_dtype_and_null_value(
                typing.get_type_hints(self.data_class)[name]
            )
        index = np.asarray(index)
        return np.where(index != -9999, np.take(values, index, mode="clip"), fill_value)

    def to_dict(self):
        fields = typing.get_type_hints(self.data_class)
        return {field: getattr(self, field) for field in fields}

    def __getitem__(self, idx) -> T:
        """Create a masked copy of this array dataclass"""
        args = {}
        for field in typing.get_type_hints(self.data_class):
            args[field] = getattr(self, field)[idx]
        for scalar_arg in self.scalars:
            args[scalar_arg] = getattr(self, scalar_arg)
        return self.__class__(**args)

    def reorder(self, idx) -> None:
        """Reorder self by given index, inplace.

        Note that this skips self.id: the records are renumbered.
        """
        for field in typing.get_type_hints(self.data_class):
            if field == "id":
                continue
            setattr(self, field, getattr(self, field)[idx])

    def reorder_by(self, attr: str, **kwargs) -> None:
        """Reorder self by given column, inplace.

        Note that this skips self.id: the records are renumbered.
        """
        idx = np.argsort(getattr(self, attr), **kwargs)
        self.reorder(idx)

    def split_in_two(self, by):
        """Split self in to 2 arrays having the first and second occurence of 'by'.

        Returned are index arrays that index into self.

        If 'by' contains only unique values, the second array will be empty. If any value
        occurs more than 2 times, a 'TooManyExist' error is raised, having a 'values'
        attribute of the values that occur more than 2 times.
        """
        if len(by) != len(self):
            raise ValueError("'by' must have the same length as 'self'")
        values, idx1, counts = np.unique(by, return_index=True, return_counts=True)
        if np.any(counts > 2):
            raise TooManyExist("", values=values[counts > 2])
        assert not np.any(counts > 2)
        idx2 = np.delete(np.arange(len(self)), idx1)
        return idx1, idx2


class TooManyExist(KeyError):
    def __init__(self, msg, values):
        self.values = values
        super().__init__(msg)


def replace(arr, mapping, check_present=False):
    """Return array with its values replaced according to ``mapping``.

    If ``check_present`` is False, it is assumed that all elements in the array are keys
    of ``mapping``.
    """
    keys, values = np.array(sorted(mapping.items())).T
    indices = np.digitize(arr, keys, right=True)
    if check_present and not np.all(keys[indices] == arr):
        raise ValueError("Not all values are present in the replacement dict")
    return values[indices]


class DoesNotExist(KeyError):
    def __init__(self, msg, values, indices):
        self.values = values
        self.indices = indices
        super().__init__(msg)


def search(a, v, mask=None, assume_ordered=False, check_exists=True):
    """Find indices where `v` occurs in `a`.

    Args:
      a (1D array_like): Input array. If assume_ordered is True, then it must be sorted
        in ascending order.
      v (array_like): Values to find in a.
      mask (1D array_like): A (boolean or index) mask to apply to a before searching.
      assume_ordered (bool, optional): Whether to assume a is ordered, default False.
      check_exists (bool, optional): Check whether the values in v exist in a,
        default True.

    Returns:
      1D array with the same shape as `v`

    Raises:
      If check_exists is True or if the (masked) ``a`` is empty, this function raises a
      KeyError with added ``values`` and ``indices`` attributes corresponding to the
      missing values and the indices of them into ``v``.
    """
    v = np.asarray(v)
    if v.size == 0:
        return np.empty_like(v, dtype=int)

    if mask is not None:
        mask = np.asarray(mask)
        if mask.dtype == bool:
            mask = np.where(mask)[0]
        a = np.take(a, mask)

    # If there is no array to search in: raise directly
    if len(a) == 0:
        if check_exists:
            raise DoesNotExist(
                "search encountered missing elements",
                values=v,
                indices=np.arange(v.shape[0]),
            )
        return np.full(len(v), -9999, dtype=int)

    if assume_ordered:
        ind = np.searchsorted(a, v)
    else:
        sorter = np.argsort(a)
        ind = np.take(sorter, np.searchsorted(a, v, sorter=sorter), mode="clip")

    missing = np.take(a, ind, mode="clip") != v
    if check_exists and missing.any():
        raise DoesNotExist(
            "search encountered missing elements",
            values=np.compress(missing, v),
            indices=np.where(missing)[0],
        )
    elif missing.any():
        if np.isscalar(ind):
            ind = -9999
        else:
            ind[missing] = -9999

    # Map to original (unmasked) indices
    if mask is not None:
        ind = np.where(ind != -9999, np.take(mask, ind, mode="clip"), -9999)

    return ind
