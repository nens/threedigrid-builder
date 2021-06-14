from enum import IntEnum

import numpy as np
import sys
import typing


USE_PEP560 = sys.version_info[:3] >= (3, 7, 0)
if USE_PEP560:
    from typing import _GenericAlias
else:
    from typing import _Union
    from typing import TupleMeta

__all__ = ["array_of", "is_tuple_type", "is_int_enum", "unpack_optional_type"]


def is_tuple_type(_type):
    """Test if the type is a indexed typing.Tuple type.

    Examples:
        is_tuple_type(Tuple[int, int]) == True
    """
    if USE_PEP560:
        return isinstance(_type, _GenericAlias) and _type.__origin__ is tuple
    else:
        return type(_type) is TupleMeta


def is_union_type(_type):
    """Test if the type is an typing.Union type.

    Examples:
        is_union_type(Union[int, str]) == True
    """
    if USE_PEP560:
        return isinstance(_type, _GenericAlias) and _type.__origin__ is typing.Union
    else:
        return type(_type) is _Union


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


class ArrayDataClass:
    """A dataclass with fields ("columns") that are numpy arrays of equal size."""

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

        # convert each field to an array and set it to self
        fields = self.data_class.__annotations__
        for name, elem_type in fields.items():
            if name == "id":
                continue
            value = kwargs.pop(name, None)
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

    def __add__(self, other):
        """Concatenate two array dataclasses of equal type."""
        if self.__class__ is not other.__class__:
            raise TypeError(
                f"Cannot concatenate {self} with {other} as they are not of "
                f"equal types."
            )
        new_fields = {
            name: np.concatenate((getattr(self, name), getattr(other, name)))
            for name in self.data_class.__annotations__.keys()
        }
        return self.__class__(**new_fields)

    def id_to_index(self, id):
        """Find the index of records with given id.

        Note:
            It is not checked if the id actually exists! In that case
            behaviour of numpy.searchsorted is produced: the returned index is
            the index where to insert id to maintain order in self.id.

        Args:
            id (int or array_like): The id(s) to find in self.id

        Returns:
            int or array_like: the indexes into self.id
        """
        # some timings:
        # len(id)    len(self.id)   timing (microseconds)
        # 1000       2000           44
        # 1000       10000          62
        return np.searchsorted(self.id, id)  # inverse of self.id[index]

    def index_to_id(self, index):
        """Find the id of records with given index.

        Note:
            Index should be 0 <= index < len(self).

        Args:
            index (int or array_like): The indices to return from self.id

        Returns:
            int or array_like: the ids from self.id
        """
        # some timings:
        # len(id)    len(self.id)   timing (microseconds)
        # 1000       2000           5.2
        # 1000       10000          5.2
        return np.take(self.id, index)  # same as self.id[index]

    def to_dict(self):
        fields = self.data_class.__annotations__
        return {field: getattr(self, field) for field in fields}

    def reorder_by(self, idx):
        """Reorder self by given idx, inplace.

        Note that this skips self.id, effectively the records are renumbered.
        """
        for field in self.data_class.__annotations__:
            if field == "id":
                continue
            setattr(self, field, getattr(self, field)[idx])


class array_of:
    """A decorator to create an array dataclass from a normal dataclass.

    The decorated class will have the same attributes as the provided dataclass
    but then transformed into 1D numpy arrays.

    Supported types are:
    - float (becomes 1D numpy array of numpy.float64)
    - int (becomes 1D numpy array of numpy.int32)
    - bool (becomes 1D numpy array of numpy._bool)
    - Enum (becomes 1D numpy array of numpy.int32)
    - Tuple[<type>] (becomes 2D numpy array of <type>)

    Other types are not converted to numpy dtypes and kept as python object.

    The id field is required.

    Example:
    >>> class MyRecord:
    ...     id: int
    ...
    >>> @array_of(MyRecord)
    >>> class MyRecords:
    >>>     pass
    ...
    >>> MyRecords(id=[1, 2, 3])
    < array of MyRecord (len:3) >
    """

    def __init__(self, data_class):
        """Validate the dataclass passed into the decorator"""
        fields = data_class.__annotations__

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

        self.data_class = data_class

    def __call__(self, cls):
        """Transform a class into a ArrayDataClass subclass.

        Note that the normal decorator usage is equivalent to::

        >>> class MyRecords:
        ...     pass
        ...
        >>> MyRecords = array_of(MyRecord)(MyRecords)
        """

        class Wrapper(ArrayDataClass, cls):
            pass

        Wrapper.data_class = self.data_class

        # transfers __name__, __doc__, etc. from cls to Wrapper
        Wrapper.__doc__ = cls.__doc__
        Wrapper.__name__ = cls.__name__
        Wrapper.__qualname__ = cls.__qualname__
        Wrapper.__module__ = cls.__module__
        return Wrapper
