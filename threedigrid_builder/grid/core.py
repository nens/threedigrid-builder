from enum import Enum
from typing import Tuple
import numpy as np


def _to_ndarray(value, elem_type, name, length):
    """Cast value to numpy array, with some type mappings.
    
    Args:
        value (iterable or None): to convert
        elem_type (type): one of int, float, str, Enum, typing.Tuple
        name (str): the name of the field, for error formatting
        length (int): the expected number of elements
    """
    if value is None:
        return
    # Tuple[...] becomes 2D Array
    if issubclass(elem_type, Tuple):
        sub_types = elem_type.__args__
        n_columns = len(sub_types)
        elem_type = sub_types[0]
    else:
        n_columns = None

    # cast python to numpy dtypes
    if elem_type is int or issubclass(elem_type, Enum):
        elem_type = np.int32
    elif elem_type is float:
        elem_type = np.float64
    elif elem_type is str:
        elem_type = object
    else:
        raise RuntimeError(f"Unknown type for array_of: {elem_type}")

    arr = np.asarray(value, dtype=elem_type)
    if n_columns is None and arr.ndim != 1:
        raise ValueError(
            f"Expected a one dimensional array for {name}, got {arr.ndim}."
        )
    elif n_columns is not None and arr.ndim != 2:
        raise ValueError(
            f"Expected a two dimensional array for {name}, got {arr.ndim}."
        )
    elif arr.ndim == 2 and arr.shape[1] != n_columns:
        raise ValueError(f"Expected {n_columns} columns in {name}, got {arr.shape[1]}.")

    if arr.shape[0] != length:
        raise ValueError(f"Expected {length} rows in {name}, got {arr.shape[0]}.")

    return arr


class array_of:
    """A decorator to create an array of given datastructure.
    
    The decorated class will have the same attributes, but then transformed
    into 1D numpy arrays.

    Supported types are:
    - float (cast to numpy.float64)
    - int (cast to numpy.int32)
    - Enum (cast to numpy.int32)
    - Tuple[...] (cast to numpy.float64)

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

    def __init__(self, _type):
        # perform some checks on _type
        annotations = _type.__annotations__
        for name, elem_type in annotations.items():
            # Tuple[...] becomes 2D Array
            if issubclass(elem_type, Tuple):
                sub_types = elem_type.__args__
                n_columns = len(sub_types)
                if n_columns == 0:
                    raise TypeError("Tuple cannot contain 0 elements")
                elem_type = sub_types[0]
                if not all(x is elem_type for x in sub_types):
                    raise TypeError("Tuple cannot contain mixed dtypes in array_of")

        # check if elem_type is supported
        if elem_type not in (int, float, str) and not issubclass(elem_type, Enum):
            raise ValueError(f"Incompatible field type for array_of: {elem_type}")

        self._type = _type

    def __call__(self, cls):
        class Wrapper(cls):
            def __init__(self, **kwargs):
                length = len(kwargs["id"])

                for name, elem_type in self._type.__annotations__.items():
                    value = kwargs.pop(name, None)
                    setattr(self, name, _to_ndarray(value, elem_type, name, length))
                super().__init__(**kwargs)

            def __repr__(self):
                return f"<array of {self._type.__name__} (len:{len(self.id)})>"

        Wrapper._type = self._type
        return Wrapper
