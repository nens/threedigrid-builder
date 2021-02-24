from enum import IntEnum
from typing import Tuple
import numpy as np


def get_type_class(typ):
    try:
        # Python 3.5 / 3.6
        return typ.__extra__
    except AttributeError:
        # Python 3.7
        return typ.__origin__


def _to_ndarray(value, elem_type, expected_length):
    """Cast value to numpy array, with some type mappings.

    Args:
        value (iterable or None): iterable to convert to ndarray
        elem_type (type)
        expected_length (int or None): the expected number of rows
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
    if elem_type is int or get_type_class(elem_type) is tuple:
        dtype = np.int32
    elif elem_type is float:
        dtype = np.float64
    elif elem_type is bool:
        dtype = bool
    else:
        dtype = object

    arr = np.asarray(value, dtype=dtype, order="F")
    if n_columns is None and arr.ndim != 1:
        raise ValueError(f"Expected a one dimensional array, got {arr.ndim}.")
    elif n_columns is not None and arr.ndim != 2:
        raise ValueError(f"Expected a two dimensional array, got {arr.ndim}.")
    elif arr.ndim == 2 and arr.shape[1] != n_columns:
        raise ValueError(f"Expected {n_columns} columns, got {arr.shape[1]}.")

    if expected_length is not None and arr.shape[0] != expected_length:
        raise ValueError(f"Expected {expected_length} rows, got {arr.shape[0]}.")

    return arr


class array_of:
    """A decorator to create an array dataclass from a normal dataclass.

    The decorated class will have the same attributes as the provided dataclass
    but then transformed into 1D numpy arrays.

    Supported types are:
    - float (cast to numpy.float64)
    - int (cast to numpy.int32)
    - bool (cast to numpy._bool)
    - Enum (cast to numpy.int32)
    - Tuple[...] (cast to numpy.float64)

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

        # check the subtypes of Tuples
        fields = data_class.__annotations__
        for name, elem_type in fields.items():
            if get_type_class(elem_type) is tuple:
                sub_types = elem_type.__args__
                n_columns = len(sub_types)
                if n_columns == 0:
                    raise TypeError("Tuple cannot contain 0 elements")
                elem_type = sub_types[0]
                if not all(x is elem_type for x in sub_types):
                    raise TypeError("Tuple cannot contain mixed dtypes")

        self.data_class = data_class

    def __call__(self, cls):
        """Wrap a class by subclassing it.

        Note that the below example is equal to the normal decorator usage
        illustrated above.

        >>> class MyRecords:
        ...     pass
        ...
        >>> MyRecords = array_of(MyRecord)(MyRecords)
        """

        class Wrapper(cls):
            def __init__(self, **kwargs):
                """Initialize the array dataclass

                Take the iterables in the arguments and convert them into 1D
                ndarrays according to the types in the normal dataclass.
                """
                # id is required; get its length to use it for validation
                try:
                    expected_length = len(kwargs["id"])
                except KeyError:
                    raise TypeError("missing required keyword argument 'id")

                # convert each field to an array and set it to self
                fields = self.data_class.__annotations__
                for name, elem_type in fields.items():
                    value = kwargs.pop(name, None)
                    try:
                        arr = _to_ndarray(value, elem_type, expected_length)
                    except ValueError as e:
                        raise ValueError(f"Error parsing {name} to array: {e}")
                    setattr(self, name, arr)

                # call the init of the wrapped class
                super().__init__(**kwargs)

            def __repr__(self):
                return f"<array of {self._type.__name__} (len:{len(self.id)})>"

        Wrapper.data_class = self.data_class

        # transfers __name__, __doc__, etc. from cls to Wrapper
        Wrapper.__doc__ = cls.__doc__
        Wrapper.__name__ = cls.__name__
        Wrapper.__qualname__ = cls.__qualname__
        Wrapper.__module__ = cls.__module__
        return Wrapper
