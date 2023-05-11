import inspect
import itertools
import logging
import sys

_log = logging.getLogger(__name__)


def chunkify(iterable, size: int):
    it = iter(iterable)
    while True:
        chunk = tuple(itertools.islice(it, size))
        if not chunk:
            break
        yield chunk


def get_size(obj, seen=None):
    """Recursively finds size of objects in bytes"""
    size = sys.getsizeof(obj)
    if seen is None:
        seen = set()
    obj_id = id(obj)
    if obj_id in seen:
        return 0
    # Important mark as seen *before* entering recursion to gracefully handle
    # self-referential objects
    seen.add(obj_id)
    if hasattr(obj, "__dict__"):
        for cls in obj.__class__.__mro__:
            if "__dict__" in cls.__dict__:
                d = cls.__dict__["__dict__"]
                if inspect.isgetsetdescriptor(d) or inspect.ismemberdescriptor(d):
                    size += get_size(obj.__dict__, seen)
                break
    if isinstance(obj, dict):
        size += sum((get_size(v, seen) for v in obj.values()))
        size += sum((get_size(k, seen) for k in obj))
    elif hasattr(obj, "__iter__") and not isinstance(obj, (str, bytes, bytearray)):
        try:
            size += sum((get_size(i, seen) for i in obj))
        except TypeError:
            _log.exception("Unable to get size of %r. This may lead to incorrect sizes. Please report this error.", obj)
    if hasattr(obj, "__slots__"):  # can have __slots__ with __dict__
        size += sum(get_size(getattr(obj, s), seen) for s in obj.__slots__ if hasattr(obj, s))  # type: ignore

    return size
