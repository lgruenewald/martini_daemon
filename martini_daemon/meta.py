# Python metaprogramming helpers
import functools
from typing import Callable


def alias(aliases: dict[str, str]) -> Callable:
    """
    Wrap a function with this to provide a dictionary of aliases for
    named arguments.
    """
    def decorator(func: Callable):
        @functools.wraps(func)
        def inner(*args, **kwargs):
            for k, v in aliases.items():
                if kwargs.get(k) is not None:
                    kwargs[v] = kwargs.get(k)
                    del kwargs[k]
            return func(*args, **kwargs)
        return inner
    return decorator
