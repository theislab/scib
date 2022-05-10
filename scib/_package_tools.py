import inspect
import warnings
from functools import wraps

warnings.simplefilter("default")  # or 'always'


def wrap_func_naming(func, name):
    """
    Decorator that adds a `DeprecationWarning` and a name to `func`.
    """

    @wraps(func)
    def wrapper(*args, **kwargs):
        warnings.warn(
            f"Mixed case function naming is deprecated for '{name}'. "
            f"Please use '{func.__name__}' instead.",
            DeprecationWarning,
        )
        return func(*args, **kwargs)

    wrapper.__name__ = name
    return wrapper


def rename_func(function, new_name):
    if callable(function):
        function = wrap_func_naming(function, new_name)
    setattr(inspect.getmodule(function), new_name, function)
