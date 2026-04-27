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
            stacklevel=2,
        )
        return func(*args, **kwargs)

    wrapper.__name__ = name
    return wrapper


def rename_func(function, new_name):
    if callable(function):
        function = wrap_func_naming(function, new_name)
    setattr(inspect.getmodule(function), new_name, function)


def renamed_arg(old_name, new_name, *, pos_0: bool = False):
    """
    Taken from: https://github.com/scverse/scanpy/blob/214e05bdc54df61c520dc563ab39b7780e6d3358/scanpy/_utils/__init__.py#L130C1-L157C21
    """

    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            if old_name in kwargs:
                f_name = func.__name__
                pos_str = (
                    (
                        f" at first position. Call it as `{f_name}(val, ...)` "
                        f"instead of `{f_name}({old_name}=val, ...)`"
                    )
                    if pos_0
                    else ""
                )
                msg = (
                    f"In function `{f_name}`, argument `{old_name}` "
                    f"was renamed to `{new_name}`{pos_str}."
                )
                warnings.warn(msg, FutureWarning, stacklevel=3)
                if pos_0:
                    args = (kwargs.pop(old_name), *args)
                else:
                    kwargs[new_name] = kwargs.pop(old_name)
            return func(*args, **kwargs)

        return wrapper

    return decorator
