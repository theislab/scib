import warnings
from functools import wraps

warnings.simplefilter('default')  # or 'always'


def wrap_func_naming(func, name):
    """
    Decorator that adds a `DeprecationWarning` and a name to `func`.
    """

    @wraps(func)
    def wrapper(*args, **kwargs):
        warnings.warn(
            "Mixed case function naming is deprecated."
            "Please use the snake_case version of this file.",
            DeprecationWarning
        )
        return func(*args, **kwargs)

    wrapper.__name__ = name
    return wrapper
