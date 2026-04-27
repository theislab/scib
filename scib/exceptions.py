class OptionalDependencyNotInstalled(ModuleNotFoundError):
    def __init__(self, exception, module_name=None):
        if module_name is None:
            module_name = exception.name
        original = str(exception) if exception is not None else ""
        self.message = (
            f"\n'{module_name}' is an optional dependency and not installed by default. "
            f"Please make sure you install it manually."
        )
        if original:
            # include the original error message to aid debugging
            self.message += f"\nRaised from: {original}"
        super().__init__(self.message)


class RLibraryNotFound(ModuleNotFoundError):
    def __init__(self, exception, lib_name=None):
        original = str(exception) if exception is not None else ""
        name = lib_name or getattr(exception, "name", None) or "R library"
        self.message = (
            f"\nProblem loading {name}. This usually means the R package is not installed "
            f"or the R environment is not available to Python (rpy2)."
        )
        self.message += (
            "\nTry installing the package in R (e.g. install.packages('pkg')) "
            "or confirm your `rpy2`/R setup."
        )
        if original:
            self.message += f"\nRaised from: {original}"
        super().__init__(self.message)
