class OptionalDependencyNotInstalled(ModuleNotFoundError):
    def __init__(self, exception):
        self.message = (
            f"\n{exception.name} is an optional dependency and not installed by default. "
            f"Please make sure you install it manually."
        )
        super().__init__(self.message)


class RLibraryNotFound(ModuleNotFoundError):
    def __init__(self, exception):
        self.message = f"\nproblem loading library: {exception}"
        super().__init__(self.message)
