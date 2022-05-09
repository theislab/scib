class IntegrationMethodNotFound(ModuleNotFoundError):
    def __init__(self, exception):
        self.message = (
            f"\n{exception.name} is not installed by default, "
            "please make sure you install it manually."
        )
        super().__init__(self.message)


class RLibraryNotFound(ModuleNotFoundError):
    def __init__(self, exception):
        self.message = f"\nproblem loading library: {exception}"
        super().__init__(self.message)
