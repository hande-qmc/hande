"""Simple, useful callables when selecting."""


def do_nothing(*args):
    """Do nothing."""
    pass


class RaiseValueError:
    """Raise ValueError."""

    def __init__(self, message):
        self._message = message

    def __call__(self, *args):
        raise ValueError(self._message)
