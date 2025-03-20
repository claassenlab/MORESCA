from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("MORESCA")
except PackageNotFoundError:
    pass
