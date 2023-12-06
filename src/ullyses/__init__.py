from ._release import set_release
from importlib.metadata import version

__version__ = version(__name__)

__release__ = set_release()