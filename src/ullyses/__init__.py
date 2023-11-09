import sys

# Import toml library to read tool table from pyproject.toml
if sys.version_info.major == 3 and sys.version_info.minor < 11:
    import tomli
    loader = tomli.load
else:
    import tomllib
    loader = tomllib.load

# Load pyproject.toml
with open('pyproject.toml', 'rb') as f:
    project_file = loader(f)

# Pull release string from ullyses tool table
__release__ = project_file['tool']['ullyses']['release']