import pkg_resources  # part of setuptools

program_version = pkg_resources.require("anonymise")[0].version
