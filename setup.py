"""Native-extension build configuration for the PEP 517 setuptools backend."""

from __future__ import annotations

import sys
from pathlib import Path

from setuptools import Extension, setup


def native_extensions() -> list[Extension]:
    import numpy

    if sys.platform == "win32":
        c_compile_args = ["/O2", "/openmp"]
        cpp_compile_args = ["/O2", "/openmp", "/std:c++14"]
        link_args: list[str] = []
    else:
        c_compile_args = ["-O3", "-fopenmp", "-std=c99"]
        cpp_compile_args = ["-O3", "-fopenmp", "-std=c++11"]
        link_args = ["-fopenmp"]
    common = {
        "include_dirs": [numpy.get_include()],
        "define_macros": [("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION")],
    }
    namespace = "hercules.core"
    return [
        Extension(
            f"{namespace}.utils.math_utils",
            ["hercules/core/utils/math_utils.pyx"],
            extra_compile_args=c_compile_args,
            extra_link_args=link_args,
            **common,
        ),
        Extension(
            f"{namespace}.model.vi.e_step",
            ["hercules/core/model/vi/e_step.pyx"],
            extra_compile_args=c_compile_args,
            extra_link_args=link_args,
            **common,
        ),
        Extension(
            f"{namespace}.model.vi.e_step_cpp",
            ["hercules/core/model/vi/e_step_cpp.pyx"],
            language="c++",
            extra_compile_args=cpp_compile_args,
            extra_link_args=link_args,
            **common,
        ),
    ]


extensions = native_extensions()
if extensions:
    from Cython.Build import cythonize

    extensions = cythonize(
        extensions,
        build_dir=str(Path("build") / "cython"),
        compiler_directives={
            "language_level": 3,
            "embedsignature": True,
            "boundscheck": False,
            "wraparound": False,
            "nonecheck": False,
            "cdivision": True,
        },
    )

setup(ext_modules=extensions)
