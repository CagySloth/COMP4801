from setuptools import setup, Extension
from distutils.sysconfig import customize_compiler
from Cython.Build import cythonize
from setuptools.command.build_ext import build_ext as _build_ext

import os

HERE = os.path.abspath(os.path.dirname(__file__))


def CppExtension(name, sources):
    # Prepend our local paths
    full_sources = [os.path.join(HERE, s) for s in sources]
    return Extension(
        name,
        sources=full_sources,
        language="c++",
        # extra_compile_args=["-std=c++11", "-Werror=return-type", "-Werror=narrowing"],
        extra_compile_args=["-std=c++11"],
        include_dirs=[os.path.join(HERE, "src")],
        undef_macros=["NDEBUG"],
    )


extensions = [
    # 1) Core DP + ReadSet wrapper
    CppExtension(
        "whatshap.core",
        sources=[
            "whatshap/core.pyx",
            "src/backwardcolumniterator.cpp",
            "src/binomial.cpp",
            "src/caller.cpp",
            "src/columnindexingiterator.cpp",
            "src/columnindexingscheme.cpp",
            "src/columniterator.cpp",
            "src/entry.cpp",
            "src/genotype.cpp",
            "src/genotypecolumncostcomputer.cpp",
            "src/genotypedistribution.cpp",
            "src/genotypedptable.cpp",
            "src/genotyper.cpp",
            "src/graycodes.cpp",
            "src/indexset.cpp",
            "src/multinomial.cpp",
            "src/pedigree.cpp",
            "src/pedigreecolumncostcomputer.cpp",
            "src/pedigreedptable.cpp",
            "src/pedigreepartitions.cpp",
            "src/pedmecheuristic.cpp",
            "src/phredgenotypelikelihoods.cpp",
            "src/read.cpp",
            "src/readset.cpp",
            "src/transitionprobabilitycomputer.cpp",
            "src/hapchat/balancedcombinations.cpp",
            "src/hapchat/basictypes.cpp",
            "src/hapchat/binomialcoefficient.cpp",
            "src/hapchat/hapchatcolumniterator.cpp",
            "src/hapchat/hapchatcore.cpp",
        ],
    ),

    # 2) Read selection Cython extension (no extra C++)
    CppExtension(
        "whatshap.readselect",
        sources=["whatshap/readselect.pyx"],
    ),

    # 3) Priority queue Cython extension
    CppExtension(
        "whatshap.priorityqueue",
        sources=["whatshap/priorityqueue.pyx"],
    ),
]


class BuildExt(_build_ext):
    def build_extensions(self):
        customize_compiler(self.compiler)
        # Remove the C compiler flag that annoys C++ builds
        try:
            self.compiler.compiler_so.remove("-Wstrict-prototypes")
        except (AttributeError, ValueError):
            pass
        super().build_extensions()


setup(
    name="whcore",
    version="0.0.1",
    packages=["whatshap"],
    package_dir={"whatshap": "whatshap"},
    ext_modules=cythonize(extensions),
    cmdclass={"build_ext": BuildExt},
)
