using BinaryBuilder, Pkg

name = "polymers"
version = v"0.3.7"

# Collection of sources required to build tar
sources = [
    GitSource(
        "https://github.com/sandialabs/polymers.git",
        "3cc31d8d27bca6b7cdc1b20aabec9f5110467808",
    ),
]

# Bash recipe for building across all platforms
script = raw"""
cd $WORKSPACE/srcdir/polymers
cargo build --release --features extern
"""

# Some platforms disabled for now due issues with rust and musl cross compilation. See #1673.
platforms = [supported_platforms()[2]]

# The products that we will ensure are always built
products = [ExecutableProduct("polymers", :polymers)]

# Dependencies that must be installed before this package can be built
dependencies = [Dependency("DocStringExtensions")]

# Build the tarballs, and possibly a `build.jl` as well.
build_tarballs(
    ARGS,
    name,
    version,
    sources,
    script,
    platforms,
    products,
    dependencies;
    compilers = [:c, :rust],
    preferred_gcc_version = v"7",
    lock_microarchitecture = false,
    julia_compat = "1.9",
)
