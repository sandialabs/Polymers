using BinaryBuilder, Pkg

name = "polymers"
version = v"0.3.7"

# Collection of sources required to build tar
sources = [
    GitSource(
        "https://github.com/sandialabs/Polymers.git",
        "3cc31d8d27bca6b7cdc1b20aabec9f5110467808",
    ),
]

# Bash recipe for building across all platforms
script = raw"""
cd $WORKSPACE/srcdir/polymers

if [[ "${target}" == *-mingw* ]]; then
    export RUSTFLAGS="-Clink-args=-L${libdir}"
fi

cargo build --release --features extern
install -Dvm 755 "target/${rust_target}/release/polymers${exeext}" "${bindir}/polymers${exeext}"
"""

# Some platforms disabled for now due issues with rust and musl cross compilation. See #1673.
platforms = supported_platforms(; experimental = true)
# We dont have all dependencies for armv6l
filter!(p -> arch(p) != "armv6l", platforms)
# Rust toolchain for i686 Windows is unusable
filter!(p -> !Sys.iswindows(p) || arch(p) != "i686", platforms)
platforms = expand_cxxstring_abis(platforms)

# The products that we will ensure are always built
products = [ExecutableProduct("polymers", :polymers)]

# Dependencies that must be installed before this package can be built
dependencies = [
    Dependency("DocStringExtensions"),
]

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
