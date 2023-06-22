using BinaryBuilder, Pkg

name = "polymers"
version = v"0.3.7"
sources = [
    GitSource(
        "https://github.com/sandialabs/polymers.git",
        "3cc31d8d27bca6b7cdc1b20aabec9f5110467808",
    ),
]
script = raw"""
cd $WORKSPACE/srcdir/polymers
cargo build --release --features extern
install -Dvm 755 "target/${rust_target}/release/libpolymers.so" "${libdir}/libpolymers.so"
"""
platforms = [supported_platforms()[2]]
products = [LibraryProduct("libpolymers", :libpolymers)]
dependencies = [Dependency("DocStringExtensions")]
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
