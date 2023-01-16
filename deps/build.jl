using Pkg

Pkg.activate(".")

run(`cargo build --features julia`)