cd("..")
if Sys.isunix()
    run(`cargo build --features extern`)
elseif Sys.iswindows()
    run(`cmd /c cargo build --features extern`)
end
