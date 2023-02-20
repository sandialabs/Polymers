module Polymers

if Sys.isunix()
    const PATHSEP = "/"
elseif Sys.iswindows()
    const PATHSEP = "\\"
end

const PROJECT_ROOT = string(dirname(@__FILE__), PATHSEP, "..", PATHSEP)

include(string("physics", PATHSEP, "mod.jl"))

end
