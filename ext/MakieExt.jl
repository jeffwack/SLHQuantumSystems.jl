module MakieExt

using SLHQuantumSystems
using Makie

__init__() = println("Plotting extension loaded (backend: $(Makie.current_backend()))")

function SLHQuantumSystems.bode()
    return scatter([1,2,3],[3,2,1])
end

end
