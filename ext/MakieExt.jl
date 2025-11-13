module MakieExt

using SLHQuantumSystems
using Makie

__init__() = println("Plotting extension loaded (backend: $(Makie.current_backend()))")

function SLHQuantumSystems.bode(sys::StateSpace,input,output,freq)
    tfs = fresponse_allIO(sys,freq)
    j = first(findall(name->name==output[1],sys.outputs))+output[2]-1
    i = first(findall(name->name==input[1],sys.inputs))+input[2]-1

    tf = tfs[i,j]

    fig = Figure(size = (800,500))
    title = Label(fig[1,1], "Frequency Response of $(sys.name)")
    magax = Axis(fig[2,1],width = 700,height = 300,xscale=log10,yscale=log10,xticksvisible=false,xticklabelsvisible=false)
    phaseax = Axis(fig[3,1],width = 700, height = 100,xscale=log10)
    linkxaxes!(magax,phaseax)

    scatter!(magax,freq,abs.(tf),label = "$input -> $output")
    scatter!(phaseax,freq,angle.(tf))

    axislegend(magax)

    return fig
end

end
