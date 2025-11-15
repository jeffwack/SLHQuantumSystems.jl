using SLHQuantumSystems

cavA = cavity("A")
cavB = cavity("B")

sys = concatenate([cavA,cavB],"sys")

sys = feedbackreduce(sys,"B_out","A_in")
#=
println(sys.parameters)
println(sys.operators)
println(sys.H)
=#
