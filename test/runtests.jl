using SLHQuantumSystems
using Test
using Aqua
using SecondQuantizedAlgebra
using Symbolics
using LinearAlgebra
using ControlSystemsBase
using ControlSystemsBase: nstates, ninputs, noutputs

# Test systems on FockSpace subspaces, with ports in and out, or in1, in2, … and out1, out2, …:
# fock has one unnamed subspace, with operator a; fock2 has subspaces a and b, with operators a_a and b_a
ports(n, dir) = n == 1 ? [dir] : ["$dir$j" for j in 1:n]
fock(name, S, L, H) = SLH(name, [OpticalMode("")], ports(length(L), "in"), ports(length(L), "out"), S, L, H)
const TWOMODES = [GenericSubspace("a", FockSpace), GenericSubspace("b", FockSpace)]
fock2(name, S, L, H) = SLH(name, TWOMODES, ports(length(L), "in"), ports(length(L), "out"), S, L, H)

@testset "SLHQuantumSystems.jl" begin

    @testset "Aqua" begin
        # Ambiguities are skipped, since most reported ones are between methods of our dependencies.
        # Base.conj(::Complex{Num}) is deliberate piracy (see abcd.jl); remove it from treat_as_own along with that method.
        Aqua.test_all(SLHQuantumSystems; ambiguities = false, piracies = (treat_as_own = [conj],))
    end

    @testset "SLH constructor" begin
        # Create test system
        hf = FockSpace(:cavity)
        @qnumbers a::Destroy(hf)
        @variables ω κ g

        H = ω * a' * a
        L = [√κ * a]
        S = [1]
        modes = [GenericSubspace("", FockSpace)]
        sys = SLH("test", modes, ["in"], ["out"], S, L, H)

        # SLH construction
        @test sys.name == "test"
        @test sys.inputs == ["in"]
        @test sys.outputs == ["out"]
        @test sys.S == [1]
        @test length(sys.L) == 1

        # Ports and S must match the channels of L, and port names must be unique
        @test_throws "L has 1 channels, but there are 2 inputs" SLH("t", modes, ["i1", "i2"], ["o"], [1], [√κ * a], H)
        @test_throws "S must be 1×1" SLH("t", modes, ["i"], ["o"], [1 0; 0 1], [√κ * a], H)
        @test_throws "S must be 2×2" SLH("t", modes, ["i1", "i2"], ["o1", "o2"], [1, 1], [√κ * a, √κ * a'], H)
        @test_throws "input names must be unique" SLH("t", modes, ["i", "i"], ["o1", "o2"], [1 0; 0 1], [√κ * a, √κ * a'], H)
        @test_throws "output names must be unique" SLH("t", modes, ["i1", "i2"], ["o", "o"], [1 0; 0 1], [√κ * a, √κ * a'], H)

        # With no operators there are no subspaces: a beamsplitter is pure scattering, D = S in quadratures
        @variables θ
        bs = SLH("bs", Subspace[], ["in1", "in2"], ["out1", "out2"], [cos(θ) im * sin(θ); im * sin(θ) cos(θ)], [0, 0], 0)
        @test isempty(bs.subspaces)
        q = substitute(QuantumStateSpace(bs), Dict(θ => 0.3))
        @test nstates(q) == 0 && eltype(q.A) == Float64
        T = transfer_function(q, [1.0, 2.0])
        @test T["out1_x", "in1_x"] ≈ fill(cos(0.3), 2) && T["out1_x", "in2_p"] ≈ fill(-sin(0.3), 2)

        # Two systems of the same name would share ports and parameters, even when their modes are named
        m = [GenericSubspace("m", FockSpace)]
        ma, = operators(m, "m")
        named = SLH("t", m, ["in"], ["out"], [1], [√κ * ma], ω * ma' * ma)
        @test_throws "must have unique names" concatenate("top", named, named)

        # Names given by the user cannot contain '_', which separates the levels of promoted names
        @test_throws "system name \"my_t\"" SLH("my_t", modes, ["in"], ["out"], S, L, H)
        @test_throws "port name \"in_x\"" SLH("t", modes, ["in_x"], ["out"], S, L, H)
        @test_throws "subspace name \"my_m\"" OpticalMode("my_m")
        @test_throws "system name \"my_top\"" concatenate("my_top", sys)
        @test_throws "system name \"my_top\"" concatenate("my_top", QuantumStateSpace(sys))
        @test subspace_names(concatenate("top", named)) == ["t_m"]

    end

    @testset "SLH composition" begin
        cavityA = cavity("A")
        cavityB = cavity("B")

        # Test that concatenate runs without error
        combined = concatenate("S", cavityA, cavityB)
        @test isa(combined, SLH)

        # Test that feedbackreduce runs without error
        cascaded = feedbackreduce(combined, "A_out" => "B_in")
        @test isa(cascaded, SLH)

        # Feeding a cavity's output back into its own input is a loop with unit gain, which has no solution
        @test_throws ErrorException feedbackreduce(cavityA, "out" => "in")
        @test_throws ErrorException feedbackreduce(combined, "A_out" => "A_in")

        # An unknown port name is an error naming the ports, and a swapped pair says so
        @test_throws "no output port named \"out\"" feedbackreduce(combined, "out" => "B_in")
        @test_throws "no input port named \"B_inn\"" feedbackreduce(combined, "A_out" => "B_inn")
        @test_throws "the order is output => input" feedbackreduce(combined, "B_in" => "A_out")

        # concatenate promotes the parameters of S along with those of L and H, and leaves numbers alone
        hf = FockSpace(:cavity)
        @qnumbers a::Destroy(hf)
        @variables θ κ ω bs_θ bs2_θ
        bs(name) = fock(name, [cos(θ)], [sqrt(κ) * sin(θ) * a], ω * a' * a)
        promoted = concatenate("top", bs("bs"))
        @test isequal(promoted.S[1, 1], cos(bs_θ))
        @test any(v -> isequal(v, bs_θ), Symbolics.get_variables(promoted.L[1]))
        two = concatenate("top", bs("bs"), bs("bs2"))
        @test isequal(two.S, [cos(bs_θ) 0; 0 cos(bs2_θ)])
        @test combined.S == [1 0; 0 1]
    end

    @testset "SLH to ABCD" begin
        qss = QuantumStateSpace(cavity("A"))
        @test isa(qss, QuantumStateSpace)
    end

    @testset "operator validation" begin
        hf = FockSpace(:cavity)
        @qnumbers a::Destroy(hf)
        @variables ω κ λ m

        # PhaseSpace operators are rejected by the annihilation-creation builder
        ps = PhaseSpace(:mass)
        x = Position(ps, :x)
        p = Momentum(ps, :p)
        mass(L, H) = SLH("mass", [MechanicalMode("", PhaseSpace)], ["in"], ["out"], [1], L, H)
        @test_throws "use the quadrature basis" QuantumStateSpace(mass([1im * λ * x], p * p / (2m)), AnnihilationCreationBasis())

        # Every operator must act on a declared subspace, of its kind, and be named after it
        @test_throws "which is not a FockSpace" mass([κ * a], p * p / (2m))
        hp = FockSpace(:one) ⊗ FockSpace(:two)
        b = Destroy(hp, :b, 2)
        @test_throws "only 1 subspaces are declared" fock("bad", [1], [κ * b], ω * b' * b)
        @test_throws "must be named a" fock("bad", [1], [κ * Destroy(hf, :b)], 0)

        # A declared subspace that no operator touches is allowed, and is decoupled
        sys = SLH("spare", [GenericSubspace("", FockSpace), GenericSubspace("spare", FockSpace)], ["in"], ["out"], [1], [κ * a], ω * a' * a)
        qss = QuantumStateSpace(sys, AnnihilationCreationBasis())
        @test size(qss.A) == (4, 4)
        @test all(iszero, qss.A[3:4, :]) && all(iszero, qss.A[:, 3:4])
        @test all(iszero, qss.B[3:4, :])
        @test all(iszero, qss.C[:, 3:4])

        # A constant in L is a coherent drive, which the linear model rejects; a zero constant is fine
        @variables β
        @test_throws ErrorException QuantumStateSpace(fock("drive", [1], [β + κ * a], ω * a' * a), AnnihilationCreationBasis())
        @test_throws ErrorException QuantumStateSpace(
            SLH("drive", [GenericSubspace("", FockSpace)], ["i1", "i2"], ["o1", "o2"], [1 0; 0 1], [κ * a, β], ω * a' * a),
            AnnihilationCreationBasis()
        )
        @test isa(QuantumStateSpace(fock("nodrive", [1], [0 + κ * a], ω * a' * a), AnnihilationCreationBasis()), QuantumStateSpace)

        # A c-number H only shifts the energy, so H = 0 or H = ω builds a pure-loss mode in either basis
        for H in (0, ω), basis in (AnnihilationCreationBasis(), QuadratureBasis())
            @test isequal(QuantumStateSpace(fock("loss", [1], [sqrt(κ) * a], H), basis).A, QuantumStateSpace(fock("loss", [1], [sqrt(κ) * a], 0 * a), basis).A)
        end

        # A bare operator is coerced to an iterable QAdd
        bare = fock("t", [1], [a], ω * a' * a)
        @test all(Lj -> Lj isa SecondQuantizedAlgebra.QAdd, bare.L)
        @test isequal(
            QuantumStateSpace(bare, AnnihilationCreationBasis()).A,
            QuantumStateSpace(fock("t", [1], [1 * a], ω * a' * a), AnnihilationCreationBasis()).A
        )
    end

    # These pin the conventions in the docs (LinearQuantumNetworks.md, "Conventions")
    @testset "decay-rate convention" begin
        @variables κ Δ
        κ0, Δ0 = 0.3, 0.7
        qss = substitute(QuantumStateSpace(cavity("A")), Dict(κ => κ0, Δ => Δ0))

        # L = √κ a puts the poles at -κ/2, so κ is the energy decay rate (FWHM)
        @test sort(eigvals(qss.A); by = imag) ≈ [-κ0 / 2 - im * Δ0, -κ0 / 2 + im * Δ0]
        @test qss.B ≈ -sqrt(κ0) * I(2)
        @test qss.C ≈ sqrt(κ0) * I(2)
        @test qss.D ≈ I(2)

        # Full substitution gives concrete eltypes; partial substitution keeps the symbols
        @test all(isconcretetype ∘ eltype, (qss.A, qss.B, qss.C, qss.D))
        partial = substitute(QuantumStateSpace(cavity("A")), Dict(κ => κ0))
        @test any(x -> !(x isa Number), partial.A)
    end

    @testset "names" begin
        @variables κ Δ
        κ0, Δ0 = 0.3, 0.7

        # An unnamed mode takes the system's name on concatenate, as its ports do
        @test subspace_names(cavity("A")) == [""]
        @test state_names(QuantumStateSpace(cavity("A"))) == [:x, :p]
        @test state_names(QuantumStateSpace(cavity("A"), AnnihilationCreationBasis())) == [:a, :adag]
        both = concatenate("S", cavity("A"), cavity("B"))
        @test subspace_names(both) == ["A", "B"]
        @test state_names(QuantumStateSpace(both)) == [:A_x, :A_p, :B_x, :B_p]

        # input_names and output_names label the columns and rows of D, two per port, as
        # ControlSystemsBase expects; input_ports and output_ports name the ports alone
        qboth = QuantumStateSpace(both)
        @test input_ports(both) == input_ports(qboth) == ["A_in", "B_in"]
        @test output_ports(qboth) == ["A_out", "B_out"]
        @test input_names(qboth) == ["A_in_x", "A_in_p", "B_in_x", "B_in_p"]
        @test output_names(QuantumStateSpace(both, AnnihilationCreationBasis())) == ["A_out_a", "A_out_adag", "B_out_a", "B_out_adag"]
        @test length(input_names(qboth)) == ninputs(qboth)
        @test length(output_names(qboth)) == noutputs(qboth)
        @test input_names(qboth, 3) == "B_in_x"
        @test output_names(qboth, 2) == "A_out_p"
        @test state_names(qboth, 4) == :B_p

        # Mode names must be unique, which catches concatenating two systems of the same name
        @test_throws ErrorException concatenate("S", cavity("A"), cavity("A"))
        hp = FockSpace(:one) ⊗ FockSpace(:two)
        dup = [GenericSubspace("", FockSpace), GenericSubspace("", FockSpace)]
        @test_throws "subspace names must be unique" SLH("dup", dup, ["in"], ["out"], [1], [κ * Destroy(hp, :a, 1)], 0)

        # Strings name ports and Symbols name states; a name without its suffix selects both entries
        @variables A_κ A_Δ B_κ B_Δ
        qss = substitute(QuantumStateSpace(both), Dict(A_κ => κ0, A_Δ => Δ0, B_κ => κ0, B_Δ => Δ0))
        ω = [0.1, 1.0, 3.0]
        T = transfer_function(qss, ω)
        G = freqresp(qss, ω)
        @test T["A_out_p", "A_in_x"] == G[2, 1, :]
        @test T["A_out", "A_in"] == G[1:2, 1:2, :]
        @test T[:A_x, "A_in_x"] ≈ [-sqrt(κ0) * (κ0 / 2 + im * w) / ((im * w + κ0 / 2)^2 + Δ0^2) for w in ω]
        @test T["A_out_x", :A_x] ≈ sqrt(κ0) * T[:A_x, :A_x]
        @test_throws ErrorException T[:A_out_p, :A_in_x]
        @test_throws ErrorException T["A_ou", "A_in_x"]

        # ss drops to a plain StateSpace, for the functions that build a new system
        plain = ss(qss)
        @test plain isa StateSpace{Continuous, Float64}
        @test ssdata(plain) == ssdata(qss)
        @test convert(StateSpace, qss) == plain
        @test freqresp(plain * plain, ω) ≈ stack(G[:, :, k]^2 for k in eachindex(ω))
        @test ss(QuantumStateSpace(both)) isa StateSpace
    end

    @testset "subspaces, operators and parameters" begin
        # The operators and parameters of a subspace are named after it, by the types of its space and of itself
        subs = [OpticalMode("cavity"), MechanicalMode("mirror", PhaseSpace), GenericSubspace("atom", NLevelSpace, 2)]
        @test hilbert(subs) == FockSpace(:cavity) ⊗ PhaseSpace(:mirror) ⊗ NLevelSpace(:atom, 2)
        a, = operators(subs, "cavity")
        x, p = operators(subs, "mirror")
        σ = operators(subs, "atom")
        @test string.((a, x, p)) == ("cavity_a", "mirror_x", "mirror_p") && size(σ) == (2, 2)
        @test string.(parameters(subs, "cavity")) == ["cavity_L"] && string.(parameters(subs, "mirror")) == ["mirror_m"]
        @test isempty(parameters(subs, "atom"))
        @test_throws "no subspace named \"mass\"" operators(subs, "mass")

        # concatenate promotes the subspace and everything named after it together
        m, = parameters(subs, "mirror")
        @variables κ g
        arm = SLH("arm", subs[1:2], ["in"], ["out"], [1], [sqrt(κ) * a], p^2 / (2m) - g * x * (a' + a))
        ifo = concatenate("ifo", arm)
        @test subspace_names(ifo) == ["arm_cavity", "arm_mirror"]
        @test string.(operators(ifo)) == ["arm_cavity_a", "arm_mirror_x", "arm_mirror_p"]
        @test issubset(Set(operators(ifo)), Set(SLHQuantumSystems.qops(ifo.H)))
        @test string.(parameters(ifo)) == ["arm_cavity_L", "arm_g", "arm_mirror_m", "arm_κ"]
        mm, = parameters(ifo, "arm_mirror")
        @test any(v -> isequal(v, mm), Symbolics.get_variables(ifo.H))
        @test subspace(ifo, "arm_mirror") isa MechanicalMode{PhaseSpace}
        @test_throws "no subspace named \"mirror\"" subspace(ifo, "mirror")

        # The state space keeps its subspaces, and so its operators and parameters
        q = substitute(QuantumStateSpace(ifo), Dict(mm => 2.0))
        @test string.(operators(q, "arm_mirror")) == ["arm_mirror_x", "arm_mirror_p"]
        @test string.(parameters(q)) == ["arm_cavity_L", "arm_g", "arm_mirror_m", "arm_κ"]
    end

    @testset "transfer function" begin
        @variables s κ Δ A_κ A_Δ B_κ B_Δ
        vals = Dict(A_κ => 0.3, A_Δ => 0.7, B_κ => 0.5, B_Δ => -0.2)
        ω = [0.4, 2.0]
        # A symbolic entry at s = iω, with the parameters in vals
        at(G, w) = Symbolics.value(substitute(G, merge(Dict{Any, Any}(vals), Dict(s => im * w)); fold = Val(true)))

        # The symbolic transfer function agrees with the numeric one at s = iω, between ports and states
        sys = QuantumStateSpace(feedbackreduce(concatenate("S", cavity("A"), cavity("B")), "A_out" => "B_in"))
        T = transfer_function(sys, s)
        Tω = transfer_function(substitute(sys, vals), ω)
        for i in ["B_out_x", "B_out_p", :A_x, :B_p], j in ["A_in_x", "A_in_p", :A_p]
            @test [at(T[i, j], w) for w in ω] ≈ Tω[i, j]
        end
        @test size(T["B_out", "A_in"]) == (2, 2)

        # A short name selects only its own quadratures, not names that extend it by more words
        subs2 = [OpticalMode(""), MechanicalMode("mirror")]
        (a2,), (b2,) = operators(subs2, ""), operators(subs2, "mirror")
        arm = SLH("arm", subs2, ["in", "inaux"], ["out", "outaux"],
                  [1 0; 0 1], [sqrt(0.3) * a2, sqrt(0.5) * b2], 0.2 * (a2' * b2 + b2' * a2))
        Tarm = transfer_function(QuantumStateSpace(concatenate("net", arm)), ω)
        @test size(Tarm["arm_out", "arm_in"]) == (2, 2, 2)
        @test size(Tarm[:arm, "arm_inaux_x"]) == (2, 1, 2)

        # The cavity's phase response to amplitude, and an entry that is identically zero
        cav = transfer_function(QuantumStateSpace(cavity("A")), s)
        @test isequal(simplify(cav["out_p", "in_x"] - Δ * κ / (s^2 + κ * s + Δ^2 + κ^2 / 4)), 0)
        unconnected = transfer_function(QuantumStateSpace(concatenate("S", cavity("A"), cavity("B"))), s)
        @test Symbolics.value(unconnected["B_out_x", "A_in_x"]) == 0

        # Numeric and partly numeric systems give a polynomial ratio in s too
        num = transfer_function(substitute(sys, vals), s)
        @test [at(num["B_out_p", "A_in_x"], w) for w in ω] ≈ Tω["B_out_p", "A_in_x"]
        part = transfer_function(substitute(sys, Dict(A_κ => 0.3, B_Δ => -0.2)), s)
        @test [at(part[:B_x, "A_in_p"], w) for w in ω] ≈ Tω[:B_x, "A_in_p"]

        # A system with no modes is its scattering matrix
        @variables θ
        splitter = SLH("bs", Subspace[], ["in1", "in2"], ["out1", "out2"], [cos(θ) im * sin(θ); im * sin(θ) cos(θ)], [0, 0], 0)
        bs = transfer_function(QuantumStateSpace(splitter), s)
        @test isequal(bs["out1_x", "in1_x"], cos(θ))

        # The symbolic method needs the quadrature basis, real matrices, and an s that is not a parameter
        @test_throws "quadrature basis" transfer_function(QuantumStateSpace(cavity("A"), AnnihilationCreationBasis()), s)
        @test_throws "is a parameter" transfer_function(QuantumStateSpace(cavity("A")), κ)
        hf = FockSpace(:c)
        @qnumbers a::Destroy(hf)
        @test_throws "complex entry" transfer_function(QuantumStateSpace(fock("nh", [1], [sqrt(κ) * a], Δ * a' * a')), s)

        # The numeric method needs every parameter, and folds constants without a substitute
        @test_throws "substitute values" transfer_function(sys, ω)
        lit = QuantumStateSpace(fock("lit", [1], [sqrt(0.3) * a], 0.7 * a' * a))
        @test transfer_function(lit, ω)["out_p", "in_x"] ≈ transfer_function(substitute(QuantumStateSpace(cavity("A")), Dict(κ => 0.3, Δ => 0.7)), ω)["out_p", "in_x"]
    end

    @testset "single-mode squeezing convention" begin
        hf = FockSpace(:sqz)
        @qnumbers a::Destroy(hf)
        @variables κ ϵ
        # Collett & Gardiner: threshold at ϵ = κ/2
        sqz = fock("sqz", [1], [sqrt(κ) * a], (1im / 2) * ϵ * (a'^2 - a^2))
        κ0, ϵ0 = 1.0, 0.3
        qss = substitute(QuantumStateSpace(sqz), Dict(κ => κ0, ϵ => ϵ0))

        # The ½ in H: a†² contributes 2c to Ω₊, so the gain is ϵ
        @test qss.A ≈ [ϵ0 - κ0 / 2 0; 0 -ϵ0 - κ0 / 2]

        ω = [0.1, 1.0, 3.0]
        # The vacuum spectrum of out_p, ½ in each input quadrature
        vacuum_p(sys, ω) = (T = transfer_function(sys, ω); (abs2.(T["out_p", "in_x"]) .+ abs2.(T["out_p", "in_p"])) ./ 2)
        Sp = vacuum_p(qss, ω)
        @test Sp ≈ [0.5 * (w^2 + (κ0 / 2 - ϵ0)^2) / (w^2 + (κ0 / 2 + ϵ0)^2) for w in ω]

        # The output stays a minimum-uncertainty state at every frequency
        G = freqresp(qss, ω)
        @test abs2.(G[1, 1, :]) .* abs2.(G[2, 2, :]) ≈ ones(length(ω))

        # At threshold the squeezed quadrature goes to zero noise as ω → 0
        thr = substitute(QuantumStateSpace(sqz), Dict(κ => κ0, ϵ => κ0 / 2))
        @test only(vacuum_p(thr, [1.0e-3])) < 1.0e-5

        # The component library follows the same convention
        lib = substitute(QuantumStateSpace(squeezing_cavity("A")), Dict(κ => κ0, ϵ => ϵ0))
        @test lib.A ≈ qss.A
        @test lib.B ≈ qss.B && lib.C ≈ qss.C
    end

    @testset "coefficient accumulation (C1)" begin
        hf = FockSpace(:cavity)
        @qnumbers a::Destroy(hf)
        @variables ω κ m
        κ0, ω0, m0 = 0.3, 0.7, 4.0
        vals = Dict(κ => κ0, ω => ω0, m => m0)

        # ω(a†a + aa†) = 2ω a†a + ω: the like terms merge and the constant is dropped
        sys = fock("r3", [1], [sqrt(κ) * a], ω * (a' * a + a * a'))
        @test substitute(QuantumStateSpace(sys, AnnihilationCreationBasis()), vals).A ≈ diagm([-2im * ω0 - κ0 / 2, 2im * ω0 - κ0 / 2])

        # A non-Hermitian H builds every row from i[H, ⋅]: H = ω a†² gives ȧ = -2iω a† and ȧ† = 0
        sys = fock("nh", [1], [sqrt(κ) * a], ω * a' * a')
        ladder = substitute(QuantumStateSpace(sys, AnnihilationCreationBasis()), vals)
        @test ladder.A ≈ [-κ0/2 -2im * ω0; 0 -κ0/2]
        @test substitute(QuantumStateSpace(QuantumStateSpace(sys, AnnihilationCreationBasis()), QuadratureBasis()), vals).A ≈ substitute(QuantumStateSpace(sys), vals).A

        # A free mass written in ladder operators: H = p²/2m has a†a, a†², a² and constant
        # terms, which must combine to ẋ = p/m, ṗ = 0 (it came out as a spring before C1)
        p = 1im * (a' - a) / sqrt(2)
        sys = fock("mass", [1], [sqrt(κ) * a], p * p / (2m))
        @test substitute(QuantumStateSpace(sys), vals).A ≈ [-κ0/2 1/m0; 0 -κ0/2]

        # L = √κ a + √κ a couples with strength 2√κ
        sys = fock("r5", [1], [sqrt(κ) * a + sqrt(κ) * a], ω * a' * a)
        @test substitute(QuantumStateSpace(sys, AnnihilationCreationBasis()), vals).C ≈ 2 * sqrt(κ0) * I(2)
    end

    @testset "complex parameters" begin
        hf = FockSpace(:cavity)
        @qnumbers a::Destroy(hf)
        @variables η::Number ω κ

        # concatenate keeps a complex parameter complex, so conj survives and H stays Hermitian
        drive = fock("d", [1], [sqrt(κ) * a], η * a' + conj(η) * a + ω * a' * a)
        promoted = concatenate("top", drive)
        @variables d_η::Number
        @test any(v -> isequal(v, d_η), Symbolics.get_variables(promoted.H))
        @test isequal(promoted.H, promoted.H')

        # A complex beamsplitter coupling gives the same numeric A built directly or concatenated
        bc = [GenericSubspace("b", FockSpace), GenericSubspace("c", FockSpace)]
        (b,), (c,) = operators(bc, "b"), operators(bc, "c")
        @variables g::Number
        bs = SLH(
            "bs", bc, ["i1", "i2"], ["o1", "o2"],
            [1 0; 0 1], [sqrt(κ) * b, sqrt(κ) * c], g * b' * c + conj(g) * b * c'
        )
        g0, κ0 = 0.3 + 0.4im, 1.0
        direct = substitute(QuantumStateSpace(bs, AnnihilationCreationBasis()), Dict(g => g0, κ => κ0)).A
        @variables bs_g::Number bs_κ
        viaconcat = substitute(
            QuantumStateSpace(concatenate("top", bs), AnnihilationCreationBasis()),
            Dict(bs_g => g0, bs_κ => κ0)
        ).A
        # ḃ = -ig c and ċ = -ig* b, and the rows for b† and c† carry the conjugates
        @test direct[1, 3] ≈ -im * g0 && direct[3, 1] ≈ -im * conj(g0)
        @test direct[2, 4] ≈ im * conj(g0) && direct[4, 2] ≈ im * g0
        @test viaconcat ≈ direct

        # The quadrature A of a physical system is real: ẋ_b = Im(g) x_c + Re(g) p_c, ...
        quad = substitute(QuantumStateSpace(bs), Dict(g => g0, κ => κ0)).A
        @test eltype(quad) == Float64
        gr, gi = real(g0), imag(g0)
        @test quad ≈ [
            -κ0/2 0 gi gr;
            0 -κ0/2 -gr gi;
            -gi gr -κ0/2 0;
            -gr -gi 0 -κ0/2
        ]
    end

    @testset "phase space" begin
        # Commutators are preserved: AΘ + ΘAᵀ + BΘBᵀ = 0 (James, Nurdin & Petersen 2008, Theorem 2.1)
        Θ(k) = kron(I(k), [0 1; -1 0])
        realizable(q) = isapprox(q.A * Θ(nstates(q) ÷ 2) + Θ(nstates(q) ÷ 2) * q.A' + q.B * Θ(ninputs(q) ÷ 2) * q.B', zeros(nstates(q), nstates(q)); atol = 1.0e-12)

        # The quadrature builder agrees with the annihilation-creation builder followed by conversion
        c, = operators([OpticalMode("")], "")
        (a,), (b,) = operators(TWOMODES, "a"), operators(TWOMODES, "b")
        @variables κ Δ ϵ ω g::Number A_κ A_Δ B_κ B_Δ
        vals = Dict(κ => 0.3, Δ => 0.7, ϵ => 0.1, ω => 1.3, g => 0.3 + 0.4im, A_κ => 0.3, A_Δ => 0.7, B_κ => 0.5, B_Δ => -0.2)
        systems = [
            cavity("A"),
            squeezing_cavity("A"),
            feedbackreduce(concatenate("S", cavity("A"), cavity("B")), "A_out" => "B_in"),
            fock("sqz", [1], [sqrt(κ) * c], (1im / 2) * ϵ * (c'^2 - c^2)),
            fock2("tms", [1 0; 0 1], [sqrt(κ) * a, sqrt(κ) * b], 1im * ϵ * (a' * b' - a * b)),
            fock2("bs", [1 0; 0 1], [sqrt(κ) * a, sqrt(κ) * b], g * a' * b + conj(g) * a * b'),
            fock2("om", [1], [sqrt(κ) * a], ω * b' * b - ϵ * (b' + b) * (a' + a)),
            fock("scatter", [0.6 0.8im; 0.8im 0.6], [sqrt(κ) * c, sqrt(Δ) * c'], ω * c' * c),
        ]
        for sys in systems
            direct = substitute(QuantumStateSpace(sys), vals)
            converted = substitute(QuantumStateSpace(QuantumStateSpace(sys, AnnihilationCreationBasis()), QuadratureBasis()), vals)
            @test all(ssdata(direct) .≈ ssdata(converted))
            @test realizable(direct)
        end

        # Converting to the same basis returns the system, and the reverse conversion is not implemented
        q = QuantumStateSpace(cavity("A"))
        @test QuantumStateSpace(q, QuadratureBasis()) === q
        @test_throws "is not implemented" QuantumStateSpace(q, AnnihilationCreationBasis())

        # The free mass H = p²/2m, L = iλx: x is read out, and the force kicks p
        ps = PhaseSpace(:mass)
        x, p = Position(ps, :x), Momentum(ps, :p)
        @variables m λ
        m0, λ0, Ω = 2.0, 0.7, [0.5, 3.0]
        mass(ω0) = substitute(
            QuantumStateSpace(SLH("mass", [MechanicalMode("", PhaseSpace)], ["in"], ["out"], [1], [1im * λ * x], p * p / (2m) + m * ω^2 * x * x / 2)),
            Dict(m => m0, λ => λ0, ω => ω0)
        )
        free = mass(0.0)
        @test free.A ≈ [0 1 / m0; 0 0]
        @test realizable(free)
        G = freqresp(free, Ω)
        @test G[:, :, 1] ≈ [1 0; -2λ0^2 / (m0 * Ω[1]^2) 1]
        @test G[2, 1, :] ≈ -2λ0^2 ./ (m0 .* Ω .^ 2)

        # A harmonic potential moves the response smoothly away from the free mass
        @test freqresp(mass(1.3), Ω)[2, 1, :] ≈ 2λ0^2 ./ (m0 .* (1.3^2 .- Ω .^ 2))

        # Ladder and position/momentum operators mix within one system
        subs = [OpticalMode("cavity"), MechanicalMode("mirror", PhaseSpace)]
        (c,), (X, P) = operators(subs, "cavity"), operators(subs, "mirror")
        om = SLH("om", subs, ["in"], ["out"], [1], [sqrt(κ) * c], P * P / (2m) - ϵ * X * (c' + c))
        q = substitute(QuantumStateSpace(om), Dict(κ => 0.3, ϵ => 0.1, m => m0))
        @test state_names(q) == [:cavity_x, :cavity_p, :mirror_x, :mirror_p]
        @test q.A ≈ [-0.15 0 0 0; 0 -0.15 sqrt(2) * 0.1 0; 0 0 0 1 / m0; sqrt(2) * 0.1 0 0 0]
        @test realizable(q)

        # Without a length scale there is no annihilation-creation basis; linear H and constant L are rejected
        @test_throws ErrorException QuantumStateSpace(om, AnnihilationCreationBasis())
        @test_throws ErrorException QuantumStateSpace(SLH("lin", [MechanicalMode("", PhaseSpace)], ["in"], ["out"], [1], [1im * λ * x], λ * x))
        @test_throws ErrorException QuantumStateSpace(SLH("drive", [MechanicalMode("", PhaseSpace)], ["in"], ["out"], [1], [λ + 1im * λ * x], p * p / (2m)))
    end

    @testset "state-space composition" begin
        # Building the state space commutes with concatenate and feedbackreduce, names included
        @variables κ Δ θ ϵ g::Number μ
        @variables A_κ A_Δ B_κ B_Δ Q_κ Q_ϵ F_κ F_Δ bs_θ bs_κ bs_g::Number arm_κ arm_μ arm_ϵ scatter_κ scatter_Δ scatter_ϵ
        vals = Dict(
            A_κ => 0.3, A_Δ => 0.7, B_κ => 0.5, B_Δ => -0.2, Q_κ => 1.0, Q_ϵ => 0.3, F_κ => 0.2, F_Δ => 0.1,
            bs_θ => 0.3, bs_κ => 0.4, bs_g => 0.3 + 0.4im, arm_κ => 0.4, arm_μ => 2.0, arm_ϵ => 0.2,
            scatter_κ => 0.4, scatter_Δ => 0.6, scatter_ϵ => 0.1,
        )
        same(P, Q) = all(ssdata(substitute(P, vals)) .≈ ssdata(substitute(Q, vals))) &&
            P.inputs == Q.inputs && P.outputs == Q.outputs && state_names(P) == state_names(Q)

        (a,), (b,) = operators(TWOMODES, "a"), operators(TWOMODES, "b")
        cbs = fock2("bs", [1 0; 0 1], [sqrt(κ) * a, sqrt(κ) * b], g * a' * b + conj(g) * a * b')
        splitter = SLH("bs", Subspace[], ["in1", "in2"], ["out1", "out2"], [cos(θ) im * sin(θ); im * sin(θ) cos(θ)], [0, 0], 0)
        subs = [OpticalMode("cavity"), MechanicalMode("mirror", PhaseSpace)]
        (c,), (x, p) = operators(subs, "cavity"), operators(subs, "mirror")
        arm = SLH("arm", subs, ["in"], ["out"], [1], [sqrt(κ) * c], p^2 / (2μ) - ϵ * x * (c' + c))

        # (systems, connections) in both bases, and position and momentum in the quadrature basis only
        networks = [
            ([cavity("A"), cavity("B")], ["A_out" => "B_in"], (QuadratureBasis(), AnnihilationCreationBasis())),
            ([squeezing_cavity("Q"), cavity("F")], ["Q_out" => "F_in"], (QuadratureBasis(), AnnihilationCreationBasis())),
            ([cbs, cavity("A")], ["bs_out1" => "A_in", "A_out" => "bs_in2"], (QuadratureBasis(), AnnihilationCreationBasis())),
            ([cavity("A"), splitter], ["A_out" => "bs_in2", "bs_out2" => "A_in"], (QuadratureBasis(), AnnihilationCreationBasis())),
            ([cavity("F"), arm], ["F_out" => "arm_in"], (QuadratureBasis(),)),
        ]
        for (systems, connections, bases) in networks, basis in bases
            slh = concatenate("net", systems...)
            ss_ = concatenate("net", (QuantumStateSpace(sys, basis) for sys in systems)...)
            @test same(QuantumStateSpace(slh, basis), ss_)
            for connection in connections
                slh = feedbackreduce(slh, connection)
                ss_ = feedbackreduce(ss_, connection)
                @test same(QuantumStateSpace(slh, basis), ss_)
            end
            # several connections in one call are made in order
            @test same(QuantumStateSpace(feedbackreduce(concatenate("net", systems...), connections...), basis), ss_)
        end

        # A single system with a non-diagonal S, connected to itself
        d, = operators([OpticalMode("")], "")
        scatter = concatenate("top", fock("scatter", [0.6 0.8im; 0.8im 0.6], [sqrt(κ) * d, sqrt(Δ) * d'], ϵ * d' * d))
        for basis in (QuadratureBasis(), AnnihilationCreationBasis())
            @test same(QuantumStateSpace(feedbackreduce(scatter, "scatter_out1" => "scatter_in2"), basis),
                feedbackreduce(QuantumStateSpace(scatter, basis), "scatter_out1" => "scatter_in2"))
        end

        # The same errors as for SLH systems
        q = QuantumStateSpace(cavity("A"))
        @test_throws "unit gain" feedbackreduce(q, "out" => "in")
        @test_throws "the order is output => input" feedbackreduce(q, "in" => "out")
        @test_throws "unique names" concatenate("top", q, q)
        @test_throws "share a basis" concatenate("top", q, QuantumStateSpace(cavity("B"), AnnihilationCreationBasis()))
        @test_throws "cannot concatenate SLH and QuantumStateSpace" concatenate("top", cavity("B"), q)
        @test_throws "needs at least one system" concatenate("top")

        # ControlSystemsBase's compositions point to these, and to ss for plain state spaces
        @test_throws "concatenate and feedbackreduce" q * q
        @test_throws "concatenate and feedbackreduce" q + q
        @test_throws "concatenate and feedbackreduce" ControlSystemsBase.append(q, q)
    end

    @testset "two-mode squeezing convention" begin
        # Control case for the ½: an off-diagonal a†b† term is handled correctly
        (a,), (b,) = operators(TWOMODES, "a"), operators(TWOMODES, "b")
        @variables κ ϵ
        sys = fock2("tms", [1 0; 0 1], [sqrt(κ) * a, sqrt(κ) * b], 1im * ϵ * (a' * b' - a * b))
        κ0, ϵ0 = 1.0, 0.3
        A = substitute(QuantumStateSpace(sys, AnnihilationCreationBasis()), Dict(κ => κ0, ϵ => ϵ0)).A

        # state order (a, a†, b, b†): ȧ = ϵ b†, ȧ† = ϵ b, ...
        @test A ≈ [
            -κ0 / 2 0 0 ϵ0;
            0 -κ0 / 2 ϵ0 0;
            0 ϵ0 -κ0 / 2 0;
            ϵ0 0 0 -κ0 / 2
        ]
    end

end
