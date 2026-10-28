using SLHQuantumSystems
using Test
using SecondQuantizedAlgebra
using Symbolics
using LinearAlgebra
using ControlSystemsBase

@testset "SLHQuantumSystems.jl" begin

    @testset "SLH constructor" begin
        # Create test system
        hf = FockSpace(:cavity)
        @qnumbers a::Destroy(hf)
        @variables ω κ g

        H = ω * a' * a
        L = [√κ * a]
        S = [1]
        sys = SLH("test", S, L, H)

        # SLH construction
        @test sys.name == "test"
        @test sys.inputs == ["in"]
        @test sys.outputs == ["out"]
        @test sys.S == [1]
        @test length(sys.L) == 1

    end

    @testset "SLH composition" begin
        cavityA = cavity("A")
        cavityB = cavity("B")

        # Test that concatenate runs without error
        combined = concatenate([cavityA, cavityB], "S")
        @test isa(combined, SLH)

        # Test that feedbackreduce runs without error
        cascaded = feedbackreduce(combined, "A_out", "B_in")
        @test isa(cascaded, SLH)
    end

    @testset "SLH to ABCD" begin
        qss = QuantumStateSpace(cavity("A"))
        @test isa(qss, QuantumStateSpace)
    end

    @testset "operator validation" begin
        hf = FockSpace(:cavity)
        @qnumbers a::Destroy(hf)
        @variables ω κ λ m

        # PhaseSpace operators are rejected by the ladder parser, in H and in L
        ps = PhaseSpace(:mass)
        x = Position(ps, :x)
        p = Momentum(ps, :p)
        @test_throws ErrorException QuantumStateSpace(SLH("mass", [1], [κ * a], p * p / (2m)), LadderBasis())
        @test_throws ErrorException QuantumStateSpace(SLH("mass", [1], [1im * λ * x], ω * a' * a), LadderBasis())

        # Every operator must act on a declared subspace
        hp = FockSpace(:one) ⊗ FockSpace(:two)
        b = Destroy(hp, :b, 2)
        @test_throws ErrorException SLH("bad", [GenericMode("")], ["in"], ["out"], [1], [κ * b], ω * b' * b)

        # A declared subspace that no operator touches is allowed, and is decoupled
        sys = SLH("spare", [GenericMode(""), GenericMode("spare")], ["in"], ["out"], [1], [κ * a], ω * a' * a)
        qss = QuantumStateSpace(sys, LadderBasis())
        @test size(qss.A) == (4, 4)
        @test all(iszero, qss.A[3:4, :]) && all(iszero, qss.A[:, 3:4])
        @test all(iszero, qss.B[3:4, :])
        @test all(iszero, qss.C[:, 3:4])

        # A constant in L is a coherent drive, which the linear model rejects; a zero constant is fine
        @variables β
        @test_throws ErrorException QuantumStateSpace(SLH("drive", [1], [β + κ * a], ω * a' * a), LadderBasis())
        @test_throws ErrorException QuantumStateSpace(
            SLH("drive", [GenericMode("")], ["i1", "i2"], ["o1", "o2"], [1 0; 0 1], [κ * a, β], ω * a' * a),
            LadderBasis()
        )
        @test isa(QuantumStateSpace(SLH("nodrive", [1], [0 + κ * a], ω * a' * a), LadderBasis()), QuantumStateSpace)

        # A bare operator is coerced to an iterable QAdd
        bare = SLH("t", [1], [a], ω * a' * a)
        @test all(Lj -> Lj isa SecondQuantizedAlgebra.QAdd, bare.L)
        @test isequal(
            QuantumStateSpace(bare, LadderBasis()).A,
            QuantumStateSpace(SLH("t", [1], [1 * a], ω * a' * a), LadderBasis()).A
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
        @test mode_names(cavity("A")) == [""]
        @test state_names(QuantumStateSpace(cavity("A"))) == [:x, :p]
        @test state_names(QuantumStateSpace(cavity("A"), LadderBasis())) == [:a, :adag]
        both = concatenate([cavity("A"), cavity("B")], "S")
        @test mode_names(both) == ["A", "B"]
        @test state_names(QuantumStateSpace(both)) == [:A_x, :A_p, :B_x, :B_p]

        # input_names and output_names label the columns and rows of D, two per port, as
        # ControlSystemsBase expects; input_ports and output_ports name the ports alone
        qboth = QuantumStateSpace(both)
        @test input_ports(both) == input_ports(qboth) == ["A_in", "B_in"]
        @test output_ports(qboth) == ["A_out", "B_out"]
        @test input_names(qboth) == ["A_in_x", "A_in_p", "B_in_x", "B_in_p"]
        @test output_names(QuantumStateSpace(both, LadderBasis())) == ["A_out_a", "A_out_adag", "B_out_a", "B_out_adag"]
        @test length(input_names(qboth)) == ninputs(qboth)
        @test length(output_names(qboth)) == noutputs(qboth)
        @test input_names(qboth, 3) == "B_in_x"
        @test output_names(qboth, 2) == "A_out_p"
        @test state_names(qboth, 4) == :B_p

        # Mode names must be unique, which catches concatenating two systems of the same name
        @test_throws ErrorException concatenate([cavity("A"), cavity("A")], "S")
        hp = FockSpace(:one) ⊗ FockSpace(:two)
        b, c = Destroy(hp, :b, 1), Destroy(hp, :c, 2)
        @test_throws ErrorException SLH("dup", [GenericMode(""), GenericMode("")], ["in"], ["out"], [1], [κ * b], Δ * c' * c)
        @test mode_names(SLH("tms", [1 0; 0 1], [κ * b, κ * c], Δ * b' * c')) == ["mode1", "mode2"]

        # Strings name ports and Symbols name states; a name without its suffix selects both entries
        @variables A_κ A_Δ B_κ B_Δ
        qss = substitute(QuantumStateSpace(both), Dict(A_κ => κ0, A_Δ => Δ0, B_κ => κ0, B_Δ => Δ0))
        ω = [0.1, 1.0, 3.0]
        T = transfer_matrix(qss, ω)
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

    @testset "single-mode squeezing convention" begin
        hf = FockSpace(:sqz)
        @qnumbers a::Destroy(hf)
        @variables κ ϵ
        # Collett & Gardiner: threshold at ϵ = κ/2
        sqz = SLH("sqz", [1], [sqrt(κ) * a], (1im / 2) * ϵ * (a'^2 - a^2))
        κ0, ϵ0 = 1.0, 0.3
        qss = substitute(QuantumStateSpace(sqz), Dict(κ => κ0, ϵ => ϵ0))

        # The ½ in H: a†² contributes 2c to Ω₊, so the gain is ϵ
        @test qss.A ≈ [ϵ0 - κ0 / 2 0; 0 -ϵ0 - κ0 / 2]

        ω = [0.1, 1.0, 3.0]
        # The vacuum spectrum of out_p, ½ in each input quadrature
        vacuum_p(sys, ω) = (T = transfer_matrix(sys, ω); (abs2.(T["out_p", "in_x"]) .+ abs2.(T["out_p", "in_p"])) ./ 2)
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
        sys = SLH("r3", [1], [sqrt(κ) * a], ω * (a' * a + a * a'))
        @test substitute(QuantumStateSpace(sys, LadderBasis()), vals).A ≈ diagm([-2im * ω0 - κ0 / 2, 2im * ω0 - κ0 / 2])

        # A free mass written in ladder operators: H = p²/2m has a†a, a†², a² and constant
        # terms, which must combine to ẋ = p/m, ṗ = 0 (it came out as a spring before C1)
        p = 1im * (a' - a) / sqrt(2)
        sys = SLH("mass", [1], [sqrt(κ) * a], p * p / (2m))
        @test substitute(QuantumStateSpace(sys), vals).A ≈ [-κ0/2 1/m0; 0 -κ0/2]

        # L = √κ a + √κ a couples with strength 2√κ
        sys = SLH("r5", [1], [sqrt(κ) * a + sqrt(κ) * a], ω * a' * a)
        @test substitute(QuantumStateSpace(sys, LadderBasis()), vals).C ≈ 2 * sqrt(κ0) * I(2)
    end

    @testset "complex parameters" begin
        hf = FockSpace(:cavity)
        @qnumbers a::Destroy(hf)
        @variables η::Number ω κ

        # concatenate keeps a complex parameter complex, so conj survives and H stays Hermitian
        drive = SLH("d", [1], [sqrt(κ) * a], η * a' + conj(η) * a + ω * a' * a)
        promoted = concatenate([drive], "top")
        @variables d_η::Number
        @test any(v -> isequal(v, d_η), Symbolics.get_variables(promoted.H))
        @test isequal(promoted.H, promoted.H')

        # A complex beamsplitter coupling gives the same numeric A built directly or concatenated
        hp = FockSpace(:one) ⊗ FockSpace(:two)
        b = Destroy(hp, :b, 1)
        c = Destroy(hp, :c, 2)
        @variables g::Number
        bs = SLH(
            "bs", [GenericMode("b"), GenericMode("c")], ["i1", "i2"], ["o1", "o2"],
            [1 0; 0 1], [sqrt(κ) * b, sqrt(κ) * c], g * b' * c + conj(g) * b * c'
        )
        g0, κ0 = 0.3 + 0.4im, 1.0
        direct = substitute(QuantumStateSpace(bs, LadderBasis()), Dict(g => g0, κ => κ0)).A
        @variables bs_g::Number bs_κ
        viaconcat = substitute(
            QuantumStateSpace(concatenate([bs], "top"), LadderBasis()),
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

        # The quadrature builder agrees with the ladder builder followed by toquadrature
        hp = FockSpace(:one) ⊗ FockSpace(:two)
        a, b = Destroy(hp, :a, 1), Destroy(hp, :b, 2)
        @variables κ Δ ϵ ω g::Number A_κ A_Δ B_κ B_Δ
        vals = Dict(κ => 0.3, Δ => 0.7, ϵ => 0.1, ω => 1.3, g => 0.3 + 0.4im, A_κ => 0.3, A_Δ => 0.7, B_κ => 0.5, B_Δ => -0.2)
        systems = [
            cavity("A"),
            squeezing_cavity("A"),
            feedbackreduce(concatenate([cavity("A"), cavity("B")], "S"), "A_out", "B_in"),
            SLH("sqz", [1], [sqrt(κ) * a], (1im / 2) * ϵ * (a'^2 - a^2)),
            SLH("tms", [1 0; 0 1], [sqrt(κ) * a, sqrt(κ) * b], 1im * ϵ * (a' * b' - a * b)),
            SLH("bs", [1 0; 0 1], [sqrt(κ) * a, sqrt(κ) * b], g * a' * b + conj(g) * a * b'),
            SLH("om", [1], [sqrt(κ) * a], ω * b' * b - ϵ * (b' + b) * (a' + a)),
            SLH("scatter", [0.6 0.8im; 0.8im 0.6], [sqrt(κ) * a, sqrt(Δ) * a'], ω * a' * a),
        ]
        for sys in systems
            direct = substitute(QuantumStateSpace(sys), vals)
            vialadder = substitute(toquadrature(QuantumStateSpace(sys, LadderBasis())), vals)
            @test all(ssdata(direct) .≈ ssdata(vialadder))
            @test realizable(direct)
        end

        # The free mass H = p²/2m, L = iλx: x is read out, and the force kicks p
        ps = PhaseSpace(:mass)
        x, p = Position(ps, :x), Momentum(ps, :p)
        @variables m λ
        m0, λ0, Ω = 2.0, 0.7, [0.5, 3.0]
        mass(ω0) = substitute(
            QuantumStateSpace(SLH("mass", [MechanicalMode("")], ["in"], ["out"], [1], [1im * λ * x], p * p / (2m) + m * ω^2 * x * x / 2)),
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
        hm = FockSpace(:cavity) ⊗ PhaseSpace(:mirror)
        c, X, P = Destroy(hm, :a, 1), Position(hm, :x, 2), Momentum(hm, :p, 2)
        om = SLH("om", [OpticalMode("cavity"), MechanicalMode("mirror")], ["in"], ["out"], [1], [sqrt(κ) * c], P * P / (2m) - ϵ * X * (c' + c))
        q = substitute(QuantumStateSpace(om), Dict(κ => 0.3, ϵ => 0.1, m => m0))
        @test state_names(q) == [:cavity_x, :cavity_p, :mirror_x, :mirror_p]
        @test q.A ≈ [-0.15 0 0 0; 0 -0.15 sqrt(2) * 0.1 0; 0 0 0 1 / m0; sqrt(2) * 0.1 0 0 0]
        @test realizable(q)

        # Without a length scale there is no ladder basis; linear H and constant L are rejected
        @test_throws ErrorException QuantumStateSpace(om, LadderBasis())
        @test_throws ErrorException QuantumStateSpace(SLH("lin", [MechanicalMode("")], ["in"], ["out"], [1], [1im * λ * x], λ * x))
        @test_throws ErrorException QuantumStateSpace(SLH("drive", [MechanicalMode("")], ["in"], ["out"], [1], [λ + 1im * λ * x], p * p / (2m)))
    end

    @testset "two-mode squeezing convention" begin
        # Control case for the ½: an off-diagonal a†b† term is handled correctly
        hp = FockSpace(:one) ⊗ FockSpace(:two)
        a = Destroy(hp, :a, 1)
        b = Destroy(hp, :b, 2)
        @variables κ ϵ
        sys = SLH("tms", [1 0; 0 1], [sqrt(κ) * a, sqrt(κ) * b], 1im * ϵ * (a' * b' - a * b))
        κ0, ϵ0 = 1.0, 0.3
        A = substitute(QuantumStateSpace(sys, LadderBasis()), Dict(κ => κ0, ϵ => ϵ0)).A

        # state order (a, a†, b, b†): ȧ = ϵ b†, ȧ† = ϵ b, ...
        @test A ≈ [
            -κ0 / 2 0 0 ϵ0;
            0 -κ0 / 2 ϵ0 0;
            0 ϵ0 -κ0 / 2 0;
            ϵ0 0 0 -κ0 / 2
        ]
    end

end
