using SLHQuantumSystems
using Test
using SecondQuantizedAlgebra
using Symbolics
using LinearAlgebra
using ControlSystems

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
        sys = SLH("spare", [GenericMode(""), GenericMode("")], ["in"], ["out"], [1], [κ * a], ω * a' * a)
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
        @test length(SLHQuantumSystems.fresponse_state2output(qss, [0.1, 1.0], 1, 1)) == 2
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
        Sp = real(spectral_density(qss, ω)["out_p", "out_p"])
        @test Sp ≈ [0.5 * (w^2 + (κ0 / 2 - ϵ0)^2) / (w^2 + (κ0 / 2 + ϵ0)^2) for w in ω]

        # The output stays a minimum-uncertainty state at every frequency
        G = freqresp(qss, ω)
        @test abs2.(G[1, 1, :]) .* abs2.(G[2, 2, :]) ≈ ones(length(ω))

        # At threshold the squeezed quadrature goes to zero noise as ω → 0
        thr = substitute(QuantumStateSpace(sqz), Dict(κ => κ0, ϵ => κ0 / 2))
        @test real(spectral_density(thr, [1.0e-3])["out_p", "out_p"][1]) < 1.0e-5

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
            "bs", [GenericMode(""), GenericMode("")], ["i1", "i2"], ["o1", "o2"],
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
