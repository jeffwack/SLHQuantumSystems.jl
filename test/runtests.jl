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
