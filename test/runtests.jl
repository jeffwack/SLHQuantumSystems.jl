using SLHQuantumSystems
using Test
using SecondQuantizedAlgebra
using Symbolics

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

end
