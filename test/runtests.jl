using SLHQuantumSystems
using Test
using SecondQuantizedAlgebra
using Symbolics
using LinearAlgebra
using PhysicalConstants.CODATA2018: ReducedPlanckConstant as ℏ_SI

@testset "SLHQuantumSystems.jl" begin

    @testset "get_additive_terms" begin
        # Create test Hilbert space and operators
        hf = FockSpace(:cavity)
        @qnumbers a::Destroy(hf)

        # Create symbolic parameters
        @cnumbers g ω κ

        # Single term (no addition)
        single_term = g * a
        terms = get_additive_terms(single_term)
        @test length(terms) == 1
        @test isequal(terms[1], g * a)
        @test isequal(sum(terms), single_term)

        # Simple addition of two terms
        two_terms = g * a + ω * a'
        terms = get_additive_terms(two_terms)
        @test length(terms) == 2
        @test isequal(sum(terms), two_terms)

        # Multiple terms with different operators
        multi_terms = g * a + ω * a' + κ * a' * a
        terms = get_additive_terms(multi_terms)
        @test length(terms) == 3
        @test isequal(sum(terms), multi_terms)

        # Nested additions
        nested = (g * a + ω * a') + κ * a' * a
        terms = get_additive_terms(nested)
        @test length(terms) == 3
        @test isequal(sum(terms), nested)

        # Single number
        single_num = 5
        terms = get_additive_terms(single_num)
        @test length(terms) == 1
        @test terms[1] == 5
        @test isequal(sum(terms), single_num)

        # Complex multiplication term (commutes to two terms due to [a,a']=1)
        complex_mult = g * ω * a * a'
        terms = get_additive_terms(complex_mult)
        @test length(terms) == 2  # g*ω*(a′*a) + g*ω due to commutation
        @test isequal(sum(terms), complex_mult)

        # Zero term
        zero_term = 0 * a
        terms = get_additive_terms(zero_term)
        @test length(terms) == 1
        @test isequal(sum(terms), zero_term)
    end

    @testset "get_qnumbers" begin
        # Create test Hilbert space and operators
        hf = FockSpace(:cavity)
        @qnumbers a::Destroy(hf)

        # Create symbolic parameters
        @cnumbers g ω κ

        # Single operator
        qsyms = get_qnumbers(a)
        @test length(qsyms) == 1
        @test a in qsyms

        # Multiple operators
        expr = g * a + ω * a'
        qsyms = get_qnumbers(expr)
        @test length(qsyms) == 2
        @test a in qsyms
        @test a' in qsyms

        # Just parameters (no quantum operators)
        param_expr = g + ω * κ
        qsyms = get_qnumbers(param_expr)
        @test length(qsyms) == 0

        # Complex expression
        complex_expr = g * a' * a + κ * a
        qsyms = get_qnumbers(complex_expr)
        @test length(qsyms) == 2
        @test a in qsyms
        @test a' in qsyms
    end

    @testset "get_cnumbers" begin
        # Create test Hilbert space and operators
        hf = FockSpace(:cavity)
        @qnumbers a::Destroy(hf)

        # Create symbolic parameters
        @cnumbers g ω κ

        # Single parameter
        numsyms = get_cnumbers(g)
        @test length(numsyms) == 1
        @test g in numsyms

        # Multiple parameters
        expr = g * ω + κ
        numsyms = get_cnumbers(expr)
        @test length(numsyms) == 3
        @test g in numsyms
        @test ω in numsyms
        @test κ in numsyms

        # Parameters with operators
        mixed_expr = g * a + ω * a'
        numsyms = get_cnumbers(mixed_expr)
        @test length(numsyms) == 2
        @test g in numsyms
        @test ω in numsyms

        # Just operators (no parameters)
        op_expr = a + a'
        numsyms = get_cnumbers(op_expr)
        @test length(numsyms) == 0
    end


    @testset "SLH struct and basic operations" begin
        # Create test system
        hf = FockSpace(:cavity)
        @qnumbers a::Destroy(hf)
        @cnumbers ω κ g

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

    @testset "cascade" begin
        # Create two test systems (following the docs quick start pattern)
        hilb = FockSpace("cavity")
        @qnumbers a::Destroy(hilb)
        @cnumbers ω κ

        # Define system components
        H = ω * a' * a
        L = [√κ * a]
        S = [1]

        cavityA = SLH("A", S, L, H)
        cavityB = SLH("B", S, L, H)

        # Test that concatenate runs without error
        combined = concatenate([cavityA, cavityB], "chain")
        @test isa(combined, SLH)

        # Test that feedbackreduce runs without error
        cascaded = feedbackreduce(combined, "A_out", "B_in")
        @test isa(cascaded, SLH)
    end

    @testset "optomechanical quadrature" begin
        # Builds the optomechanical SLH from examples/calibrated_strain.jl and
        # verifies that, with g = 0, every diagonal output spectral density
        # equals the vacuum value 1/2. This is the vacuum shot-noise sanity
        # check — it exercises the unitless mechanical quadrature convention
        # end to end (A matrix in rad/s, (1/√2) transform for both modes).
        hilb = FockSpace(:cavity) ⊗ FockSpace(:mirror)
        subspaces = [OpticalMode(""), MechanicalMode("")]

        opt_sub, mech_sub = subspaces
        a = Destroy(hilb, operatornames(opt_sub)[1], 1)
        b = Destroy(hilb, operatornames(mech_sub)[1], 2)

        @variables ω l κ Ω m Γ g
        H = Ω * b' * b - g * (b' + b) * (a' + a)
        L_ops = [κ * a, Γ * b]
        S_mat = [1 0; 0 1]

        pdict = Dict(zip(nameof.([ω, l, κ, Ω, m, Γ, g]), [ω, l, κ, Ω, m, Γ, g]))
        opdict = Dict(zip(getfield.([a, b], :name), [a, b]))

        slh = SLH(
            "opto", subspaces, pdict, opdict,
            ["l_in", "m_in"], ["l_out", "m_out"], S_mat, L_ops, H
        )

        qss = toquadrature(QuantumStateSpace(slh))

        # LIGO-scale numerics; κ carries √ because L_op = κ·a gives amplitude
        # decay rate κ²/2, matching the calibrated_strain.jl convention.
        L_arm = 3995.0
        T_ITM = 0.014
        c_phys = 2.99792458e8
        κ_cavity = T_ITM * c_phys / (2 * L_arm)
        paramdict = Dict(
            ω => 0.0,
            l => L_arm,
            κ => sqrt(κ_cavity),
            Ω => 2π * 5.0e-6,
            m => 39.6 / 2,
            g => 0.0,
            Γ => 1.0e-23,
        )
        numeric = substitute(qss, paramdict)

        # A matrix must be finite with rad/s-scale entries (no m·Ω factors).
        A = ssdata(numeric)[1]
        @test all(isfinite, A)

        # Vacuum shot-noise: diag(S_out) ≡ 1/2 at every frequency.
        freqs = 2π .* [10.0, 100.0, 1.0e3, 1.0e4]
        sd = spectral_density(numeric, freqs)
        for i in 1:size(sd.S, 1), k in 1:length(freqs)
            @test isapprox(real(sd.S[i, i, k]), 0.5; atol = 1.0e-8)
        end
    end

    @testset "quadrature transform with complex ladder entries" begin
        # A squeezer cascaded into a filter cavity gives a ladder-basis A matrix
        # whose entries are Complex{Num} with a nonzero symbolic imaginary part
        # (the detuning terms). Symbolics.expand cannot handle those directly,
        # so toquadrature has to flatten them first.
        sys = feedbackreduce(
            concatenate([squeezing_cavity("A"), cavity("B")], "sys"),
            "A_out", "B_in"
        )

        qss = toquadrature(QuantumStateSpace(sys))
        @test isa(qss, QuantumStateSpace)

        # The quadrature-basis A matrix must come out real once evaluated.
        params = qss.parameters
        numeric = substitute(
            qss, Dict(
                params[:A_κ] => 1.0, params[:A_ϵ] => 0.1,
                params[:B_B_κ] => 2.0, params[:B_B_Δ] => 3.0,
                params[:B_B_ω] => 0.0, params[:B_B_l] => 1.0
            )
        )
        A = ssdata(numeric)[1]
        @test all(a -> isapprox(imag(Symbolics.value(a)), 0.0; atol = 1.0e-12), A)
    end

end
