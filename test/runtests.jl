# Unit tests of MuRAT.jl. The comparison with MATLAB MuRAT 4.0 on the MuRAT test
# data sets is summarised in the README (Validation).
using MuRAT
using Test
using Statistics

@testset "MuRAT.jl" begin
    @testset "MATLAB colon operator" begin
        v = MuRAT.matlab_colon(1e-5, 1e-5, 1e-1)
        @test length(v) == 10000
        @test v[1] == 1e-5 && v[end] == 0.1
        @test MuRAT.matlab_colon(0.5, 1, 3.2) == [0.5, 1.5, 2.5]
        @test MuRAT.matlab_colon(3, -1.5, 0) == [3.0, 1.5, 0.0]
        @test isempty(MuRAT.matlab_colon(1, 1, 0))
    end

    @testset "smooth (moving average)" begin
        y = collect(1.0:10.0) .^ 2
        s = MuRAT.smooth_moving(y, 5)
        @test s[1] == y[1]
        @test s[2] ≈ mean(y[1:3])
        @test s[5] ≈ mean(y[3:7])
        @test s[end] == y[end]
        @test MuRAT.smooth_moving(y, 6) == s          # even spans are reduced to odd
    end

    @testset "grid formats" begin
        x, y, z = [1.0, 2.0, 3.0], [10.0, 20.0], [0.0, -5.0, -10.0, -15.0]
        r = MuRAT.unfold_xyz(x, y, z)
        @test size(r) == (24, 3)
        @test r[1, :] == [1, 10, 0] && r[2, :] == [1, 10, -5] && r[5, :] == [1, 20, 0]
        v = collect(1.0:24.0)
        V = MuRAT.fold_field(x, y, z, v)
        @test size(V) == (2, 3, 4)
        @test MuRAT.unfold_field(V) == v
        X, Y, Z = MuRAT.meshgrid3(x, y, z)
        @test MuRAT.unfold_mesh(X, Y, Z, V)[:, 4] == v
    end

    @testset "geodesy" begin
        d, az = MuRAT.geodistance(45.9805, -122.503, 46.39, -122.0)
        @test isapprox(d, 59832.39, atol = 0.01)
        lat, lon = MuRAT.georeckon(45.9805, -122.503, d, az)
        @test isapprox(lat, 46.39, atol = 1e-9) && isapprox(lon, -122.0, atol = 1e-9)
        @test MuRAT.geodistance(10.0, 20.0, 10.0, 20.0) == (0.0, 0.0)
    end

    @testset "interpolation" begin
        x, y, z = [0.0, 1.0, 2.0], [0.0, 1.0], [2.0, 1.0, 0.0]       # decreasing z
        X, Y, Z = MuRAT.meshgrid3(x, y, z)
        V = X .+ 2 .* Y .+ 3 .* Z
        @test MuRAT.interp3_linear(x, y, z, V, 0.5, 0.25, 1.5) ≈ 0.5 + 0.5 + 4.5
        @test isnan(MuRAT.interp3_linear(x, y, z, V, 3.0, 0.0, 0.0))
        K = [1.0 2.0; 3.0 4.0]
        @test MuRAT.interp2_linear_clamped([0.0, 1.0], [0.0, 1.0], K, 0.5, 0.5) ≈ 2.5
        @test MuRAT.interp2_linear_clamped([0.0, 1.0], [0.0, 1.0], K, -1.0, 0.8) == 2.0  # nearest node
        # fillmissing(...,'nearest') results checked in MATLAB R2025b
        @test MuRAT.fill_nearest!([1.0, NaN, 3.0]) == [1.0, 3.0, 3.0]
        @test MuRAT.fill_nearest!([5.0, NaN, NaN, NaN, 9.0]) == [5.0, 5.0, 9.0, 9.0, 9.0]
        @test MuRAT.fill_nearest!([NaN, NaN, 2.0, NaN]) == [2.0, 2.0, 2.0, 2.0]
        @test all(isnan, MuRAT.fill_nearest!([NaN, NaN]))
        A = reshape(collect(1.0:27.0), 3, 3, 3)
        B = copy(A); B[:, :, 1] .= NaN
        @test MuRAT.fill_nan3(B)[:, :, 1] == A[:, :, 2]
    end

    @testset "RMS envelope (MATLAB envelope(x, n, 'rms'))" begin
        x = [1.0, -2.0, 3.0, -4.0, 5.0, -6.0]
        xm = sum(x) / 6
        c = x .- xm
        e = MuRAT.rms_envelope(x, 4)          # even window: 1 sample before, 2 after
        @test e[1] ≈ xm + sqrt(sum(abs2, c[1:3]) / 3)
        @test e[3] ≈ xm + sqrt(sum(abs2, c[2:5]) / 4)
        @test e[6] ≈ xm + sqrt(sum(abs2, c[5:6]) / 2)
        @test MuRAT.rms_envelope(x, 3)[2] ≈ xm + sqrt(sum(abs2, c[1:3]) / 3)
    end

    @testset "pchip (MATLAB interp1 'pchip', 'extrap')" begin
        # values from MATLAB R2025b: interp1([0 1 2.5 4 7],[1 3 2 2.5 9],xq,'pchip','extrap')
        x, y = [0, 1, 2.5, 4, 7], [1.0, 3, 2, 2.5, 9]
        @test MuRAT.pchip(x, y, x) ≈ y
        p = MuRAT.pchip(x, y, collect(0:0.01:1))
        @test all(diff(p) .>= -1e-12)                 # monotone where the data are
        @test MuRAT.pchip([0.0, 1.0], [0.0, 2.0], [0.5, 2.0]) ≈ [1.0, 4.0]
    end

    @testset "checkerboard and robust fit" begin
        I = MuRAT.checkerboard3d((4, 4, 4), 2)
        @test I[1, 1, 1] == I[2, 2, 1]
        @test I[1, 1, 1] != I[3, 1, 1]
        @test I[1, 1, 1] == I[1, 1, 2] && I[1, 1, 1] != I[1, 1, 3]
        x = collect(1.0:20.0)
        y = 2 .* x .+ 1
        y[5] += 50
        p = MuRAT.robust_line_fit(x, y)
        @test isapprox(p[1], 2, atol = 1e-3) && isapprox(p[2], 1, atol = 1e-2)
    end

    @testset "Qc from a synthetic coda" begin
        srate, f, Qc = 100.0, 6.0, 0.004
        t = collect(0:1/srate:40)
        E = t .^ -1.5 .* exp.(-2π * f * Qc .* t)
        E[1] = E[2]
        sp = reshape(E, :, 1)
        cs = round(Int, 15 * srate)
        ce = cs + round(Int, 10 * srate) - 1
        q, unc = MuRAT.murat_qc([f], sp, cs, ce, t[cs], srate, "NonLinear")
        @test isapprox(q[1], Qc, rtol = 0.02)
        ql, _ = MuRAT.murat_qc([f], sp, cs, ce, t[cs], srate, "Linearized")
        @test isapprox(ql[1], Qc, rtol = 0.02)
    end

    @testset "Paasschens function" begin
        t0, A, N, coda = MuRAT.paasschens(5.0, 2.5, [0.8], [0.1], 0.5, 20.0)
        @test t0 == 2.0
        @test N == length(MuRAT.matlab_colon(2.0, 0.5, 21.0))
        @test all(coda[2:end, 1] .> 0)
        @test A[1] ≈ exp(-0.1 * 2.5 * 2.0) / (4π * 25 * 2.5)
    end
end
