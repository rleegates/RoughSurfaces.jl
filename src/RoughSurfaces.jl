module RoughSurfaces

include("psds.jl")

using FFTW, Distributions, StaticArrays
import RandomNumbers.Random: seed!
import Distributions: Sampleable
import Statistics: mean, std

function fftfreq(N, dx)
	L = N*dx
	if iseven(N)
		return vcat(0:1:(N/2-1), -(N/2):1:-1) / L
	else
		return vcat(0:1:((N-1)/2), -((N-1)/2):1:-1) / L
	end
end

struct FTDomain
	N::NTuple{2,Int}
	L::NTuple{2,Float64}
	xy::NTuple{2,Vector{Float64}}
	fxy::NTuple{2,Vector{Float64}}
	function FTDomain(N::NTuple{2,Int}, L::NTuple{2,Float64})
		(Nx, Ny) = N
		@assert iseven(Nx) && iseven(Ny)
		(Lx, Ly) = L
		dx = Lx / Nx
		dy = Ly / Ny
		x = ifftshift(-Lx/2 : dx : (Lx/2 - dx))
		y = ifftshift(-Ly/2 : dy : (Ly/2 - dy))
		fx = fftfreq(Nx, dx)
		fy = fftfreq(Ny, dy)
		return new(N, L, (x,y), (fx, fy))
	end
end

cu_fft(ftdom::FTDomain, x) = fft(ifftshift(x)) * (ftdom.L[1] * ftdom.L[2]) / (ftdom.N[1] * ftdom.N[2])
cu_ifft(ftdom::FTDomain, x) = fftshift(ifft(x)) * (ftdom.N[1] * ftdom.N[2]) / (ftdom.L[1] * ftdom.L[2])

export FTDomain
struct FTWhiteNoise
    ftdom::FTDomain
    sample::Matrix{Float64}
	Z::Matrix{Complex{Float64}}
	function FTWhiteNoise(ftdom::FTDomain, distribution::Sampleable, randseed::Int)
		seed!(randseed)
		(Nx, Ny) = ftdom.N
		smple = rand(distribution, Nx, Ny)
		empμ = mean(smple)
		meancentered = smple .- empμ
        empσ = std(meancentered)
		meancenterednormalized = meancentered/empσ
		Z = cu_fft(ftdom, meancenterednormalized)
		seed!()
		return new(ftdom, meancenterednormalized, Z)
	end
end

export FTWhiteNoise

function std(ftwhit::FTWhiteNoise)
    return std(ftwhit.sample)
end

function mean(ftwhit::FTWhiteNoise)
    return mean(ftwhit.sample)
end

mutable struct PSDKernel
    ftdom::FTDomain
	S::Matrix{Float64}
	G::Matrix{Float64}
	function PSDKernel(ftdom::FTDomain, psd::AbstractPSD)
		(Nx, Ny) = ftdom.N
		(Lx, Ly) = ftdom.L
		(fx, fy) = ftdom.fxy
		S = zeros(Nx, Ny)
		for j = 1:Ny, i = 1:Nx
			fxi = fx[i]
			fyj = fy[j]
			f = sqrt(fxi*fxi + fyj*fyj)
			S[i,j] = evaluate(psd, f)
		end
		dx = Lx / Nx
		dy = Ly / Ny
		G = 1/sqrt(dx) * 1/sqrt(dy) * sqrt.(S)
		return new(ftdom, S, G)
	end
end



struct RoughSurface
    ftdom::FTDomain
	ftrep::Matrix{Complex{Float64}}
	psdkern::PSDKernel
	ftwhit::FTWhiteNoise
	function RoughSurface(ftdom::FTDomain, psd::AbstractPSD, ftwhit::FTWhiteNoise)
		(Lx, Ly) = ftdom.L
        psdkern = PSDKernel(ftdom, psd)
		ftrep = psdkern.G .* ftwhit.Z * 1/Lx * 1/Ly
		return new(ftdom, ftrep, psdkern, ftwhit)
	end
end

export RoughSurface

function surface(rs::RoughSurface)
    (xs, ys) = rs.ftdom.xy
    (Lx, Ly) = rs.ftdom.L
	return ((fftshift(xs), fftshift(ys)), real(cu_ifft(rs.ftdom, rs.ftrep*Lx*Ly)))
end

export surface

function eval_surface(rs::RoughSurface, x::SVector{2, Float64})
	Nx, Ny = rs.ftdom.N
	z = zero(T)
	dz_dx = zero(T)
	dz_dy = zero(T)
	d2z_dx2 = zero(T)
	d2z_dy2 = zero(T)
	d2z_dxdy = zero(T)
	evx = x[1]
	evy = x[2]
	fx, fy = rs.ftdom.fxy
	f1 = im * 2 * π
	f2x = f1 * evx
	f2y = f1 * evy
	expy = zero(f2x)
	checkbounds(fx, Nx)
	checkbounds(fy, Ny)
	@inbounds for j = 0:(Ny-1)
		dfy = f1 * fy[j+1]
		expy = exp(f2y * fy[j+1])
		@simd for i = 0:(Nx-1)
			dfx = f1 * fx[i+1]
			ztmp = rs.ftrep[i+1, j+1] * exp(f2x * fx[i+1]) * expy
			z += ztmp
			dz_dx_tmp = dfx * ztmp
			dz_dy_tmp = dfy * ztmp
			dz_dx += dz_dx_tmp
			dz_dy += dz_dy_tmp
			d2z_dx2 += dfx * dz_dx_tmp
			d2z_dy2 += dfy * dz_dy_tmp
			d2z_dxdy += dfy * dz_dx_tmp
		 end
	end
	return (real(z), real(dz_dx), real(dz_dy), real(d2z_dx2), real(d2z_dy2), real(d2z_dxdy))
end

export eval_surface

end
