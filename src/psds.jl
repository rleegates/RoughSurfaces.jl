abstract type AbstractPSD end

struct SelfAffinePSD <: AbstractPSD
    σ::Float64
    H::Float64
    λr::Float64
    λL::Float64
    λS::Float64
end

export SelfAffinePSD

function evaluate(psd::SelfAffinePSD, f_eval::Float64)
	# Definition of Isotropic PSD of 2D surface
	# T. Jacobs, T. Junge, L. Pastewka.
	# Quantitative characterization of surface topography using spectral analysis
	# Eqs. 7 and 47
	# https://arxiv.org/pdf/1607.03040.pdf
	# WARNING: This is an isotropic PSD of a 2D surface, the other functions here
	# are 1D PSDs. These are not the same... see the paper above!!!

	σ, H, λr, λL, λS = psd.σ, psd.H, psd.λr, psd.λL, psd.λS
	@assert 0.5 <= H <= 1.0
	fr = 1/λr # largest wavelength, above no power
	fS = 1/λS # smallest wavelength, below no power
	fL = 1/λL # largest power law wavelength, above constant power, below powerlaw
	@assert fS >= fL >= fr
	qr = 2π*fr # roll-off wavenumber
	qL = 2π*fL # lower-cutoff wavenumber
	qS = 2π*fS # upper-cutoff wavenumber
	q = 2π*f_eval
	# RMS height is equal to standard deviation for surface of mean height zero
	rms_height = σ
	@assert qS > 5.0*qL
	# Eq. 47, assuming qS >> qL >> qr
	if qL == qr
		α = 1.0
	else
		α = 1/(1+H)
	end
	C0 = 4π*α*H*(rms_height/qL)^2 * (1/qL)^(-2-2H)
	if q < qr
		return 0.
	elseif qr ≤ q < qL
		return C0 * qL^(-2-2H)
	elseif qL ≤ q < qS
		return C0 * q^(-2-2H)
	else # q ≥ qS
		return 0.
	end
end