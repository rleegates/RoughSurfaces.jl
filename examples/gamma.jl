using RoughSurfaces, Distributions, GLMakie, Printf

# Fixed random seed
randseed = 0

# Distribution Shape parameter
k = Observable(6.0)
# Standard deviation
σ = Observable(1.0)
# Hurst exponent
H = Observable(0.75)
# Power-law cut-off factor
plaw_fac = Observable(0.5)



# Domain
ftdom = FTDomain((512,512), (400.0, 400.0))

function parametric_surface(σ, H, k, plaw_fac, ftdom, randseed)
	θ = 1.0
	λr = 50.0
	λS = 0.01*λr
	λL = plaw_fac*λr
	distribution = Pareto(k)
	ftwhit = FTWhiteNoise(ftdom, distribution, randseed)
	psd = SelfAffinePSD(σ, H, λr, λL, λS)
	rough = RoughSurface(ftdom, psd, ftwhit)
	xy, z = RoughSurfaces.surface(rough)
	return (xy, z, rough.ftwhit, rough.psdkern)
end

function bounce(obs, pm, mn, repeat = 1)
	N = 100
	dx = 2π/N
	@async for _ = 1:repeat
		start = obs[]
		println()
		println("start: ", start)
		for i = 1:N
			bounce_i(obs, pm, mn, start, i, N)
			println(i,": ", obs[])
			sleep(0.1)
			yield()
		end
		println()
	end
end

function bounce_i(obs, pm, mn, start, i, N)
	dx = 2π/N
	obs[] = start + pm * sin(i*dx)
end

surfce(σ, H, k, plaw_fac) = parametric_surface(σ, H, k, plaw_fac, ftdom, randseed)
surface_observables = map(surfce, σ, H, k, plaw_fac)


x_obs = lift(surface_observables) do obs
	obs[1][1][(end-256):end]
end

y_obs = lift(surface_observables) do obs
	obs[1][2][(end-256):end]
end

z_obs = lift(surface_observables) do obs
	-obs[2][(end-256):end, (end-256):end]
end

latent_obs_vec = lift(surface_observables) do obs
	vec(-obs[3].sample)
end

z_obs_vec = lift(surface_observables) do obs
	vec(-obs[2])
end

psd_obs_vec = lift(surface_observables) do obs
	res = obs[4].S[2:256,1]
	res[res.<=0.0] .= NaN
	res
end

ttl = map((k, σ, H, plaw_fac)->@sprintf("Latent PDF parameter: %1.2f, Standard deviation: %1.2f mm, Hurst: %1.2f, Power-law cutoff wavelength: %1.2f mm", k, σ, H, plaw_fac*50.0), k, σ, H, plaw_fac)

#fig = Figure(fonts=(; regular = "CMU Serif Roman", bold = "Tenorite Regular"), size = (1400, 700))
fig = Figure(size = (1400, 700))

ax1 = Axis3(
	fig[1:11,4:14],
	title = ttl,
	#titlesize = 20,
	xlabel = L"$x$ (mm)",
	ylabel = L"$y$ (mm)",
	zlabel = L"$z$ (mm)",
	aspect = :data,
	azimuth = .35 * pi,
	elevation = 0.18 * pi,
	#protrusions = 50,
	perspectiveness = 0.0,
	zticks = WilkinsonTicks(3; k_min = 2)
)
GLMakie.surface!(x_obs, y_obs, z_obs, colormap = :terrain)
#zlims!(ax1, -10, 5)


ax2 = Axis(fig[1:4,1:2], title="Empirical surface PDF", xlabel = L"$h$ (mm)")
ylims!(ax2, 0.0, 1.0)
GLMakie.stephist!(z_obs_vec, bins=range(-3,3,step=0.01), normalization = :pdf, color = ("#FFC000", 1.0))

ax3 = Axis(fig[5:8,1:2], title="Empirical latent PDF", xlabel = L"$h$ (mm)")
ylims!(ax3, 0.0, 2.5)
GLMakie.stephist!(latent_obs_vec, bins=range(-3,3,step=0.01), normalization = :pdf, color = ("#FFC000", 1.0))

ax4 = Axis(fig[9:11,1:2], title="PSD", xscale=log10, yscale=log10, xlabel = L"$f$ (1/mm)")
GLMakie.lines!(ftdom.fxy[1][2:256], psd_obs_vec, color = ("#FFC000", 1.0))
ylims!(ax4, 1e-4, 1e3)

display(fig)



# wait(bounce(k, 3.0, 0.0))
# wait(bounce(σ, 0.5, 0.0))
# wait(bounce(H, 0.24, 0.0))
# wait(bounce(plaw_fac, 0.3, 0.0))

framerate = 30
duration = 10.0

start = k[]
rng = range(0.0, duration, step=1/framerate)
record(fig, "roughsurface_pdf_param.mp4", rng, framerate = framerate) do t
    bounce_i(k, 3.0, 0.0, start, t, duration)
end

start = σ[]
record(fig, "roughsurface_stddev.mp4", rng, framerate = framerate) do t
    bounce_i(σ, 0.5, 0.0, start, t, duration)
end

start = H[]
record(fig, "roughsurface_hurst.mp4", rng, framerate = framerate) do t
    bounce_i(H, 0.24, 0.0, start, t, duration)
end

start = plaw_fac[]
record(fig, "roughsurface_cutoff.mp4", rng, framerate = framerate) do t
    bounce_i(plaw_fac, 0.3, 0.0, start, t, duration)
end