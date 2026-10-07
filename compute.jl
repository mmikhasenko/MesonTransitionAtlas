# Precompute every transition QuarkModelTransitions can evaluate, for ten meson
# sectors, and write data/transitions.json for the site renderer.
#
#   julia --project=. compute.jl [sector ids...]      (from this repository; see run.sh)
#
# Every operator call is attempted; a failure is recorded as a gap, never patched.

using GIModel, GIModel.QuarkModelTransitions
using Printf, Dates, Serialization, LinearAlgebra, TOML

const OUT = joinpath(@__DIR__, "data", "transitions.json")
params, mq = load_parameters_and_quark_masses(default_parameters_path())
const SOLVER = OscillatorSolver()

# All levels with n ≤ 3 in S, n ≤ 2 in P and D. The S+P subset is solved
# separately: in it nothing mixes with a D wave, which PhotonEmission and the
# gluonic kernel require. A process is tried on the full states first and on
# the S+P states only when the full ones are rejected.
const LEVELS_FULL = vcat(spectrum_levels(3; L_labels = ("S",)), spectrum_levels(2; L_labels = ("P", "D")))
const LEVELS_SP = filter(b -> b.L_label != "D", LEVELS_FULL)

# Strong couplings: the calibration of the strong-decay tutorial
# (h/g from the paper's S0/A fit, g from Γ(ρ → ππ) = 149.1 MeV).
const H_OVER_G = (3.291 / 1.644) / (3 - (3.291 / 1.644) / 4)
const ALPHA_EM_CHARGES = [AnnihilationTerm((:u, :u), 4 / 9), AnnihilationTerm((:d, :d), 1 / 9),
    AnnihilationTerm((:s, :s), 1 / 9), AnnihilationTerm((:c, :c), 4 / 9), AnnihilationTerm((:b, :b), 1 / 9)]

# ---------------------------------------------------------------------------
# States
# ---------------------------------------------------------------------------

"Replace averaged :q components by explicit (u ū + d d̄)/√2, as the currents need."
function expand_q(ps::PhysicalState; label = ps.label, mass_GeV = ps.mass_GeV)
    comps = Any[]
    for c in ps.components
        b = c.basis
        if b.flavors == (:q, :q)
            for f in (:u, :d)
                push!(comps, (basis = BasisState(b.n, b.L_label, b.multiplicity, b.J; flavors = (f, f)),
                    coefficient = c.coefficient / sqrt(2), wave = c.wave))
            end
        else
            push!(comps, (basis = c.basis, coefficient = c.coefficient, wave = c.wave))
        end
    end
    return PhysicalState(label, mass_GeV, comps; provenance = (source = :expanded_q,))
end

struct Rep
    sector::String
    real::String          # realization id within the sector ("+", "0", ...)
    key::String           # level key within the sector
    full::Union{Nothing,PhysicalState}
    sp::Union{Nothing,PhysicalState}
end

flavor_pairs(ps::PhysicalState) = Set(c.basis.flavors for c in ps.components)
anystate(r::Rep) = something(r.full, r.sp)

basis_label(b) = "$(b.n)^$(b.multiplicity)$(b.L_label)_$(b.J)"
key_of(s, iso::Bool) = iso ? basis_label(s.basis) * ":" * string(s.basis.flavors[1]) : s.label

function solve(f1, f2, levels)
    compute_spectrum(params, Meson(mq, f1, f2); levels, solver = SOLVER)
end

function solve_iso(levels)
    compute_isoscalar_spectrum(params, Meson(mq, :q, :q), Meson(mq, :s, :s); levels,
        pseudoscalar = PaperP2Annihilation(),
        amplitudes = Dict(("S", 3, 1) => params.annihilation.s1_A, ("P", 3, 2) => params.annihilation.a_3p2),
        solver = SOLVER)
end

"Dict key => PhysicalState for a spectrum."
function states_of(spec; iso = false, expand = false)
    d = Dict{String,PhysicalState}()
    for s in spec.states
        ps = physical_state(spec, s)
        d[key_of(s, iso)] = expand ? expand_q(ps) : ps
    end
    return d
end

# ---------------------------------------------------------------------------
# Sectors
# ---------------------------------------------------------------------------

# realizations: id => (label, flavor pair or :neutral/:iso); the first is canonical
const SECTORS = [
    (id = "ud", title = "Light isovector", quarks = "u d̄", hidden = true, symbol = :isovector,
        reals = [("+", "charged, u d̄", (:u, :d)), ("0", "neutral, (uū − dd̄)/√2", :neutral), ("-", "charged, d ū", (:d, :u))]),
    (id = "iso", title = "Light isoscalar", quarks = "nn̄ + ss̄", hidden = true, symbol = :isoscalar,
        reals = [("0", "nn̄ and ss̄ with annihilation mixing", :iso)]),
    (id = "K", title = "Strange", quarks = "u s̄", hidden = false, symbol = ("K", ""),
        reals = [("+", "K⁺-like, u s̄", (:u, :s)), ("0", "K⁰-like, d s̄", (:d, :s)), ("-", "K⁻-like, s ū", (:s, :u)), ("0bar", "K̄⁰-like, s d̄", (:s, :d))]),
    (id = "D", title = "Charm", quarks = "c q̄", hidden = false, symbol = ("D", ""),
        reals = [("0", "D⁰-like, c ū", (:c, :u)), ("+", "D⁺-like, c d̄", (:c, :d))]),
    (id = "Ds", title = "Charm-strange", quarks = "c s̄", hidden = false, symbol = ("D", "s"),
        reals = [("+", "D_s⁺, c s̄", (:c, :s))]),
    (id = "cc", title = "Charmonium", quarks = "c c̄", hidden = true, symbol = :cc,
        reals = [("0", "c c̄", (:c, :c))]),
    (id = "B", title = "Bottom", quarks = "q b̄", hidden = false, symbol = ("B", ""),
        reals = [("+", "B⁺-like, u b̄", (:u, :b)), ("0", "B⁰-like, d b̄", (:d, :b))]),
    (id = "Bs", title = "Bottom-strange", quarks = "s b̄", hidden = false, symbol = ("B", "s"),
        reals = [("0", "B_s⁰, s b̄", (:s, :b))]),
    (id = "Bc", title = "Bottom-charm", quarks = "c b̄", hidden = false, symbol = ("B", "c"),
        reals = [("+", "B_c⁺, c b̄", (:c, :b))]),
    (id = "bb", title = "Bottomonium", quarks = "b b̄", hidden = true, symbol = :bb,
        reals = [("0", "b b̄", (:b, :b))]),
]

const ONLY = isempty(ARGS) ? nothing : Set(ARGS)

t0 = time()
reps = Dict{Tuple{String,String},Dict{String,Rep}}()   # (sector, real) => key => Rep
mkpath(joinpath(@__DIR__, "cache"))
const CACHE = joinpath(@__DIR__, "cache", "states.jls")   # delete to re-solve
if isfile(CACHE)
    reps = deserialize(CACHE)
    println("loaded solved states from $CACHE"); flush(stdout)
else
    for sec in SECTORS
        for (rid, _, fl) in sec.reals
            if fl == :iso
                full, sp = solve_iso(LEVELS_FULL), solve_iso(LEVELS_SP)
                F, S = states_of(full; iso = true, expand = true), states_of(sp; iso = true, expand = true)
            elseif fl == :neutral
                uuF, ddF = solve(:u, :u, LEVELS_FULL), solve(:d, :d, LEVELS_FULL)
                uuS, ddS = solve(:u, :u, LEVELS_SP), solve(:d, :d, LEVELS_SP)
                mk(a, b) = Dict(k => superpose([v, b[k]], [1, -1]; label = k, mass_GeV = v.mass_GeV) for (k, v) in a)
                F, S = mk(states_of(uuF), states_of(ddF)), mk(states_of(uuS), states_of(ddS))
            else
                full, sp = solve(fl..., LEVELS_FULL), solve(fl..., LEVELS_SP)
                F, S = states_of(full), states_of(sp)
            end
            reps[(sec.id, rid)] = Dict(k => Rep(sec.id, rid, k, v, get(S, k, nothing)) for (k, v) in F)
            @printf("solved %-4s %-5s  %3d states  (%.0f s)\n", sec.id, rid, length(F), time() - t0); flush(stdout)
        end
    end
    serialize(CACHE, reps)
end
haskey(ENV, "SOLVE_ONLY") && exit()

# isovector neutral and isoscalar levels share masses with the canonical ones
canonical(secid) = reps[(secid, first(first(s for s in SECTORS if s.id == secid).reals)[1])]

# ---------------------------------------------------------------------------
# Names
# ---------------------------------------------------------------------------

# Names use "_{...}" for subscripts; the renderer typesets them.
const HIDDEN_NAMES = Dict(
    :isovector => Dict("1S0" => "π", "3S1" => "ρ", "1P1" => "b_{1}", "3P0" => "a_{0}", "3P1" => "a_{1}", "3P2" => "a_{2}",
        "1D2" => "π_{2}", "3D1" => "ρ", "3D2" => "ρ_{2}", "3D3" => "ρ_{3}"),
    :iso_n => Dict("1S0" => "η", "3S1" => "ω", "1P1" => "h_{1}", "3P0" => "f_{0}", "3P1" => "f_{1}", "3P2" => "f_{2}",
        "1D2" => "η_{2}", "3D1" => "ω", "3D2" => "ω_{2}", "3D3" => "ω_{3}"),
    :iso_s => Dict("1S0" => "η′", "3S1" => "φ", "1P1" => "h_{1}′", "3P0" => "f_{0}′", "3P1" => "f_{1}′", "3P2" => "f_{2}′",
        "1D2" => "η_{2}′", "3D1" => "φ", "3D2" => "φ_{2}", "3D3" => "φ_{3}"),
    :cc => Dict("1S0" => "η_{c}", "3S1" => "ψ", "1P1" => "h_{c}", "3P0" => "χ_{c0}", "3P1" => "χ_{c1}", "3P2" => "χ_{c2}",
        "1D2" => "η_{c2}", "3D1" => "ψ", "3D2" => "ψ_{2}", "3D3" => "ψ_{3}"),
    :bb => Dict("1S0" => "η_{b}", "3S1" => "Υ", "1P1" => "h_{b}", "3P0" => "χ_{b0}", "3P1" => "χ_{b1}", "3P2" => "χ_{b2}",
        "1D2" => "η_{b2}", "3D1" => "Υ", "3D2" => "Υ_{2}", "3D3" => "Υ_{3}"))
# open flavor: (J subscript, star) by spectroscopic code
const OPEN_DECO = Dict("1S0" => ("", ""), "3S1" => ("", "*"), "1P1" => ("1", ""), "3P0" => ("0", "*"),
    "3P1" => ("1", ""), "3P2" => ("2", "*"), "1D2" => ("2", ""), "3D1" => ("1", "*"), "3D2" => ("2", ""), "3D3" => ("3", "*"))

function level_name(sec, key)
    b = parse_label(key)
    nL = "($(b.n)$(b.L))"
    code = "$(b.S)$(b.L)$(b.J)"
    sym = sec.symbol
    sym == :isoscalar && (sym = b.flavor == "s" ? :iso_s : :iso_n)
    if sym isa Symbol
        sym == :cc && code == "3S1" && b.n == 1 && return "J/ψ"
        return HIDDEN_NAMES[sym][code] * nL
    end
    base, flav = sym               # e.g. ("D", "s")
    jsub, star = OPEN_DECO[code]
    sub = flav * jsub
    return base * (isempty(sub) ? "" : "_{" * sub * "}") * star * nL
end

function parse_label(key)
    m = match(r"^(\d)\^(\d)([SPD])_(\d)(?::(\w))?$", key)
    m === nothing && error("bad key $key")
    return (n = parse(Int, m[1]), S = parse(Int, m[2]), L = m[3], J = parse(Int, m[4]), flavor = something(m[5], ""))
end

# ---------------------------------------------------------------------------
# JSON (no extra dependency)
# ---------------------------------------------------------------------------

js(x::AbstractString) = "\"" * replace(x, "\\" => "\\\\", "\"" => "\\\"", "\n" => "\\n") * "\""
js(x::Symbol) = js(String(x))
js(x::Bool) = x ? "true" : "false"
js(::Nothing) = "null"
js(x::Integer) = string(x)
js(x::Real) = isfinite(x) ? @sprintf("%.6g", x) : "null"
js(x::AbstractVector) = "[" * join(map(js, x), ",") * "]"
js(x::AbstractDict) = "{" * join([js(string(k)) * ":" * js(v) for (k, v) in x], ",") * "}"
js(x::NamedTuple) = "{" * join([js(string(k)) * ":" * js(v) for (k, v) in pairs(x)], ",") * "}"

# ---------------------------------------------------------------------------
# Processes
# ---------------------------------------------------------------------------

gaps = Dict{String,Dict{String,Int}}()
const GAPLOCK = ReentrantLock()
gap!(secid, what, e) = lock(() -> _gap!(secid, what, e), GAPLOCK)
function _gap!(secid, what, e)
    msg = first(split(sprint(showerror, e), '\n'))
    msg = replace(msg, r"\d+\^\d[SPD]_\d" => "·")       # merge messages that differ only by level
    msg = length(msg) > 140 ? msg[1:prevind(msg, 141)] * "…" : msg
    d = get!(gaps, secid, Dict{String,Int}())
    d[what * ": " * msg] = get(d, what * ": " * msg, 0) + 1
end

"Try f on full states, then on S+P states. Returns (value, basis) or nothing."
function attempt(f, secid, what, reps_...; sink = gap!)
    try
        any(r -> r.full === nothing, reps_) && error("no full state")
        return (f(map(r -> r.full, reps_)...), "full")
    catch e1
        if all(r -> r.sp !== nothing, reps_)
            try
                return (f(map(r -> r.sp, reps_)...), "SP")
            catch e2
                sink(secid, what, e2)
                return nothing
            end
        end
        sink(secid, what, e1)
        return nothing
    end
end

records = Dict{String,Vector{Any}}(s.id => Any[] for s in SECTORS)
emit!(secid; kw...) = push!(records[secid], (; kw...))

active(sec) = ONLY === nothing || sec.id in ONLY
const WORKER = haskey(ENV, "STRONG_WORKER")   # strong-pass worker: skip everything else

# --- radiative, within each realization -------------------------------------
pe = PhotonEmission(mq)
for sec in SECTORS
    (active(sec) && !WORKER) || continue
    for (rid, _, _) in sec.reals
        R = reps[(sec.id, rid)]
        for (ki, ri) in R, (kf, rf) in R
            anystate(ri).mass_GeV > anystate(rf).mass_GeV || continue
            res = attempt(sec.id, "photon", ri, rf) do i, f
                me = matrix_element(f, pe, i)
                (String(me.multipole), decay_width(me))
            end
            res === nothing && continue
            (mult, w), basis = res
            st = anystate(ri)
            emit!(sec.id; cls = mult, from = ki, to = kf, width = w, real = length(sec.reals) > 1 ? rid : nothing, basis)
        end
    end
    @printf("radiative %-4s %5d records (%.0f s)\n", sec.id, length(records[sec.id]), time() - t0); flush(stdout)
end

# --- radiative across the neutral light sectors (ω → π⁰γ, ρ⁰ → ηγ, ...) ----
for (a, b) in ((("iso", "0"), ("ud", "0")), (("ud", "0"), ("iso", "0")))
    (ONLY === nothing || a[1] in ONLY) && !WORKER || continue
    for (ki, ri) in reps[a], (kf, rf) in reps[b]
        anystate(ri).mass_GeV > anystate(rf).mass_GeV || continue
        res = attempt(a[1], "photon (cross)", ri, rf) do i, f
            me = matrix_element(f, pe, i)
            (String(me.multipole), decay_width(me))
        end
        res === nothing && continue
        (mult, w), basis = res
        emit!(a[1]; cls = mult, from = ki, to = nothing, toSector = b[1], toKey = kf, emitted = "γ", width = w,
            real = a[1] == "ud" ? "0" : nothing, basis, flag = nothing)
    end
end

# --- annihilation ------------------------------------------------------------
ee = LeptonicCurrent(:electromagnetic, mq)
aa = TwoPhotonAnnihilation(mq, ALPHA_EM_CHARGES)
glu = GluonicAnnihilation(mq, M -> alpha_s_q(M))
leptonic = Dict(   # sector => (realization, CKM, [(lepton, mass)])
    "ud" => ("+", 0.97367, [("μν", 0.105658)]), "K" => ("+", 0.2243, [("μν", 0.105658)]),
    "D" => ("+", 0.221, [("μν", 0.105658)]), "Ds" => ("+", 0.975, [("μν", 0.105658), ("τν", 1.77686)]),
    "B" => ("+", 0.00382, [("μν", 0.105658), ("τν", 1.77686)]), "Bc" => ("+", 0.0408, [("τν", 1.77686)]))
for sec in SECTORS
    (active(sec) && !WORKER) || continue
    neutral = sec.hidden ? (sec.id == "ud" ? "0" : first(sec.reals)[1]) : nothing
    if neutral !== nothing
        for (k, r) in reps[(sec.id, neutral)]
            st = anystate(r)
            if st.J == 1 && st.parity == -1
                res = attempt(sec.id, "e⁺e⁻", r) do i
                    decay_width(MasslessLeptonPair(), ee, i)
            end
                res === nothing || emit!(sec.id; cls = "ee", from = k, to = nothing, emitted = "e⁺e⁻",
                    width = res[1], basis = res[2], real = sec.id == "ud" ? "0" : nothing)
            end
            if (st.J == 0 && st.parity == -1) || (st.J == 2 && st.parity == 1)
                res = attempt(sec.id, "γγ", r) do i
                    decay_width(TwoPhotonChannel(), aa, i)
            end
                res === nothing || emit!(sec.id; cls = "γγ", from = k, to = nothing, emitted = "γγ",
                    width = res[1], basis = res[2], real = sec.id == "ud" ? "0" : nothing)
            end
            if sec.id in ("cc", "bb")
                for (cls, ch) in (("gg", TwoGluonChannel()), ("ggg", ThreeGluonChannel()))
                    res = attempt(sec.id, cls, r) do i
                        decay_width(ch, glu, i)
                    end
                    res === nothing || emit!(sec.id; cls, from = k, to = nothing, emitted = cls,
                        width = res[1], basis = res[2], real = nothing)
            end
            end
        end
    end
    if haskey(leptonic, sec.id)
        rid, ckm, leptons = leptonic[sec.id]
        r = reps[(sec.id, rid)]["1^1S_0"]
        fl = anystate(r).components[1].basis.flavors
        op = LeptonicCurrent(:P_P, mq, AnnihilationTerm(fl, 2sqrt(3)))
        for (name, ml) in leptons
            res = attempt(sec.id, "ℓν", r) do i
                decay_width(LeptonNeutrinoChannel(ml; ckm), op, i)
            end
            res === nothing || emit!(sec.id; cls = "lnu", from = "1^1S_0", to = nothing, emitted = name,
                width = res[1], basis = res[2], real = length(sec.reals) > 1 ? rid : nothing, ckm)
        end
    end
end
println("annihilation done ($(round(Int, time() - t0)) s)"); flush(stdout)

# Strong only: mixed states carry radial admixtures with coefficients near zero,
# and each still costs a full overlap in every amplitude. Dropping |c| < PRUNE
# changes a width by about PRUNE relative. Not applied to the radiative and
# annihilation passes above, where it could hide a refused D-wave admixture.
#
# The waves are also sampled once onto a shared mesh. Every oscillator-wave
# overlap is a quadgk integral at rtol 1e-10 over the full Laguerre series (or a
# rebuilt operator matrix), and one amplitude needs dozens of them; on the mesh
# each is a plain sum. Checked against the oscillator overlaps for ρ → ππ,
# ρ(2D) → ωπ, a₂ → ρπ, K*(2S) → K*η, D* → Dπ: widths agree to 1e-6 relative
# at h = 0.005 GeV⁻¹, 28-67 times faster on warm calls (GIModel.jl#21).
const PRUNE = 1e-4
const MESH = collect(0.005:0.005:40.0)   # GeV⁻¹
prune(::Nothing) = nothing
prune(ps::PhysicalState) = PhysicalState(ps.label, ps.mass_GeV,
    [(basis = c.basis, coefficient = c.coefficient, wave = sample_wave(c.wave, MESH))
     for c in ps.components if abs(c.coefficient) >= PRUNE];
    provenance = ps.provenance)
for (k, R) in reps, (key, r) in R
    R[key] = Rep(r.sector, r.real, r.key, prune(r.full), prune(r.sp))
end

# --- strong: elementary emission of π, K, η, η′ -------------------------------
unit = attempt("ud", "strong calibration", canonical("ud")["1^3S_1"], canonical("ud")["1^1S_0"], reps[("ud", "0")]["1^1S_0"]) do p, s, e
    decay_width(TwoMesonChannel(s, e), PseudoscalarEmission(1.0, H_OVER_G, mq), p)
end
g = sqrt(149.1 / unit[1])
strong = PseudoscalarEmission(g, H_OVER_G * g, mq)
@printf("strong couplings: g = %.4f, h = %.4f GeV⁻¹\n", g, H_OVER_G * g)

# emitted species, in order of preference when both daughters are pseudoscalars
emitted = [
    ("π", 0, "π⁺", reps[("ud", "+")]["1^1S_0"]), ("π", 0, "π⁻", reps[("ud", "-")]["1^1S_0"]), ("π", 0, "π⁰", reps[("ud", "0")]["1^1S_0"]),
    ("K", 1, "K⁺", reps[("K", "+")]["1^1S_0"]), ("K", 1, "K⁰", reps[("K", "0")]["1^1S_0"]),
    ("K̄", 1, "K⁻", reps[("K", "-")]["1^1S_0"]), ("K̄", 1, "K̄⁰", reps[("K", "0bar")]["1^1S_0"]),
    ("η", 2, "η", reps[("iso", "0")]["1^1S_0:q"]), ("η′", 3, "η′", reps[("iso", "0")]["1^1S_0:s"]),
]
emitted_index = Dict((e[4].sector, e[4].real, e[4].key) => i for (i, e) in enumerate(emitted))
BLAS.set_num_threads(1)

function flavor_allowed(P, S, E)
    for (a, b) in flavor_pairs(P), (x, y) in flavor_pairs(S), (z, w) in flavor_pairs(E)
        (z == a && w == x && y == b) && return true      # quark a emits (a w̄), becomes w
        (w == b && z == y && x == a) && return true      # antiquark b̄ emits (z b̄), becomes z̄
    end
    return false
end

survivors = [r for ((sid, _), R) in reps if !(sid in ("cc", "bb")) for r in values(R)]

# One cache file per parent: the strong pass is long, and every amplitude
# allocates heavily, so threads end up waiting on the garbage collector.
# Separate worker processes (run_strong.sh) scale; each fills part of the cache.
const STRONG_CACHE = joinpath(@__DIR__, "cache", "strong")
mkpath(STRONG_CACHE)
strong_file(secid, key) = joinpath(STRONG_CACHE, secid * "__" * replace(key, "^" => "", ":" => "-") * ".jls")

"All decays of one parent by pseudoscalar emission, with the refused calls."
function strong_parent(secid, kp, rp)
    P = anystate(rp)
    acc = Dict{Tuple{String,String,String},Any}()   # (sector, key, species) => (width, channels, basis)
    local_gaps = Any[]
    sink = (s, what, e) -> push!(local_gaps, (s, what, e))
    tp, ncalls = time(), 0
    for rs in survivors, (species, rank, ename, re) in emitted
        S, E = anystate(rs), anystate(re)
        P.mass_GeV > S.mass_GeV + E.mass_GeV || continue
        flavor_allowed(P, S, E) || continue
        # Two pseudoscalar daughters: keep one Fig.-14 assignment, the one whose
        # emitted meson is lighter in the order π < K < η < η′ (tutorial convention).
        sidx = get(emitted_index, (rs.sector, rs.real, rs.key), 0)
        if sidx > 0
            srank = emitted[sidx][2]
            (rank < srank || (rank == srank && emitted_index[(re.sector, re.real, re.key)] >= sidx)) || continue
        end
        ncalls += 1
        res = attempt(secid, "strong", rp, rs, re; sink) do p, s, e
            decay_width(TwoMesonChannel(s, e), strong, p)
        end
        res === nothing && continue
        w, basis = res
        (rs.sector, rs.real, rs.key) == (re.sector, re.real, re.key) && (w /= 2)   # identical daughters
        w > 0 || continue
        sname = sidx > 0 ? emitted[sidx][3] : rs.key
        k = (rs.sector, rs.key, species)
        old = get(acc, k, (0.0, String[], "full"))
        acc[k] = (old[1] + w, push!(old[2], "$(rs.real)$(sname == rs.key ? "" : " " * sname) + $ename"),
            basis == "SP" ? "SP" : old[3])
    end
    @printf("  %-4s %-12s %4d calls %6.1f s\n", secid, kp, ncalls, time() - tp); flush(stdout)
    result = (acc = acc, gaps = [(s, what, first(split(sprint(showerror, e), '\n'))) for (s, what, e) in local_gaps])
    serialize(strong_file(secid, kp), result)
    return result
end

zeros_count = Dict{String,Dict{String,Int}}()

strong_parents = [(sec.id, kp, rp) for sec in SECTORS if active(sec) && !(sec.id in ("cc", "bb"))
                  for (kp, rp) in sort(collect(reps[(sec.id, first(sec.reals)[1])]); by = first)]

# Worker mode (STRONG_WORKER=i/N): fill the cache for every N-th parent, then stop.
if haskey(ENV, "STRONG_WORKER")
    i, N = parse.(Int, split(ENV["STRONG_WORKER"], "/"))
    for (j, (secid, kp, rp)) in enumerate(strong_parents)
        (j - 1) % N == i || continue
        isfile(strong_file(secid, kp)) && continue
        strong_parent(secid, kp, rp)
    end
    println("worker $i/$N done ($(round(Int, time() - t0)) s)")
    exit()
end

for sec in SECTORS
    active(sec) || continue
    sec.id in ("cc", "bb") && continue
    rid = first(sec.reals)[1]
    parents = sort(collect(reps[(sec.id, rid)]); by = first)
    results = Vector{Any}(undef, length(parents))
    Threads.@threads :dynamic for ip in eachindex(parents)
        kp, rp = parents[ip]
        f = strong_file(sec.id, kp)
        results[ip] = (kp, isfile(f) ? deserialize(f) : strong_parent(sec.id, kp, rp))
    end
    # Two kinds of exact zeros are counted, not drawn. A channel with no (L, S)
    # allowed by J and parity makes partial_width sum over nothing and throw
    # (GIModel.jl#23). Channels forbidden by G parity or isospin come out as
    # rounding noise: below 1e-29 MeV, with nothing between 1e-30 and 1e-9.
    nz = get!(zeros_count, sec.id, Dict("parity" => 0, "symmetry" => 0))
    for (kp, r) in results, (s, what, msg) in r.gaps
        if occursin("reducing over an empty collection", msg)
            nz["parity"] += 1
        else
            gap!(s, what, ErrorException(msg))
        end
    end
    for (kp, r) in results, ((ssec, skey, species), (w, chans, basis)) in r.acc
        w < 1e-20 && (nz["symmetry"] += 1; continue)
        emit!(sec.id; cls = "strong", from = kp, to = ssec == sec.id ? skey : nothing,
            toSector = ssec, toKey = skey, emitted = species, width = w, basis, channels = sort(chans), real = nothing)
    end
    @printf("strong %-4s done (%.0f s)\n", sec.id, time() - t0); flush(stdout)
end

# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

function composition(ps)
    parts = String[]
    for c in ps.components
        abs(c.coefficient) < 0.05 && continue
        b = c.basis
        f = b.flavors == (:u, :u) ? "uū" : b.flavors == (:d, :d) ? "dd̄" : b.flavors == (:s, :s) ? "ss̄" : ""
        push!(parts, @sprintf("%+.2f %s%s", c.coefficient, isempty(f) ? "" : f * " ", basis_label(b)))
    end
    return join(parts, ", ")
end

# The GIModel.jl version is the one pinned in Manifest.toml.
const GIMODEL = only(TOML.parsefile(joinpath(@__DIR__, "Manifest.toml"))["deps"]["GIModel"])

sector_json = String[]
for sec in SECTORS
    active(sec) || continue
    R = canonical(sec.id)
    levels = [(key = k, name = level_name(sec, k), mass = anystate(r).mass_GeV,
        n = parse_label(k).n, S = parse_label(k).S, L = parse_label(k).L, J = parse_label(k).J,
        slot = sec.id == "iso" ? (endswith(k, ":s") ? 1 : 0) : 0,
        composition = composition(anystate(r))) for (k, r) in sort(collect(R); by = first)]
    gl = [(message = m, count = c) for (m, c) in sort(collect(get(gaps, sec.id, Dict())); by = x -> -x[2])]
    push!(sector_json, js((id = sec.id, title = sec.title, quarks = sec.quarks, hidden = sec.hidden,
        realizations = [(id = r[1], label = r[2]) for r in sec.reals],
        levels, transitions = records[sec.id], gaps = gl,
        strong_zeros = get(zeros_count, sec.id, Dict("parity" => 0, "symmetry" => 0)))))
end

meta = (generated = string(now()), gimodel_version = GIMODEL["version"], gimodel_commit = GIMODEL["repo-rev"][1:7], gimodel_url = GIMODEL["repo-url"],
    solver = string(numerics_provenance(SOLVER)), g = g, h = H_OVER_G * g,
    levels = "S: n ≤ 3; P, D: n ≤ 2", prune = PRUNE, strong_mesh_h = MESH[2] - MESH[1], threads = Threads.nthreads(), seconds = round(Int, time() - t0))
open(OUT, "w") do io
    print(io, "{\"meta\":", js(meta), ",\"sectors\":[", join(sector_json, ","), "]}")
end
# The site loads data.js through a script tag, so it also opens from disk.
write(joinpath(@__DIR__, "site", "data.js"), "window.TRANSITION_DATA = " * read(OUT, String) * ";\n")
println("wrote $OUT in $(round(Int, time() - t0)) s")
