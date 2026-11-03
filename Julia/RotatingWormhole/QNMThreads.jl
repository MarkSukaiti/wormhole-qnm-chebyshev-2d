using DelimitedFiles
using LinearAlgebra
using GenericSchur
using Printf
using Base.Threads

cd(@__DIR__)

# ============================================================
# Numerical configuration
# ============================================================

const DECIMAL_DIGITS = 200
const BIGFLOAT_BITS = ceil(Int, (DECIMAL_DIGITS + 20) * log2(10))

setprecision(BigFloat, BIGFLOAT_BITS)

const BF = BigFloat
const CBF = Complex{BigFloat}

# Three radial resolutions used in the convergence test.
const RESOLUTIONS = (90, 95, 100)

# Number of angular collocation points.
const ANGULAR_RESOLUTION = 3

# Spin labels J.
#
# The input-file suffix (spin tag) is k = round(Int, 10J).  The physical spin encoded in the
# Maple matrices of tag k is j = J/r0^2 = j*/(1 - j*) with j* = k/10 and, for k = 10,
# j* = 0.95, i.e. j = 0, 1/9, 1/4, 3/7, 2/3, 1, 3/2, 7/3, 4, 9, 19.  The value J written in
# the output headers is the label k/10, not the physical spin.
#
# For example:
#     J = 0.0  -> suffix 0  -> j = 0
#     J = 0.1  -> suffix 1  -> j = 1/9
#     J = 0.5  -> suffix 5  -> j = 1
#     J = 1    -> suffix 10 -> j = 19
#
const SPIN_VALUES = [
    0,
    0.1,
    0.2,
    0.3,
    0.4,
    0.5,
    0.6,
    0.7,
    0.8,
    0.9,
    1
]

const EPS_ABS = BigFloat("1e-3")
const EPS_REL = BigFloat("1e-3")
const INSTABILITY_THRESHOLD = BigFloat("1e-18")

# Azimuthal sector m (the Maple generator writes one data folder per m).  Choose it here or
# on the command line:   QNM_M=2 julia -t 12 QNMThreads.jl
const M_SECTOR = parse(Int, get(ENV, "QNM_M", "1"))
const DATA_DIRECTORY = "data-m$(M_SECTOR)"
const OUTPUT_DIRECTORY = "output-m$(M_SECTOR)"

mkpath(OUTPUT_DIRECTORY)

# Prevent two different J values from using the same rounded suffix.
@assert length(unique(round.(Int, 10 .* SPIN_VALUES))) == length(SPIN_VALUES) "Two or more spin values produce the same round(10J) suffix."

# ============================================================
# Input-file handling
# ============================================================

spin_tag(J::Real) = round(Int, 10 * J)

function find_matrix_file(
    prefix::String,
    Nres::Int,
    Mres::Int,
    Jtag::Int;
    dir::String = DATA_DIRECTORY
)
    stem = "$(prefix)-$(Nres)-$(Mres)-$(Jtag)"

    candidates = (
        joinpath(dir, stem),
        joinpath(dir, stem * ".txt"),
        joinpath(dir, stem * ".mat")
    )

    for path in candidates
        if isfile(path)
            return path
        end
    end

    error(
        "Could not find matrix $(stem). Tried:\n" *
        join(["    " * path for path in candidates], "\n")
    )
end

function read_big_matrix(path::AbstractString)
    raw = readdlm(path, BigFloat)
    return Matrix{BigFloat}(raw)
end

function read_big_real_matrix(
    prefix::String,
    Nres::Int,
    Mres::Int,
    Jtag::Int;
    dir::String = DATA_DIRECTORY
)
    path = find_matrix_file(
        prefix,
        Nres,
        Mres,
        Jtag;
        dir = dir
    )

    @printf(
        "[thread %d] Reading %s\n",
        threadid(),
        path
    )

    return read_big_matrix(path)
end

function make_complex_matrix(
    real_part::AbstractMatrix{BigFloat},
    imaginary_part::AbstractMatrix{BigFloat}
)
    @assert size(real_part) == size(imaginary_part) "The real and imaginary matrix parts have different dimensions."

    matrix = Matrix{CBF}(undef, size(real_part)...)

    @inbounds for column in axes(real_part, 2)
        for row in axes(real_part, 1)
            matrix[row, column] = CBF(
                real_part[row, column],
                imaginary_part[row, column]
            )
        end
    end

    return matrix
end

function load_qep_matrices(
    Nres::Int,
    Mres::Int,
    Jtag::Int;
    dir::String = DATA_DIRECTORY
)
    M2 = read_big_real_matrix(
        "M2",
        Nres,
        Mres,
        Jtag;
        dir = dir
    )

    M1R = read_big_real_matrix(
        "M1R",
        Nres,
        Mres,
        Jtag;
        dir = dir
    )

    M1I = read_big_real_matrix(
        "M1I",
        Nres,
        Mres,
        Jtag;
        dir = dir
    )

    M0R = read_big_real_matrix(
        "M0R",
        Nres,
        Mres,
        Jtag;
        dir = dir
    )

    M0I = read_big_real_matrix(
        "M0I",
        Nres,
        Mres,
        Jtag;
        dir = dir
    )

    M2_complex = CBF.(M2)
    M1_complex = make_complex_matrix(M1R, M1I)
    M0_complex = make_complex_matrix(M0R, M0I)

    return M0_complex, M1_complex, M2_complex
end

# ============================================================
# Quadratic eigenvalue problem
# ============================================================

function qep_eigenvalues(
    M0::AbstractMatrix{Complex{T}},
    M1::AbstractMatrix{Complex{T}},
    M2::AbstractMatrix{Complex{T}}
) where {T<:Real}

    matrix_dimension = size(M0, 1)

    @assert size(M0, 2) == matrix_dimension "M0 must be square."
@assert size(M1) == (matrix_dimension, matrix_dimension) "M1 must have the same dimensions as M0."
@assert size(M2) == (matrix_dimension, matrix_dimension) "M2 must have the same dimensions as M0."

    complex_type = Complex{T}

    zero_matrix = zeros(
        complex_type,
        matrix_dimension,
        matrix_dimension
    )

    identity_matrix = Matrix{complex_type}(
        I,
        matrix_dimension,
        matrix_dimension
    )

    # The quadratic problem is
    #
    #     (M0 + omega M1 + omega^2 M2) a = 0.
    #
    # Its companion linearisation is
    #
    # [  0    I  ] [a        ]     [ I   0  ] [a        ]
    # [ -M0  -M1 ] [omega a  ] = ω [ 0   M2 ] [omega a  ].
    #
    # This assumes that M1R and M1I already contain the complete
    # coefficient multiplying omega.

    A = [
        zero_matrix       identity_matrix
        -M0               -M1
    ]

    B = [
        identity_matrix   zero_matrix
        zero_matrix       M2
    ]

    @printf(
        "[thread %d] Solving %d x %d generalised problem\n",
        threadid(),
        2 * matrix_dimension,
        2 * matrix_dimension
    )

    standard_matrix = B \ A
    eigenvalues = eigvals(standard_matrix)

    return CBF.(eigenvalues)
end

# ============================================================
# Output functions
# ============================================================

function save_complex_values(
    filename::String,
    values;
    header::String = "# Re(omega) Im(omega)"
)
    open(filename, "w") do io
        println(io, header)

        for omega in values
            println(
                io,
                real(omega),
                " ",
                imag(omega)
            )
        end
    end
end

function save_average_modes(filename::String, modes)
    open(filename, "w") do io
        println(io, "# Re(omega_average) Im(omega_average)")

        for mode in modes
            println(
                io,
                real(mode.omega_average),
                " ",
                imag(mode.omega_average)
            )
        end
    end
end

function save_mode_diagnostics(
    filename::String,
    modes;
    J::Real,
    Jtag::Int
)
    N1, N2, N3 = RESOLUTIONS

    open(filename, "w") do io
        println(io, "# Three-resolution convergence results")
        println(io, "# J = $(J)")
        println(io, "# Input-file spin suffix = $(Jtag)")
        println(io, "# N = $(N1), $(N2), $(N3)")
        println(io, "# M = $(ANGULAR_RESOLUTION)")
        println(io, "# eps_abs = $(EPS_ABS)")
        println(io, "# eps_rel = $(EPS_REL)")
        println(io, "#")
        println(
            io,
            "# Re(omega_av) Im(omega_av) ",
            "Re(omega_$(N1)) Im(omega_$(N1)) ",
            "Re(omega_$(N2)) Im(omega_$(N2)) ",
            "Re(omega_$(N3)) Im(omega_$(N3)) ",
            "d12 d23 d13 max_distance tolerance"
        )

        for mode in modes
            println(
                io,
                real(mode.omega_average), " ",
                imag(mode.omega_average), " ",
                real(mode.omega1), " ",
                imag(mode.omega1), " ",
                real(mode.omega2), " ",
                imag(mode.omega2), " ",
                real(mode.omega3), " ",
                imag(mode.omega3), " ",
                mode.distance12, " ",
                mode.distance23, " ",
                mode.distance13, " ",
                mode.maximum_distance, " ",
                mode.tolerance
            )
        end
    end
end

# ============================================================
# Spectrum filtering and matching
# ============================================================

function finite_sorted_modes(values)
    modes = CBF[]

    for omega in values
        if isfinite(real(omega)) && isfinite(imag(omega))
            push!(modes, CBF(omega))
        end
    end

    sort!(
        modes,
        by = omega -> (
            real(omega),
            imag(omega)
        )
    )

    return modes
end

function positive_imaginary_modes(
    values;
    threshold::BigFloat = INSTABILITY_THRESHOLD
)
    modes = CBF[]

    for omega in values
        if isfinite(real(omega)) &&
           isfinite(imag(omega)) &&
           imag(omega) > threshold

            push!(modes, CBF(omega))
        end
    end

    sort!(
        modes,
        by = omega -> (
            real(omega),
            imag(omega)
        )
    )

    return modes
end

function nearest_unused_mode(
    target::CBF,
    modes::Vector{CBF},
    used::BitVector
)
    best_index = 0
    best_distance = BigFloat(Inf)

    for index in eachindex(modes)
        if used[index]
            continue
        end

        distance = abs(modes[index] - target)

        if distance < best_distance
            best_distance = distance
            best_index = index
        end
    end

    return best_index, best_distance
end

function convergence_tolerance(
    omega1::CBF,
    omega2::CBF,
    omega3::CBF
)
    scale = max(
        BigFloat(1),
        abs(omega1),
        abs(omega2),
        abs(omega3)
    )

    return max(
        EPS_ABS,
        EPS_REL * scale
    )
end

function match_three_resolutions(
    values1,
    values2,
    values3
)
    modes1 = finite_sorted_modes(values1)
    modes2 = finite_sorted_modes(values2)
    modes3 = finite_sorted_modes(values3)

    used2 = falses(length(modes2))
    used3 = falses(length(modes3))

    accepted = NamedTuple[]

    for omega1 in modes1
        index2, distance12 = nearest_unused_mode(
            omega1,
            modes2,
            used2
        )

        if index2 == 0
            continue
        end

        omega2 = modes2[index2]

        index3, distance23 = nearest_unused_mode(
            omega2,
            modes3,
            used3
        )

        if index3 == 0
            continue
        end

        omega3 = modes3[index3]

        distance13 = abs(omega1 - omega3)

        maximum_distance = max(
            distance12,
            distance23,
            distance13
        )

        tolerance = convergence_tolerance(
            omega1,
            omega2,
            omega3
        )

        if maximum_distance <= tolerance
            used2[index2] = true
            used3[index3] = true

            omega_average = (
                omega1 +
                omega2 +
                omega3
            ) / BigFloat(3)

            push!(
                accepted,
                (
                    omega1 = omega1,
                    omega2 = omega2,
                    omega3 = omega3,
                    omega_average = omega_average,
                    distance12 = distance12,
                    distance23 = distance23,
                    distance13 = distance13,
                    maximum_distance = maximum_distance,
                    tolerance = tolerance
                )
            )
        end
    end

    sort!(
        accepted,
        by = mode -> (
            real(mode.omega_average),
            imag(mode.omega_average)
        )
    )

    return accepted
end

# ============================================================
# Solve one resolution
# ============================================================

function solve_resolution(
    Nres::Int,
    Mres::Int,
    J::Real,
    Jtag::Int,
    output_directory::String
)
    # Explicitly set the BigFloat precision inside each thread.
    return setprecision(BigFloat, BIGFLOAT_BITS) do
        println()
        println(
            "[thread $(threadid())] Starting " *
            "N=$(Nres), M=$(Mres), J=$(J), suffix=$(Jtag)"
        )

        M0, M1, M2 = load_qep_matrices(
            Nres,
            Mres,
            Jtag;
            dir = DATA_DIRECTORY
        )

        eigenvalues = qep_eigenvalues(M0, M1, M2)

        raw_file = joinpath(
            output_directory,
            "raw-modes-N$(Nres)-M$(Mres)-J$(Jtag).txt"
        )

        save_complex_values(
            raw_file,
            eigenvalues;
            header = (
                "# Raw eigenvalues for " *
                "N=$(Nres), M=$(Mres), J=$(J), suffix=$(Jtag)"
            )
        )

        println(
            "[thread $(threadid())] Completed " *
            "N=$(Nres), M=$(Mres), J=$(J)"
        )

        return eigenvalues
    end
end

# ============================================================
# Solve and validate one spin value
# ============================================================

function solve_spin(J::Real)
    Jtag = spin_tag(J)

    output_directory = joinpath(
        OUTPUT_DIRECTORY,
        "J$(Jtag)"
    )

    mkpath(output_directory)

    println()
    println("============================================================")
    println("Starting spin J=$(J), input-file suffix=$(Jtag)")
    println("============================================================")

    tasks = Dict{Int,Task}()

    # Launch the three spectral resolutions concurrently.
    for Nres in RESOLUTIONS
        let Nlocal = Nres
            tasks[Nlocal] = Threads.@spawn begin
                solve_resolution(
                    Nlocal,
                    ANGULAR_RESOLUTION,
                    J,
                    Jtag,
                    output_directory
                )
            end
        end
    end

    results = Dict{Int,Vector{CBF}}()

    # Wait for all three resolutions to finish.
    for Nres in RESOLUTIONS
        results[Nres] = fetch(tasks[Nres])
    end

    N1, N2, N3 = RESOLUTIONS

    spectrum1 = results[N1]
    spectrum2 = results[N2]
    spectrum3 = results[N3]

    # --------------------------------------------------------
    # Validate the complete spectra
    # --------------------------------------------------------

    stable_modes = match_three_resolutions(
        spectrum1,
        spectrum2,
        spectrum3
    )

    save_average_modes(
        joinpath(
            output_directory,
            "stable-modes-J$(Jtag).txt"
        ),
        stable_modes
    )

    save_mode_diagnostics(
        joinpath(
            output_directory,
            "stable-modes-diagnostics-J$(Jtag).txt"
        ),
        stable_modes;
        J = J,
        Jtag = Jtag
    )

    # --------------------------------------------------------
    # Save converged damped QNMs
    # --------------------------------------------------------

    stable_damped_modes = [
        mode for mode in stable_modes
        if imag(mode.omega_average) < 0
    ]

    save_average_modes(
        joinpath(
            output_directory,
            "stable-damped-modes-J$(Jtag).txt"
        ),
        stable_damped_modes
    )

    save_mode_diagnostics(
        joinpath(
            output_directory,
            "stable-damped-modes-diagnostics-J$(Jtag).txt"
        ),
        stable_damped_modes;
        J = J,
        Jtag = Jtag
    )

    # --------------------------------------------------------
    # Positive-imaginary-part diagnostic
    # --------------------------------------------------------

    positive1 = positive_imaginary_modes(spectrum1)
    positive2 = positive_imaginary_modes(spectrum2)
    positive3 = positive_imaginary_modes(spectrum3)

    save_complex_values(
        joinpath(
            output_directory,
            "positive-candidates-N$(N1)-J$(Jtag).txt"
        ),
        positive1
    )

    save_complex_values(
        joinpath(
            output_directory,
            "positive-candidates-N$(N2)-J$(Jtag).txt"
        ),
        positive2
    )

    save_complex_values(
        joinpath(
            output_directory,
            "positive-candidates-N$(N3)-J$(Jtag).txt"
        ),
        positive3
    )

    stable_positive_modes = match_three_resolutions(
        positive1,
        positive2,
        positive3
    )

    save_average_modes(
        joinpath(
            output_directory,
            "stable-positive-imaginary-modes-J$(Jtag).txt"
        ),
        stable_positive_modes
    )

    save_mode_diagnostics(
        joinpath(
            output_directory,
            "stable-positive-imaginary-diagnostics-J$(Jtag).txt"
        ),
        stable_positive_modes;
        J = J,
        Jtag = Jtag
    )

    println()
    println("Completed J=$(J), suffix=$(Jtag)")
    println("Stable modes: ", length(stable_modes))
    println("Stable damped modes: ", length(stable_damped_modes))
    println(
        "Stable positive-imaginary candidates: ",
        length(stable_positive_modes)
    )
    println("Results directory: ", output_directory)

    # Free the large matrices and eigenspectra before the next spin.
    empty!(results)
    GC.gc()

    return nothing
end

# ============================================================
# Main calculation
# ============================================================


println("Julia threads available: ", nthreads())

if nthreads() < length(RESOLUTIONS)
    println(
        "Warning: only $(nthreads()) Julia thread(s) are available. " *
        "The three resolutions for each spin will share the available threads."
    )
end

# Each spin launches one task per radial resolution. Keep enough spin values
# active to fill all Julia threads, while bounding the number of large spectra
# held in memory at once.
const CONCURRENT_SPINS = min(
    length(SPIN_VALUES),
    max(1, cld(nthreads(), length(RESOLUTIONS)))
)

println("Concurrent spin values: ", CONCURRENT_SPINS)
println(
    "Maximum resolution tasks in flight: ",
    CONCURRENT_SPINS * length(RESOLUTIONS)
)

next_spin = Threads.Atomic{Int}(1)

@sync for worker_id in 1:CONCURRENT_SPINS
    Threads.@spawn begin
        while true
            # atomic_add! returns the value before incrementing it.
            spin_index = Threads.atomic_add!(next_spin, 1)

            if spin_index > length(SPIN_VALUES)
                break
            end

            J = SPIN_VALUES[spin_index]

            try
                println(
                    "[spin worker $(worker_id)] Assigned J=$(J), " *
                    "suffix=$(spin_tag(J))"
                )

                solve_spin(J)
            catch error_value
                println(stderr)
                println(
                    stderr,
                    "Calculation failed for J=$(J), suffix=$(spin_tag(J))."
                )

                showerror(
                    stderr,
                    error_value,
                    catch_backtrace()
                )

                println(stderr)
                GC.gc()
            end
        end
    end
end

println()
println("All spin values have been processed.")
