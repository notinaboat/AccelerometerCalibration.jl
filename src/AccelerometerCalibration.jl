"""
# AccelerometerCalibration.jl

Uses [Optim.jl](https://github.com/JuliaNLSolvers/Optim.jl) to fit calibration
parameters (offset, scale, rotation) to a stream of accelerometer data.

## Examples

Single channel:

```julia
cal = Calibration()
for s in sensor.get_samples()
    update_calibration!(cal, s)
    push!(calibrated_samples, cal * s)
end
```

Multiple channels
(assumed to be fixed together and moved in unison for calibration):

```julia
calv = [Calibration() for i in 1:4]
while true
    sv = [sensor1.get_sample(),
          sensor2.get_sample(),
          sensor3.get_sample(),
          sensor4.get_sample()]
    update_calibration!(calv, vs)
    push!(calibrated_samples, [cal * s for s in sv])
end
```

![](test/plot.gif)
"""
module AccelerometerCalibration

using Statistics: mean
using StaticArrays
using DataStructures: CircularBuffer, isfull
using Rotations
using Optim


"Threshold for stable input signal (unit g)."
const stable_th = 0.2

"Threshold for separation between calibration points (unit g))."
const cal_th = 0.4

"Averaging window size."
const window_size = 10

"Time limit for calibration update (unit s)."
const default_time_limit = 1/10


"""
Accelerometer Data Calibration

A sliding window is used to find periods of stable input (sensor not moving).

A cloud of points measured at different sensor orientations is used to
estimate (x,y,z) zero offset and scaling error.

Points are expected to have magnitude ~1g, so the calibration parameters
are chosen to minimise: `sum((1 - hypot(x, y, z))^2 for (x, y, z) in points`.

See [`evaluate_point`](@ref)
"""
mutable struct Calibration
    window::CircularBuffer{SVector{3,Float64}}
    points::Vector{SVector{3,Float64}}
    offset::SVector{3, Float64}
    scale::SVector{3, Float64}
    rotation::RotXYZ
    time_limit::Float64
    callback::Function
    Calibration() = new(CircularBuffer{SVector{3,Float64}}(window_size),
                        Vector{SVector{3, Float64}}(),
                        @SVector(zeros(3)),
                        @SVector(ones(3)),
                        one(RotXY),
                        default_time_limit,
                        ()->nothing)
end

import Base: *
*(c::Calibration, v) = c.rotation * ((v .* c.scale) - c.offset)

Base.length(c::Calibration) = length(c.points)

function Base.empty!(c::Calibration)
    empty!(c.window)
    empty!(c.points)
    empty!(c.points)
    c.offset = zeros(3)
    c.scale = ones(3)
end

function Base.push!(c::Calibration, p)
    push!(c.points, p)
    if length(c.points) > 2
        calibrate!(c)
    end
    nothing
end


"""
Find (x,y,z) offset and scaling that best fits a unit-sphere to points.
"""
function calibrate!(cal; time_limit=cal.time_limit)

    @assert length(cal.points) > 2

    # Prepare parameter vector `p` and upper and lower bounds...
    p = zeros(3 + 3)
    p_offset = 1:3
    p_scale = 4:6
    lower = [-1.0,-1.0,-1.0,  0.9, 0.9, 0.9]
    upper = [ 1.0, 1.0, 1.0,  1.1, 1.1, 1.1]

    # Load current paramter values.
    p[p_offset] = cal.offset
    p[p_scale] = cal.scale

    # Function to minimize.
    f = p->sum((1 - hypot(((v .* p[p_scale]) - p[p_offset])...))^2
              for v in cal.points)

    # Run optimisation.
    r = optimize(f, p, ParticleSwarm(lower,upper,10),
                       Optim.Options(time_limit = time_limit))
    p = Optim.minimizer(r)

    # Update calibration.
    cal.offset = p[p_offset]
    cal.scale = p[p_scale]

    cal.callback()
    nothing
end


"""
Find rotation that best aligns multiple channels.
"""
function calibrate_rotation!(cals; time_limit=cals[1].time_limit)

    l = length(cals[1])

    target_points = [mean(c.points[i] for c in cals) for i in 1:l]

    for c in cals
        # Prepare parameter vector `p` and upper and lower bounds...
        r = c.rotation
        p = [r.theta1, r.theta2, r.theta3]
        lower = [-π/16,-π/16,-π/16]
        upper = [ π/16, π/16, π/16]

        # Apply scaling and offset to calibration points.
        cal_points = [((p .* c.scale) - c.offset) for p in c.points]

        # Function to minimize.
        err = (i, p) -> target_points[i] - (RotXYZ(p...) * cal_points[i])
        f = (p) -> mapreduce(x->x^2, +, vcat((err(i, p) for i in 1:l)...))

        # Run optimisation.
        r = optimize(f, p, ParticleSwarm(lower,upper,3),
                           Optim.Options(time_limit = time_limit))
        p = Optim.minimizer(r)

        # Update calibration.
        c.rotation = RotXYZ(p...)

        c.callback()
    end
end


"""
Does `cal.window` contain useful new calibration data?
"""
function buffer_contains_new_point(cal)

    # Store point in window buffer and wait until window is full.
    w = cal.window
    if !isfull(w)
        return false
    end
    p = mean(w)

    #  All samples in the window must be close to the mean of the window."
    if any(v->mapreduce(x->x^2, +, v) > stable_th^2, w .- [p])
        return false
    end

    #  New calibration point must not be too close to an existirng one.
    if any(q->hypot((p - q)...) < cal_th, cal.points)
        return false
    end

    return true
end


"""
    update_calibration!(::Calibration, sample)

Update calibration with a with a new raw accelerometer sample.
Call `Calibration.callback` if the calibration is changed.
"""
update_calibration!(cal::Calibration, p) = update_calibration!([cal], [p])

"""
    update_calibration!([::Calibration, ...], [sample, ...])

Update a Vector of Calibrations from a vector of samples.
Calibrate rotation to align the samples.
"""
function update_calibration!(cal, p)
    buffer_point!.(cal, p)
    if all(buffer_contains_new_point.(cal))
        push!.(cal, p)
        calibrate_rotation!(cal)
    end
    nothing
end

buffer_point!(cal, p) = push!(cal.window, p)



end # module
