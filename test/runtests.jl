using Test
using AccelerometerCalibration
using AccelerometerCalibration: Calibration, calibrate!, update_calibration!
using AccelerometerCalibration: calibrate_rotation!
using Serialization
using Plots
using AccelerometerCalibrationPlots

@testset "AccelerometerCalibration" begin

    # Load raw acceleromter recording...
    raw_points = deserialize("raw_points.serialized")
    channel_count = length(raw_points[1])

    # Create a Calibration for each channel.
    cal_v = [Calibration() for i in 1:channel_count]

    # Set up callback to re-plot when a new calibration point is found.
    anim = Animation()
    AccelerometerCalibrationPlots.reset()
    cal_v[4].callback = () -> begin
        AccelerometerCalibrationPlots.calplot(cal_v)
        frame(anim)
    end

    cal_error(cal) = sum((1 - hypot((cal * v)...))^2 for v in cal.points)

    for (i, sample) in enumerate(raw_points)
        global point_count
        update_calibration!(cal_v, sample)
        if i == 20
            @test any(x->x>0.2, cal_error.(cal_v))
        end
    end

    for cal in cal_v
        err = cal_error(cal)
        @show err
        @test err < 0.02
    end

    # Calibrate again without time limit
    for cal in cal_v
        calibrate!(cal, time_limit=Inf)
    end
    calibrate_rotation!(cal_v, time_limit=Inf)
    for cal in cal_v
        err = cal_error(cal)
        @show err
        @test err < 0.02
    end

    AccelerometerCalibrationPlots.calplot(cal_v)
    frame(anim)

    gif(anim, "plot.gif", fps=5)
end
