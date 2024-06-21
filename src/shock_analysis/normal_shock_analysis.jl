#TODO: develop this

function check_shock_consistency(shock_positions_over_time)
    time_steps = length(shock_positions_over_time)
    # Initialize the first frame shock
    detected_shock_timestep = -1
    glitch_counter = 0
    for (i, frame) in enumerate(shock_positions_over_time)
        if !isempty(frame)
            # Flag the first frame shock is detected
            if detected_shock_timestep != -1
                gap_size = i - detected_shock_timestep
                # The gap size that happens because of inconsistent shock detection
                # shouldn't be too small nor too large since shock wave can dissipate
                # and we shouldn't confuse this with inconsistency
                if gap_size >= 2 && gap_size < 10
                    glitch_counter += 1
                end
            end
            detected_shock_timestep = i
        end

        if time_steps < 100
            if glitch_counter > 5
                println("The shock detection algorithm is inconsistent")
                println("$glitch_counter")
                return false
            end
        else 
            if glitch_counter > max(round(0.05 * time_steps))
                println("The shock detection algorithm is inconsistent")
                println("$glitch_counter")
                return false
            end
        end
    end
    return true
end