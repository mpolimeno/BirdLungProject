module DataParsing

export extract_xcoordinatedata_from_COMSOL
export extract_ycoordinatedata_from_COMSOL
export extract_horizontalvelocity_data_from_COMSOL
export extract_verticalvelocity_data_from_COMSOL
export extract_horizontalacceleration_data_from_COMSOL
export extract_verticalacceleration_data_from_COMSOL

# get x-coordinate from file
function extract_xcoordinatedata_from_COMSOL(
    filename::String,
    t_initial::Float64,
    time_step::Float64,
    t_final::Float64,
    number_of_spatial_points::Int
)

    # COMSOL saves the relevant data in a single line
    # we will restructure it as an array that is n_nodes x n_timesteps
    timeVector = collect(t_initial:time_step:t_final)
    nTimes = length(timeVector)
    nNodes = number_of_spatial_points

    x_coordinateData = zeros(nNodes, nTimes)

    open(filename, "r") do io
        lines = eachline(io)
        for line in lines
            line = strip(line)

            if startswith(line, "% x (m) @ t=")
                # Extract time value
                tVal = parse(Float64, split(line, "=")[end])

                # Find matching time index
                idx = findfirst(t -> abs(t - tVal) < 1e-12, timeVector)
                idx === nothing && continue

                # Read next non-comment, non-empty line
                dataLine = ""
                while isempty(dataLine) && !eof(io)
                    dataLine = strip(readline(io))
                    if startswith(dataLine, "%")
                        dataLine = ""
                    end
                end

                # Parse x_coordinate values
                xVals = parse.(Float64, split(dataLine))

                if length(xVals) != nNodes
                    error("Expected $nNodes x_coordinate values at t = $tVal, got $(length(xVals)).")
                end

                x_coordinateData[:, idx] .= xVals
            end
        end
    end

    return x_coordinateData
end

#get y-coordinate from file
function extract_ycoordinatedata_from_COMSOL(
    filename::String,
    t_initial::Float64,
    time_step::Float64,
    t_final::Float64,
    number_of_spatial_points::Int
)

    # COMSOL saves the relevant data in a single line
    # we will restructure it as an array that is n_nodes x n_timesteps
    timeVector = collect(t_initial:time_step:t_final)
    nTimes = length(timeVector)
    nNodes = number_of_spatial_points

    y_coordinateData = zeros(nNodes, nTimes)

    open(filename, "r") do io
        lines = eachline(io)
        for line in lines
            line = strip(line)

            if startswith(line, "% y (m) @ t=")
                # Extract time value
                tVal = parse(Float64, split(line, "=")[end])

                # Find matching time index
                idx = findfirst(t -> abs(t - tVal) < 1e-12, timeVector)
                idx === nothing && continue

                # Read next non-comment, non-empty line
                dataLine = ""
                while isempty(dataLine) && !eof(io)
                    dataLine = strip(readline(io))
                    if startswith(dataLine, "%")
                        dataLine = ""
                    end
                end

                # Parse y_coordinate values
                xVals = parse.(Float64, split(dataLine))

                if length(xVals) != nNodes
                    error("Expected $nNodes y_coordinate values at t = $tVal, got $(length(xVals)).")
                end

                y_coordinateData[:, idx] .= xVals
            end
        end
    end

    return y_coordinateData
end

# get horizontal velocity u from file
function extract_horizontalvelocity_data_from_COMSOL(
    filename::String,
    t_initial::Float64,
    time_step::Float64,
    t_final::Float64,
    number_of_spatial_points::Int
)

    # COMSOL saves the relevant data in a single line
    # we will restructure it as an array that is n_nodes x n_timesteps
    timeVector = collect(t_initial:time_step:t_final)
    nTimes = length(timeVector)
    nNodes = number_of_spatial_points

    u_velocityData = zeros(nNodes, nTimes)

    open(filename, "r") do io
        lines = eachline(io)
        for line in lines
            line = strip(line)

            if startswith(line, "% u (m/s) @ t=")
                # Extract time value
                tVal = parse(Float64, split(line, "=")[end])

                # Find matching time index
                idx = findfirst(t -> abs(t - tVal) < 1e-12, timeVector)
                idx === nothing && continue

                # Read next non-comment, non-empty line
                dataLine = ""
                while isempty(dataLine) && !eof(io)
                    dataLine = strip(readline(io))
                    if startswith(dataLine, "%")
                        dataLine = ""
                    end
                end

                # Parse u_velocity values
                uVals = parse.(Float64, split(dataLine))

                if length(uVals) != nNodes
                    error("Expected $nNodes u_velocity values at t = $tVal, got $(length(uVals)).")
                end

                u_velocityData[:, idx] .= uVals
            end
        end
    end

    return u_velocityData
end

# get vertical velocity v
function extract_verticalvelocity_data_from_COMSOL(
    filename::String,
    t_initial::Float64,
    time_step::Float64,
    t_final::Float64,
    number_of_spatial_points::Int
)

    # COMSOL saves the relevant data in a single line
    # we will restructure it as an array that is n_nodes x n_timesteps
    timeVector = collect(t_initial:time_step:t_final)
    nTimes = length(timeVector)
    nNodes = number_of_spatial_points

    v_velocityData = zeros(nNodes, nTimes)

    open(filename, "r") do io
        lines = eachline(io)
        for line in lines
            line = strip(line)

            if startswith(line, "% v (m/s) @ t=")
                # Extract time value
                tVal = parse(Float64, split(line, "=")[end])

                # Find matching time index
                idx = findfirst(t -> abs(t - tVal) < 1e-12, timeVector)
                idx === nothing && continue

                # Read next non-comment, non-empty line
                dataLine = ""
                while isempty(dataLine) && !eof(io)
                    dataLine = strip(readline(io))
                    if startswith(dataLine, "%")
                        dataLine = ""
                    end
                end

                # Parse v_velocity values
                vVals = parse.(Float64, split(dataLine))

                if length(vVals) != nNodes
                    error("Expected $nNodes v_velocity values at t = $tVal, got $(length(vVals)).")
                end

                v_velocityData[:, idx] .= vVals
            end
        end
    end

    return v_velocityData
end

# get pressure from file
function extract_pressuredata_from_COMSOL(
    filename::String,
    t_initial::Float64,
    time_step::Float64,
    t_final::Float64,
    number_of_spatial_points::Int
)

    # COMSOL saves the relevant data in a single line
    # we will restructure it as an array that is n_nodes x n_timesteps
    timeVector = collect(t_initial:time_step:t_final)
    nTimes = length(timeVector)
    nNodes = number_of_spatial_points

    pressureData = zeros(nNodes, nTimes)

    open(filename, "r") do io
        lines = eachline(io)
        for line in lines
            line = strip(line)

            if startswith(line, "% p (Pa) @ t=")
                # Extract time value
                tVal = parse(Float64, split(line, "=")[end])

                # Find matching time index
                idx = findfirst(t -> abs(t - tVal) < 1e-12, timeVector)
                idx === nothing && continue

                # Read next non-comment, non-empty line
                dataLine = ""
                while isempty(dataLine) && !eof(io)
                    dataLine = strip(readline(io))
                    if startswith(dataLine, "%")
                        dataLine = ""
                    end
                end

                # Parse pressure values
                pVals = parse.(Float64, split(dataLine))

                if length(pVals) != nNodes
                    error("Expected $nNodes pressure values at t = $tVal, got $(length(pVals)).")
                end

                pressureData[:, idx] .= pVals
            end
        end
    end

    return pressureData
end

# get horizontal acceleration ut
function extract_horizontalacceleration_data_from_COMSOL(
    filename::String,
    t_initial::Float64,
    time_step::Float64,
    t_final::Float64,
    number_of_spatial_points::Int
)

    # COMSOL saves the relevant data in a single line
    # we will restructure it as an array that is n_nodes x n_timesteps
    timeVector = collect(t_initial:time_step:t_final)
    nTimes = length(timeVector)
    nNodes = number_of_spatial_points

    u_accelerationData = zeros(nNodes, nTimes)

    open(filename, "r") do io
        lines = eachline(io)
        for line in lines
            line = strip(line)

            if startswith(line, "% ut (m/s^2) @ t=")
                # Extract time value
                tVal = parse(Float64, split(line, "=")[end])

                # Find matching time index
                idx = findfirst(t -> abs(t - tVal) < 1e-12, timeVector)
                idx === nothing && continue

                # Read next non-comment, non-empty line
                dataLine = ""
                while isempty(dataLine) && !eof(io)
                    dataLine = strip(readline(io))
                    if startswith(dataLine, "%")
                        dataLine = ""
                    end
                end

                # Parse u_acceleration values
                utVals = parse.(Float64, split(dataLine))

                if length(utVals) != nNodes
                    error("Expected $nNodes u_acceleration values at t = $tVal, got $(length(utVals)).")
                end

                u_accelerationData[:, idx] .= utVals
            end
        end
    end

    return u_accelerationData
end

# get vertical acceleration
function extract_verticalacceleration_data_from_COMSOL(
    filename::String,
    t_initial::Float64,
    time_step::Float64,
    t_final::Float64,
    number_of_spatial_points::Int
)

    # COMSOL saves the relevant data in a single line
    # we will restructure it as an array that is n_nodes x n_timesteps
    timeVector = collect(t_initial:time_step:t_final)
    nTimes = length(timeVector)
    nNodes = number_of_spatial_points

    v_accelerationData = zeros(nNodes, nTimes)

    open(filename, "r") do io
        lines = eachline(io)
        for line in lines
            line = strip(line)

            if startswith(line, "% vt (m/s^2) @ t=")
                # Extract time value
                tVal = parse(Float64, split(line, "=")[end])

                # Find matching time index
                idx = findfirst(t -> abs(t - tVal) < 1e-12, timeVector)
                idx === nothing && continue

                # Read next non-comment, non-empty line
                dataLine = ""
                while isempty(dataLine) && !eof(io)
                    dataLine = strip(readline(io))
                    if startswith(dataLine, "%")
                        dataLine = ""
                    end
                end

                # Parse v_acceleration values
                vtVals = parse.(Float64, split(dataLine))

                if length(vtVals) != nNodes
                    error("Expected $nNodes v_acceleration values at t = $tVal, got $(length(vtVals)).")
                end

                v_accelerationData[:, idx] .= vtVals
            end
        end
    end

    return v_accelerationData
end

end