# This file contains code of various filters to smothen data

"""
    savitzky_golay_filter_53(data::Vector{Float64})

Applies a 5-point, 3rd-order Savitzky-Golay filter to the input data vector.
"""
function savitzky_golay_filter_53(data::Vector{Float64})

    # 5-point, 3rd-order Savitzky-Golay filter coefficients
    half_window = div(5, 2)
    coeffs = [-3, 12, 17, 12, -3] ./ 35  # Coefficients for cubic SG filter with window size 5

    filtered_data = similar(data)
    n = length(data)

    for i in 1:n
        acc = 0.0
        for j in -half_window:half_window
            idx = clamp(i + j, 1, n) # let the index stay within bounds
            acc += coeffs[j + half_window + 1] * data[idx]
        end
        filtered_data[i] = acc
    end

    return filtered_data
end

"""
    moving_average_filter(data::Vector{Float64}, window_size::Int)

Applies a moving average filter to the input data vector with the specified window size.
"""
function moving_average_filter(data::Vector{Float64}, window_size::Int)

    half_window = div(window_size, 2)
    filtered_data = similar(data)
    n = length(data)

    for i in 1:n
        acc = 0.0
        count = 0
        for j in -half_window:half_window
            idx = i + j
            if idx >= 1 && idx <= n
                acc += data[idx]
                count += 1
            end
        end
        filtered_data[i] = acc / count
    end

    return filtered_data
end
