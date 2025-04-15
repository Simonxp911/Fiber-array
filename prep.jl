

"""
Prepare the atomic coherences for plotting (in the case of no phonons)
"""
function prep_times_σTrajectories(xTrajectories, N)
    times = xTrajectories[:, 1]
    
    σTrajectories = unpack_σFromx.(eachrow(xTrajectories[:, 2:end]))
    
    σTrajectories = [[σTrajectories[t][i] for t in eachindex(times)] for i in 1:N]
    return times, σTrajectories
end


"""
Prepare the atomic coherences and the atom-phonon correlations for plotting
"""
function prep_times_σBαTrajectories(xTrajectories, N)
    times = xTrajectories[:, 1]
    
    σBαTrajectories = unpack_σBαFromx.(eachrow(xTrajectories[:, 2:end]))
    
    σTrajectories  =  [[σBαTrajectories[t][1][i]    for t in eachindex(times)] for i in 1:N]
    BαTrajectories = [[[σBαTrajectories[t][2][α][i] for t in eachindex(times)] for i in 1:N^2] for α in 1:3]
    return times, σTrajectories, BαTrajectories
end