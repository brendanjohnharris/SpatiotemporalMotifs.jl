
dddir = datadir("calculations")
files = readdir(dddir, join = true)

map(files) do file
    f = jldopen(file, "a+")
    if haskey(f, "trials") && haskey(f, "rectified_change_times")
        try
            global trials = f["trials"]
        catch e
            return close(f)
        end
        if "rectified_change_times" in names(trials)
            return close(f)
        end
        @info "$file"
        rct = f["rectified_change_times"]
        trials.rectified_change_times = rct[2:(end - 2)]

        delete!(f, "trials")
        f["trials"] = trials
    end
    return close(f)
end
