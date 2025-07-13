using Git: git
using Scratch: @get_scratch!

const ybc_url = "https://gitlab.com/cycyustc/ybc_tables.git"
# Looks like raw spectra might be hosted here
# https://gitlab.com/cycyustc/spec_ybc

scratch_dir = ""
ybc_path = ""
"""Tracks git sparse checkout list so we don't have to call slow git functions (`git sparse-checkout list`)."""
const sparse_checkout_list = Vector{String}()
"""Tracks the valid photometric filter systems for YBC."""
const systems = Vector{String}()

function __init__()
    global scratch_dir = @get_scratch!("YBC") # This will create the YBC dir
    global ybc_path = _ybc_path(scratch_dir)
    for i in ("iYBC", "rYBC", "YBC")
        # mkdir errors if dir exists, mkpath does not
        mkpath(joinpath(ybc_path, i))
    end
    # Read the Git sparse checkout list and write contents into sparse_checkout_list
    for f in split(readchomp(`$(git()) -C $ybc_path sparse-checkout list`), "\n")
        push!(sparse_checkout_list, f)
    end
    # Read subdirectories of ybc_tables/YBC and write into systems
    for f in split(readchomp(`$(git()) -C $(joinpath(ybc_path, "YBC")) ls-tree -d master --name-only`), "\n") # $(joinpath(repo, prefix))
        push!(systems, f)
    end
end

"""
    _ybc_path(path::String = scratch_dir)

Ensure that `path` is initialized properly to contain the YBC git repository. Returns path to the directory containing the initialized YBC Git repository.
"""
function _ybc_path(path::String = scratch_dir)
    # Path where git repo will be after checkout
    full_path = joinpath(path, "ybc_tables")
    # Check if directory is already initialized -- if so, check that 
    # correct repo is checked out
    # We can use git -C <path> to perform operations in other directories, 
    # rather than changing directories.
    try
        # First check that path is in a git tree
        cmd = `$(git()) -C $full_path rev-parse --is-inside-work-tree`
        run(pipeline(cmd; stdout=devnull, stderr=devnull))
        # Now check that it is reflecting the correct directory
        url = readchomp(`$(git()) -C $full_path remote get-url origin`)
        if url != ybc_url
            # # If wrong directory checked out, delete .git configuration
            # rm(joinpath(path, ".git"))
            error("Wrong repository $url checked out in YBC scratch directory $path; expected url $ybc_url. If $path is an acceptable place to store the YBC files, please remove its `.git` subdirectory to enable a fresh clone.")
        end
    catch
        # Clone repo skipping file contents (--filter=blob:none) and not checking out
        @info "Initializing YBC repository in $path"
        run(`$(git()) -C $path clone --filter=blob:none --no-checkout $ybc_url`)
        run(`$(git()) -C $(joinpath(path, "ybc_tables")) lfs install --local`)

        # Perform a sparse checkout (only `add`ed files will be checked out),
        # interpreting sparse paths as entire directories, rather than complex patterns (--cone)
        run(`$(git()) -C $(joinpath(path, "ybc_tables")) sparse-checkout init --cone`)
        # Now checkout
        run(`$(git()) -C $(joinpath(path, "ybc_tables")) checkout master`)
    end
    return full_path
end

function update_tables(path::String = ybc_path)
    run(`$(git()) -C $path fetch origin`)
    run(`$(git()) -C $path pull origin master`)
    return nothing
end

"""
    pull_table(f::AbstractString, prefix::AbstractString="YBC")

Pull files from YBC repository corresponding to filter system `f` which must correspond to a valid subdirectory in `joinpath(ybc_tables, prefix)`. Available `prefix` entries (as of 2025-07-09) are "YBC" (for standard BCs), "rYBC" (BC tables for rotating stars), and "iYBC" (the limb darkening coefficients with Kurucz libraries).

Returns absolute path to the pulled files.
"""
function pull_table(f::AbstractString, prefix::AbstractString = "YBC")
    f = String(f)
    prefix = String(prefix)
    # Get path to repo
    repo = ybc_path
    # Check that prefix is valid
    check_prefix(prefix)

    # Check that f is a supported system
    if f ∉ systems
        throw(ArgumentError("Requested filter system $f invalid; available systems are $systems."))
    end

    # Add requested filter system to sparse-checkout list
    name = joinpath(prefix, f)
    # Running sparse-checkout add is slow -- faster to check if entry already there,
    # then if not, run add command
    if name ∉ sparse_checkout_list
        run(`$(git()) -C $repo sparse-checkout add $name`)
        push!(sparse_checkout_list, name)
    end

    return joinpath(repo, prefix, f)
end

"""
    remove_table(f::AbstractString, prefix::AbstractString = "YBC")

Remove table `joinpath(prefix, f)` from the Git sparse-checkout list for YBC -- this will uninstall the related data files.
"""
function remove_table(f::AbstractString, prefix::AbstractString = "YBC")
    # Remove invalid system from sparse-checkout list
    # To do this need to `set` with list that does not include invalid entry
    # run(`$(git()) -C $repo sparse-checkout rm $(prefix * "/" * f)`)
    repo = ybc_path
    name = joinpath(prefix, f)
    # installed = split(readchomp(`$(git()) -C $repo sparse-checkout list`), "\n")
    idxs = sparse_checkout_list .== name
    # good = filter(!=(name), sparse_checkout_list)
    good = sparse_checkout_list[.~idxs]
    run(`$(git()) -C $repo sparse-checkout set $good`)
    run(`$(git()) -C $repo sparse-checkout reapply`)
    deleteat!(sparse_checkout_list, idxs)
    return nothing
end

