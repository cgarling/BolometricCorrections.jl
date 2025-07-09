using Git: git
using Scratch: @get_scratch!

const ybc_url = "https://gitlab.com/cycyustc/ybc_tables.git"

scratch_dir = ""

function __init__()
    global scratch_dir = @get_scratch!("YBC") # This will create the YBC dir
end

"""
    ybc_path(path::String = scratch_dir)
Ensure that `path` is initialized properly to contain the YBC git repository. Returns path to the directory containing the initialized YBC Git repository.
"""
function ybc_path(path::String = scratch_dir)
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
        # Perform a sparse checkout (only `add`ed files will be checked out),
        # interpreting sparse paths as entire directories, rather than complex patterns (--cone)
        run(`$(git()) -C $(joinpath(path, "ybc_tables")) sparse-checkout init --cone`)
        # Now checkout
        run(`$(git()) -C $(joinpath(path, "ybc_tables")) checkout master`)
    end
    return full_path
end

function _clean_list(prefix::AbstractString="YBC")

end
"""
    pull_table(f::AbstractString, prefix::AbstractString="YBC")
Pull files from YBC repository corresponding to filter system `f` which must correspond to a valid subdirectory in `ybc_tables/prefix`. Available `prefix` entries (as of 2025-07-09) are "YBC" (for standard BCs), "rYBC" (BC tables for rotating stars), and "iYBC" (the limb darkening coefficients with Kurucz libraries.)
"""
function pull_table(f::AbstractString, prefix::AbstractString = "YBC")
    f = String(f)
    prefix = String(prefix)
    # Get path to repo
    repo = ybc_path()
    # Check that prefix is valid -- YBC, iYBC, and rYBC should be valid.
    # YBC is the normal BC tables, rYBC
    subdirs = split(readchomp(`$(git()) -C $repo ls-tree -d master --name-only`), "\n")
    if prefix ∉ subdirs
        throw(ArgumentError("prefix $prefix invalid; available prefixes are $subdirs."))
    end
    # Add requested filter system to sparse-checkout list
    run(`$(git()) -C $repo sparse-checkout add $(prefix * "/" * f)`)

    # Ensure that requested filter system `f` is in the remote directory list
    # Does not work until you add at least one filter to create the `prefix` folder,
    # so we'll issue the add command first, then check that the request was valid,
    # and if not we will remove the invalid entry.
    systems = split(readchomp(`$(git()) -C $(joinpath(repo, prefix)) ls-tree -d master --name-only`), "\n")
    if f ∉ systems
        # Remove invalid system from sparse-checkout list
        # To do this need to `set` with list that does not include invalid entry
        # run(`$(git()) -C $repo sparse-checkout rm $(prefix * "/" * f)`)
        installed = split(readchomp(`$(git()) -C $repo sparse-checkout list`), "\n")
        good = filter(!=(prefix*"/"*f), installed)
        run(`$(git()) -C $repo sparse-checkout set $good`)
        run(`$(git()) -C $repo sparse-checkout reapply`)
        throw(ArgumentError("Requested filter system $f invalid; available systems are $systems."))
    end
end