# Extract xz compressed tar (.txz)
function unpack_txz(fname::AbstractString, dir::AbstractString)
    txz = open(fname)
    tar = XzDecompressorStream(txz)
    Tar.extract(tar, String(dir))
    close(tar)
end

# Given the filename of a MIST BC set, return a vector containing
# the column labels to be used when unpacking
# const mist_file_info = Dict("CFHTugriz.txz" => (name="CFHT_MegaCam",
#                                                 header=["Teff", "logg", "FeH", "Av", "Rv", "u", "CaHK", "g", "r", "i_new", "i_old", "z"]))

# Read data from unpacked MIST BC file into Table
function read_mist_bc(fname::AbstractString, header) # ::AbstractVector{<:AbstractString})
    return CSV.read(fname, Table;
                    comment="#", ignorerepeated=true, delim=' ', header=header)
end

function parse_mist_header(fname::AbstractString)
    for (linenum, line) in enumerate(eachline(fname))
        if linenum == 6
            # CFHT has logg and [Fe/H] run into each other, so we
            # let first part be constant, only read special the filters in the header
            # header = filter(!isempty, split(split(line, "Rv")[2], ' '))
            # return vcat(["Teff", "logg", "feh", "Av", "Rv"], convert(Vector{String}, header))
            header = filter(!isempty, split(split(line, string(last(_mist_dependents)))[2], ' '))
            return vcat(SVector(string.(_mist_dependents)), convert(Vector{String}, header))
        end
    end
end

# Get [Fe/H] from name of MIST BC file
function mist_feh(fname::AbstractString)
    fname = basename(fname) # Get file name part of path
    feh = parse(Int, fname[5:7]) / 100
    if fname[4] == 'm'
        feh *= -1
    end
    return feh
end

function custom_unpack(fname::AbstractString) # Path to datadep
    # println(fname)
    # ~/.julia/scratchspaces/124859b0-ceae-595e-8997-d05f6a7a8dfe/datadeps
    # the directory in scratchspaces is the UUID of datadeps
    fpath = dirname(fname)
    out_dir = joinpath(fpath, "files")
    if isdir(out_dir)
        rm(out_dir; force=true, recursive=true)
    end
    unpack_txz(fname, out_dir)
    
    # Process into a single compressed CSV for easy reads
    # Use name of datadep (last directory of path) for file name
    bfname = splitpath(fname)[end-1] 
    all_files = readdir(out_dir; join=true)
    # Process files in sorted order of [Fe/H] to ensure final file
    # has sorted [Fe/H]
    feh_vals = mist_feh.(all_files)
    idxs = sortperm(feh_vals) 
    header = parse_mist_header(first(all_files))
    bigtable = reduce(vcat, read_mist_bc(all_files[i], header) for i in idxs)
    CSV.write(joinpath(fpath, splitext(bfname)[1]*".gz"), bigtable; compress=true)
    # Done with original .txz and extracted files so remove
    rm(fname)
    rm(out_dir; force=true, recursive=true)
    
    # # Process into single HDF5 file for easy future reads
    # bfname = basename(fname) # Get the file name part of a path
    # # file_info = mist_file_info[bfname]
    # all_files = readdir(out_dir; join=true)
    # header = parse_mist_header(first(all_files))
    # # Find first element of file_info that is a photometric filter;
    # # comes after Rv column
    # filter_index = findfirst(==("Rv"), header) + 1
    # # Open HDF5 file to write in data
    # fid = HDF5.h5open(joinpath(fpath, splitext(bfname)[1]*".hdf5"), "w")
    # # Loop over each file (which indexes [Fe/H]) and load
    # for file in all_files
    #     feh = mist_feh(file) # Get [Fe/H] from file name
    #     table = read_mist_bc(file, header)
    #     @argcheck feh ≈ first(table.feh) # Test that filename [Fe/H] agrees with table
    #     # Create group for this [Fe/H]
    #     g = HDF5.create_group(fid, @sprintf("%.2f", feh))
    #     # rvs = [@sprintf("%.2f", Rv) for Rv in unique(table.Rv)]
    #     # Avs = [@sprintf("%.2f", Av) for Av in unique(table.Av)]
    #     Avs = unique(table.Av)
    #     for Av in Avs
    #         subtable = filter(row -> row.Av == Av, table)
    #         println(length(subtable))
    #         h = HDF5.create_group(g, @sprintf("%.2f", Av))
    #         h["Teff", deflate=3] = subtable.Teff
    #         h["logg", deflate=3] = subtable.logg
    #         h["data", deflate=3] = Tables.matrix(subtable)[:, filter_index:end]
    #     end
    # end
    # HDF5.close(fid)
    
end



# See https://waps.cfa.harvard.edu/MIST/model_grids.html#bolometric 
"Available filter sets for MIST bolometric corrections."
const BC_sets = ("CFHT_MegaCam", "HST_ACS_HRC", "HST_ACS_WFC", "HST_WFPC2", "HST_WFC3", "JWST")

function __init__()
    register(DataDep("MIST_CFHT_MegaCam", "MIST bolometric corrections for CFHT/MegaCam",
                     "https://waps.cfa.harvard.edu/MIST/BC_tables/CFHTugriz.txz",
                     "cf10e060b65591cb43fc55544260a75ec559c28d9454a4f681a197af26aa0fcd";
                     post_fetch_method=custom_unpack))
    register(DataDep("MIST_HST_ACS_HRC", "MIST bolometric corrections for HST ACS/HRC",
                     "https://waps.cfa.harvard.edu/MIST/BC_tables/HST_ACSHR.txz",
                     "065e5235a55b24541a027000f358a08036e7305fb105865497e8e8c452fa512d";
                     post_fetch_method=custom_unpack))
    register(DataDep("MIST_HST_ACS_WFC", "MIST bolometric corrections for HST ACS WFC",
                     "https://waps.cfa.harvard.edu/MIST/BC_tables/HST_ACSWF.txz",
                     "5e48f529d245dae122f3ac86c6524f34cf6b0952816cd4a567d78c41e3439571";
                     post_fetch_method=custom_unpack))
    register(DataDep("MIST_HST_WFPC2", "MIST bolometric corrections for HST WFPC2",
                     "https://waps.cfa.harvard.edu/MIST/BC_tables/HST_WFPC2.txz",
                     "941f18684c741979d53a6b4a4083a8cd913a6fbefe03c5bfc5ee179590d14410";
                     post_fetch_method=custom_unpack))
    register(DataDep("MIST_HST_WFC3", "MIST bolometric corrections for HST WFC3",
                     "https://waps.cfa.harvard.edu/MIST/BC_tables/HST_WFC3.txz",
                     "39fdf9a2606b6f66f5388c2a02bdb18bb1c32402325484c917b5cd564bcd124d";
                     post_fetch_method=custom_unpack))
    register(DataDep("MIST_JWST", "MIST bolometric corrections for JWST (pre-launch, not updated)",
                     "https://waps.cfa.harvard.edu/MIST/BC_tables/JWST.txz",
                     "580d13c2dfe6086f8593828170f573fd653f6a51051d0437b602034cb5b49779";
                     post_fetch_method=custom_unpack))

    # Load datadeps for development
    # datadep"MIST_CFHT_MegaCam"
    # datadep"MIST_HST_ACS_HRC"
    # datadep"MIST_HST_ACS_WFC"
    # datadep"MIST_HST_WFPC2"
    # datadep"MIST_HST_WFC3"
    datadep"MIST_JWST"

end
                     
