# Extract xz compressed tar (.txz)
function unpack_txz(fname::AbstractString, dir::AbstractString)
    # txz = open(fname)
    # tar = XzDecompressorStream(txz)
    # Tar.extract(tar, String(dir))
    # close(tar)
    open(XzDecompressorStream, fname) do stream
        Tar.extract(stream, String(dir))
    end
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
const BC_sets = ("CFHT_MegaCam", "DECam", "GALEX", "HST_ACS_HRC", "HST_ACS_WFC",
                 "HST_WFC3", "HST_WFPC2", "JWST", "LSST", "PanStARRS", "SDSS",
                 "SkyMapper", "Subaru/HSC", "IPHAS", "Spitzer/IRAC", "SPLUS",
                 "Swift", "UBVRIplus", "UKIDSS", "UVIT", "VISTA", "Washington",
                 "WFIRST", "WISE")

function __init__()
    prefix = "https://waps.cfa.harvard.edu/MIST/BC_tables/"
    # Actually just distribute zeropoints with source
    # register(DataDep("MIST_zeropoints", "Zeropoint conversion table for MIST BCs",
    #                  prefix*"zeropoints.txt"))
    # datadep"MIST_zeropoints" # We always want to download these
    register(DataDep("MIST_CFHT_MegaCam", "MIST bolometric corrections for CFHT/MegaCam",
                     prefix*"CFHTugriz.txz",
                     "cf10e060b65591cb43fc55544260a75ec559c28d9454a4f681a197af26aa0fcd";
                     post_fetch_method=custom_unpack))
    register(DataDep("MIST_DECam", "MIST bolometric corrections for the Dark Energy Camera (DECam)",
                     prefix*"DECam.txz",
                     "6da654fed0e5bcf26fb196411e94c4c8b4306abc05f1dfddfab697848c3d41a7";
                     post_fetch_method=custom_unpack))
    register(DataDep("MIST_GALEX", "MIST bolometric corrections for GALEX",
                     prefix*"GALEX.txz",
                     "c390b58e3d1d15fd6c276ead3d08030ee12b36d03b7413296bb9a1caeda48766";
                     post_fetch_method=custom_unpack))
    register(DataDep("MIST_HST_ACS_HRC", "MIST bolometric corrections for HST ACS High-Resolution Channel",
                     prefix*"HST_ACSHR.txz",
                     "065e5235a55b24541a027000f358a08036e7305fb105865497e8e8c452fa512d";
                     post_fetch_method=custom_unpack))
    register(DataDep("MIST_HST_ACS_WFC", "MIST bolometric corrections for HST ACS Wide Field Channel",
                     prefix*"HST_ACSWF.txz",
                     "5e48f529d245dae122f3ac86c6524f34cf6b0952816cd4a567d78c41e3439571";
                     post_fetch_method=custom_unpack))
    register(DataDep("MIST_HST_WFC3", "MIST bolometric corrections for HST WFC3",
                     prefix*"HST_WFC3.txz",
                     "39fdf9a2606b6f66f5388c2a02bdb18bb1c32402325484c917b5cd564bcd124d";
                     post_fetch_method=custom_unpack))
    register(DataDep("MIST_HST_WFPC2", "MIST bolometric corrections for HST WFPC2",
                     prefix*"HST_WFPC2.txz",
                     "941f18684c741979d53a6b4a4083a8cd913a6fbefe03c5bfc5ee179590d14410";
                     post_fetch_method=custom_unpack))
    register(DataDep("MIST_JWST", "MIST bolometric corrections for the James Webb Space Telescope (JWST, pre-launch, added 02/15/2017)",
                     prefix*"JWST.txz",
                     "580d13c2dfe6086f8593828170f573fd653f6a51051d0437b602034cb5b49779";
                     post_fetch_method=custom_unpack))
    register(DataDep("MIST_LSST", "MIST bolometric corrections for the Vera Rubin Observatory which will conduct the Legacy Survey of Space and Time (LSST, added 02/15/2017)",
                     prefix*"LSST.txz",
                     "6546d0e444a5c363a8649f5f5515b1f569161ae79c1d21ef2525aa169c57791e";
                     post_fetch_method=custom_unpack))
    register(DataDep("MIST_PanSTARRS", "MIST bolometric corrections for PanSTARRS",
                     prefix*"PanSTARRS.txz",
                     "f841169b5dacb8063111795212b6e84a4977fbaeb5cc1e6afc9f4e3579e89a15";
                     post_fetch_method=custom_unpack))
    register(DataDep("MIST_SDSS", "MIST bolometric corrections for the Sloan Digital Sky Survey (SDSS)",
                     prefix*"SDSSugriz.txz",
                     "f2d6eb565f851ee52a45530a498817470c3567c7af245fb4228f5d65ee3f7d78";
                     post_fetch_method=custom_unpack))
    register(DataDep("MIST_SkyMapper", "MIST bolometric corrections for SkyMapper",
                     prefix*"SkyMapper.txz",
                     "7e54b341235d9e0ab824c6c77a22348f0c194a997d28faf26c910e9996b35de2";
                     post_fetch_method=custom_unpack))
    register(DataDep("MIST_HSC", "MIST bolometric corrections for Subaru Hyper Suprime-Cam (added 01/29/2019)",
                     prefix*"HSC.txz",
                     "5682c9bc306a871961b836a9a95cb7c934de254bf86678d2b6b0f9ec4fac37fe";
                     post_fetch_method=custom_unpack))
    register(DataDep("MIST_IPHAS", "MIST bolometric corrections for the INT Photometric H-α Survey (IPHAS)",
                     prefix*"IPHAS.txz",
                     "f3f9f1d9cb21ea9e6d503df0e56c528bf30c3c32a0ea82f9c96580012685be9d";
                     post_fetch_method=custom_unpack))
    register(DataDep("MIST_Spitzer", "MIST bolometric corrections for the Spitzer IRAC",
                     prefix*"SPITZER.txz",
                     "a65a25d256fe77db5400ce0acd10f70ffd80b2ece70344bbc9b6af0bea0188cb";
                     post_fetch_method=custom_unpack))
    register(DataDep("MIST_SPLUS", "MIST bolometric corrections for the S-PLUS survey (added 12/27/2020)",
                     prefix*"SPLUS.txz",
                     "9da093595b0c2c87c3b1fa7cb4dec418059f0e24520c3534a79ef586357c6cdf";
                     post_fetch_method=custom_unpack))
    register(DataDep("MIST_Swift", "MIST bolometric corrections for the Swift Observatory",
                     prefix*"Swift.txz",
                     "be5db86e0f92b2329278202703962b28dfe9134b251408bc4bea48719d2a0b40";
                     post_fetch_method=custom_unpack))
    register(DataDep("MIST_UBVRIplus", "MIST bolometric corrections for the Bessell filters, 2MASS, Kepler, Hipparcos, Gaia (DR2/MAW/EDR3), Tycho, and Tess",
                     prefix*"UBVRIplus.txz",
                     "524dc5d0fe7c9993ce9d07d32063fda8ce3dad09139739fff14868e89ae622c4";
                     post_fetch_method=custom_unpack))
    register(DataDep("MIST_UKIDSS", "MIST bolometric corrections for the UKIDSS survey",
                     prefix*"UKIDSS.txz",
                     "995464124c4b347eaaf164baf8a9a43213ae55c75dd9bd9d0b400632fa71a578";
                     post_fetch_method=custom_unpack))
    register(DataDep("MIST_UVIT", "MIST bolometric corrections for the Ultra-Violet Imaging Telescope (UVIT)",
                     prefix*"UVIT.txz",
                     "f9cba3a8636d8221c639bb894bf27c7defc66eb9e72ec33533bec1fc4480c559";
                     post_fetch_method=custom_unpack))
    register(DataDep("MIST_VISTA", "MIST bolometric corrections for the Visible and Infrared Survey Telescope for Astronomy (VISTA)",
                     prefix*"VISTA.txz",
                     "3895d64e3958739db7ba7afd45d17947094d583cd86a5aa382513cc444744d8a";
                     post_fetch_method=custom_unpack))
    register(DataDep("MIST_Washington", "MIST bolometric corrections for Washington, Strömgren, and DDO51 filters",
                     prefix*"WashDDOuvby.txz",
                     "1c0e783020c7496c31deec5313e31f0af50e7b4452da9a1062d1e4401c16925d";
                     post_fetch_method=custom_unpack))
    register(DataDep("MIST_WFIRST", "MIST bolometric corrections for the Nancy Grace Roman Telescope (formerly WFIRST)",
                     prefix*"WFIRST.txz",
                     "b8830f926858426fcfa3cc47a101911a48806daabd64817d68788a94a33f6df6";
                     post_fetch_method=custom_unpack))    
    register(DataDep("MIST_WISE", "MIST bolometric corrections for the Wide-field Infrared Explorer (WISE)",
                     prefix*"WISE.txz",
                     "c4d6dc726fca9e0d228e660a70bf4e96b9e8c2e4abde99673621a5ef3aa75dff";
                     post_fetch_method=custom_unpack))
    # Load datadeps for development
    # datadep"MIST_CFHT_MegaCam"
    # datadep"MIST_DECam"
    # datadep"MIST_GALEX"
    # datadep"MIST_HST_ACS_HRC"
    # datadep"MIST_HST_ACS_WFC"
    # datadep"MIST_HST_WFC3"
    # datadep"MIST_HST_WFPC2"
    # datadep"MIST_JWST"
    # datadep"MIST_LSST"
    # datadep"MIST_PanSTARRS"
    # datadep"MIST_SDSS"
    # datadep"MIST_SkyMapper"
    # datadep"MIST_HSC"
    # datadep"MIST_IPHAS"
    # datadep"MIST_Spitzer"
    # datadep"MIST_SPLUS"
    # datadep"MIST_Swift"
    # datadep"MIST_UBVRIplus"
    # datadep"MIST_UKIDSS"
    # datadep"MIST_UVIT"
    # datadep"MIST_VISTA"
    # datadep"MIST_Washington"
    # datadep"MIST_WFIRST"
    # datadep"MIST_WISE"
end
                     
