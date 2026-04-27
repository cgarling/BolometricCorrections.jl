# Extract xz compressed tar (.txz)
function unpack_txz(fname::AbstractString, dir::AbstractString)
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

function parse_mist_v2_header(fname::AbstractString)
    for (linenum, line) in enumerate(eachline(fname))
        if linenum == 4
            # Header line: "# lgTef  logg  Fe_H a_Fe   Av   Rv  <filters...>"
            # Split on last dependent (Rv) to isolate filter names
            header = filter(!isempty, split(split(line, string(last(_mist_v2_dependents)))[2], ' '))
            return vcat(SVector(string.(_mist_v2_dependents)), convert(Vector{String}, header))
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

# Get [Fe/H] and [α/Fe] from name of MIST v2.5 BC file
# Format: feh+0.00_afe+0.0.GALEX
function mist_v2_feh_afe(fname::AbstractString)
    fname = basename(fname)
    m = match(r"feh([+-][0-9.]+)_afe([+-][0-9.]+)\.", fname)
    feh = parse(Float64, m[1])
    afe = parse(Float64, m[2])
    return (feh, afe)
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

function custom_unpack_v2(fname::AbstractString) # Path to datadep
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
    # Sort by (feh, afe) to ensure final file has sorted [Fe/H] and [α/Fe]
    feh_afe_vals = mist_v2_feh_afe.(all_files)
    idxs = sortperm(feh_afe_vals)
    header = parse_mist_v2_header(first(all_files))
    bigtable = reduce(vcat, read_mist_bc(all_files[i], header) for i in idxs)
    CSV.write(joinpath(fpath, splitext(bfname)[1]*".gz"), bigtable; compress=true)
    # Done with original .txz and extracted files so remove
    rm(fname)
    rm(out_dir; force=true, recursive=true)
end



# See https://mist.science/model_grids.html
"Available filter sets for MIST v1.2 bolometric corrections."
const BC_sets_v1 = ("CFHT_MegaCam", "DECam", "GALEX", "HST_ACS_HRC", "HST_ACS_WFC",
                    "HST_WFC3", "HST_WFPC2", "JWST", "LSST", "PanStARRS", "SDSS",
                    "SkyMapper", "Subaru/HSC", "IPHAS", "Spitzer/IRAC", "SPLUS",
                    "Swift", "UBVRIplus", "UKIDSS", "UVIT", "VISTA", "Washington",
                    "WFIRST", "WISE")
"Available filter sets for MIST v2.5 bolometric corrections."
const BC_sets_v2 = ("CFHT_MegaCam", "DECam", "Euclid", "GALEX", "HSC", "HST_ACS_HRC",
                    "HST_ACS_SBC", "HST_ACS_WFC", "HST_WFC3", "HST_WFPC2", "IPHAS",
                    "JWST", "NIRISS", "LSST", "PanSTARRS", "RoboAO", "Roman", "SDSS",
                    "Spitzer", "SPLUS", "SkyMapper", "Swift", "UBVRIplus", "UKIDSS",
                    "UVIT", "VISTA", "Washington", "WISE")
"Available filter sets for MIST bolometric corrections. Use `BC_sets.v1` or `BC_sets.v2` to access specific versions."
const BC_sets = (v1 = BC_sets_v1, v2 = BC_sets_v2)

function __init__()
    v1_prefix = "https://mist.science/BC_tables/v1/"
    # Actually just distribute zeropoints with source
    # register(DataDep("MIST_zeropoints", "Zeropoint conversion table for MIST BCs",
    #                  v1_prefix*"zeropoints.txt"))
    # datadep"MIST_zeropoints" # We always want to download these
    register(DataDep("MIST_CFHT_MegaCam", "MIST bolometric corrections for CFHT/MegaCam",
                     v1_prefix*"CFHTugriz.txz",
                     "cf10e060b65591cb43fc55544260a75ec559c28d9454a4f681a197af26aa0fcd";
                     post_fetch_method=custom_unpack))
    register(DataDep("MIST_DECam", "MIST bolometric corrections for the Dark Energy Camera (DECam)",
                     v1_prefix*"DECam.txz",
                     "6da654fed0e5bcf26fb196411e94c4c8b4306abc05f1dfddfab697848c3d41a7";
                     post_fetch_method=custom_unpack))
    register(DataDep("MIST_GALEX", "MIST bolometric corrections for GALEX",
                     v1_prefix*"GALEX.txz",
                     "c390b58e3d1d15fd6c276ead3d08030ee12b36d03b7413296bb9a1caeda48766";
                     post_fetch_method=custom_unpack))
    register(DataDep("MIST_HST_ACS_HRC", "MIST bolometric corrections for HST ACS High-Resolution Channel",
                     v1_prefix*"HST_ACSHR.txz",
                     "065e5235a55b24541a027000f358a08036e7305fb105865497e8e8c452fa512d";
                     post_fetch_method=custom_unpack))
    register(DataDep("MIST_HST_ACS_WFC", "MIST bolometric corrections for HST ACS Wide Field Channel",
                     v1_prefix*"HST_ACSWF.txz",
                     "5e48f529d245dae122f3ac86c6524f34cf6b0952816cd4a567d78c41e3439571";
                     post_fetch_method=custom_unpack))
    register(DataDep("MIST_HST_WFC3", "MIST bolometric corrections for HST WFC3",
                     v1_prefix*"HST_WFC3.txz",
                     "39fdf9a2606b6f66f5388c2a02bdb18bb1c32402325484c917b5cd564bcd124d";
                     post_fetch_method=custom_unpack))
    register(DataDep("MIST_HST_WFPC2", "MIST bolometric corrections for HST WFPC2",
                     v1_prefix*"HST_WFPC2.txz",
                     "941f18684c741979d53a6b4a4083a8cd913a6fbefe03c5bfc5ee179590d14410";
                     post_fetch_method=custom_unpack))
    register(DataDep("MIST_JWST", "MIST bolometric corrections for the James Webb Space Telescope (JWST, pre-launch, added 02/15/2017)",
                     v1_prefix*"JWST.txz",
                     "580d13c2dfe6086f8593828170f573fd653f6a51051d0437b602034cb5b49779";
                     post_fetch_method=custom_unpack))
    register(DataDep("MIST_LSST", "MIST bolometric corrections for the Vera Rubin Observatory which will conduct the Legacy Survey of Space and Time (LSST, added 02/15/2017)",
                     v1_prefix*"LSST.txz",
                     "6546d0e444a5c363a8649f5f5515b1f569161ae79c1d21ef2525aa169c57791e";
                     post_fetch_method=custom_unpack))
    register(DataDep("MIST_PanSTARRS", "MIST bolometric corrections for PanSTARRS",
                     v1_prefix*"PanSTARRS.txz",
                     "f841169b5dacb8063111795212b6e84a4977fbaeb5cc1e6afc9f4e3579e89a15";
                     post_fetch_method=custom_unpack))
    register(DataDep("MIST_SDSS", "MIST bolometric corrections for the Sloan Digital Sky Survey (SDSS)",
                     v1_prefix*"SDSSugriz.txz",
                     "f2d6eb565f851ee52a45530a498817470c3567c7af245fb4228f5d65ee3f7d78";
                     post_fetch_method=custom_unpack))
    register(DataDep("MIST_SkyMapper", "MIST bolometric corrections for SkyMapper",
                     v1_prefix*"SkyMapper.txz",
                     "7e54b341235d9e0ab824c6c77a22348f0c194a997d28faf26c910e9996b35de2";
                     post_fetch_method=custom_unpack))
    register(DataDep("MIST_HSC", "MIST bolometric corrections for Subaru Hyper Suprime-Cam (added 01/29/2019)",
                     v1_prefix*"HSC.txz",
                     "5682c9bc306a871961b836a9a95cb7c934de254bf86678d2b6b0f9ec4fac37fe";
                     post_fetch_method=custom_unpack))
    register(DataDep("MIST_IPHAS", "MIST bolometric corrections for the INT Photometric H-α Survey (IPHAS)",
                     v1_prefix*"IPHAS.txz",
                     "f3f9f1d9cb21ea9e6d503df0e56c528bf30c3c32a0ea82f9c96580012685be9d";
                     post_fetch_method=custom_unpack))
    register(DataDep("MIST_Spitzer", "MIST bolometric corrections for the Spitzer IRAC",
                     v1_prefix*"SPITZER.txz",
                     "a65a25d256fe77db5400ce0acd10f70ffd80b2ece70344bbc9b6af0bea0188cb";
                     post_fetch_method=custom_unpack))
    register(DataDep("MIST_SPLUS", "MIST bolometric corrections for the S-PLUS survey (added 12/27/2020)",
                     v1_prefix*"SPLUS.txz",
                     "9da093595b0c2c87c3b1fa7cb4dec418059f0e24520c3534a79ef586357c6cdf";
                     post_fetch_method=custom_unpack))
    register(DataDep("MIST_Swift", "MIST bolometric corrections for the Swift Observatory",
                     v1_prefix*"Swift.txz",
                     "be5db86e0f92b2329278202703962b28dfe9134b251408bc4bea48719d2a0b40";
                     post_fetch_method=custom_unpack))
    register(DataDep("MIST_UBVRIplus", "MIST bolometric corrections for the Bessell filters, 2MASS, Kepler, Hipparcos, Gaia (DR2/MAW/EDR3), Tycho, and Tess",
                     v1_prefix*"UBVRIplus.txz",
                     "524dc5d0fe7c9993ce9d07d32063fda8ce3dad09139739fff14868e89ae622c4";
                     post_fetch_method=custom_unpack))
    register(DataDep("MIST_UKIDSS", "MIST bolometric corrections for the UKIDSS survey",
                     v1_prefix*"UKIDSS.txz",
                     "995464124c4b347eaaf164baf8a9a43213ae55c75dd9bd9d0b400632fa71a578";
                     post_fetch_method=custom_unpack))
    register(DataDep("MIST_UVIT", "MIST bolometric corrections for the Ultra-Violet Imaging Telescope (UVIT)",
                     v1_prefix*"UVIT.txz",
                     "f9cba3a8636d8221c639bb894bf27c7defc66eb9e72ec33533bec1fc4480c559";
                     post_fetch_method=custom_unpack))
    register(DataDep("MIST_VISTA", "MIST bolometric corrections for the Visible and Infrared Survey Telescope for Astronomy (VISTA)",
                     v1_prefix*"VISTA.txz",
                     "3895d64e3958739db7ba7afd45d17947094d583cd86a5aa382513cc444744d8a";
                     post_fetch_method=custom_unpack))
    register(DataDep("MIST_Washington", "MIST bolometric corrections for Washington, Strömgren, and DDO51 filters",
                     v1_prefix*"WashDDOuvby.txz",
                     "1c0e783020c7496c31deec5313e31f0af50e7b4452da9a1062d1e4401c16925d";
                     post_fetch_method=custom_unpack))
    register(DataDep("MIST_WFIRST", "MIST bolometric corrections for the Nancy Grace Roman Telescope (formerly WFIRST)",
                     v1_prefix*"WFIRST.txz",
                     "b8830f926858426fcfa3cc47a101911a48806daabd64817d68788a94a33f6df6";
                     post_fetch_method=custom_unpack))    
    register(DataDep("MIST_WISE", "MIST bolometric corrections for the Wide-field Infrared Explorer (WISE)",
                     v1_prefix*"WISE.txz",
                     "c4d6dc726fca9e0d228e660a70bf4e96b9e8c2e4abde99673621a5ef3aa75dff";
                     post_fetch_method=custom_unpack))

    # V2.5 updates (added 2024-06-17)
    v2_prefix = "https://mist.science/BC_tables/v2/"
    register(DataDep("MIST_CFHT_MegaCam_v2.5", "MIST v2.5 bolometric corrections for CFHT/MegaCam",
                     v2_prefix*"CFHTugriz.txz",
                     "92facb1702564a66ee18f85da12e24b8e01d1a3ae05849105cfd3a7aacb76187";
                     post_fetch_method=custom_unpack_v2))
    register(DataDep("MIST_DECam_v2.5", "MIST v2.5 bolometric corrections for the Dark Energy Camera (DECam)",
                     v2_prefix*"DECam.txz",
                     "5682500cbdabae5bae488966534bd5e4e2740fec76571e3935c08e7f1ab64fb1";
                     post_fetch_method=custom_unpack_v2))
    register(DataDep("MIST_Euclid_v2.5", "MIST v2.5 bolometric corrections for Euclid",
                     v2_prefix*"Euclid.txz",
                     "d92a6d870ccfe6b3ff8d36d300dae63db6e1383d1e92e40f34e5ce8b1871e815";
                     post_fetch_method=custom_unpack_v2))
    register(DataDep("MIST_GALEX_v2.5", "MIST v2.5 bolometric corrections for GALEX",
                     v2_prefix*"GALEX.txz",
                     "534e7c33cb0803f1de3399fb4921b7095ba6a959e9567981768749bc2f8b6edf";
                     post_fetch_method=custom_unpack_v2))
    register(DataDep("MIST_HSC_v2.5", "MIST v2.5 bolometric corrections for Subaru Hyper Suprime-Cam",
                     v2_prefix*"HSC.txz",
                     "58f608cbec6d60272b03b076375c6f7ac29ecb9bdc11aae8cb91462412a97570";
                     post_fetch_method=custom_unpack_v2))
    register(DataDep("MIST_HST_ACS_HRC_v2.5", "MIST v2.5 bolometric corrections for HST ACS High-Resolution Channel",
                     v2_prefix*"HST_ACS_HRC.txz",
                     "8506f8561f854d7f4c9abeca9ebd087c7d347e5037bda04207a42003201e9733";
                     post_fetch_method=custom_unpack_v2))
    register(DataDep("MIST_HST_ACS_SBC_v2.5", "MIST v2.5 bolometric corrections for HST ACS Solar Blind Channel",
                     v2_prefix*"HST_ACS_SBC.txz",
                     "f07df373c0a9fdbbacdaba870e61a5b3ed3c80b5b458c27b49eebaa631007de2";
                     post_fetch_method=custom_unpack_v2))
    register(DataDep("MIST_HST_ACS_WFC_v2.5", "MIST v2.5 bolometric corrections for HST ACS Wide Field Channel",
                     v2_prefix*"HST_ACS_WFC.txz",
                     "139825efaf12c5aa0ee8f07949ab079b642c641ec54508082625c9749b8f13ea";
                     post_fetch_method=custom_unpack_v2))
    register(DataDep("MIST_HST_WFC3_v2.5", "MIST v2.5 bolometric corrections for HST WFC3",
                     v2_prefix*"HST_WFC3.txz",
                     "c82c7fb10ff41f0b9e9a7fb47c185ba5f1a8eb370e6851c7bdaee0745a50fdd8";
                     post_fetch_method=custom_unpack_v2))
    register(DataDep("MIST_HST_WFPC2_v2.5", "MIST v2.5 bolometric corrections for HST WFPC2",
                     v2_prefix*"HST_WFPC2.txz",
                     "4ab28fae743a65b51cdc6296a7e6d884c63ccfaba7b6800201773df0914fe90d";
                     post_fetch_method=custom_unpack_v2))
    register(DataDep("MIST_IPHAS_v2.5", "MIST v2.5 bolometric corrections for the INT Photometric H-α Survey (IPHAS)",
                     v2_prefix*"IPHAS.txz",
                     "8b251212b2338c885ea21cbb051980b95faab533dea1fa4ae77157e2f4d04466";
                     post_fetch_method=custom_unpack_v2))
    register(DataDep("MIST_JWST_v2.5", "MIST v2.5 bolometric corrections for the James Webb Space Telescope (JWST) NIRCam",
                     v2_prefix*"JWST.txz",
                     "dbefa0879fb2889c03d117ee877c1dd94ab68c8ed1d24b5a644a7f59248e8c74";
                     post_fetch_method=custom_unpack_v2))
    register(DataDep("MIST_NIRISS_v2.5", "MIST v2.5 bolometric corrections for the James Webb Space Telescope (JWST) NIRISS",
                     v2_prefix*"NIRISS.txz",
                     "2a4d7da891001e169f404db158646cd0b150c7e698a7a5547e745e8a9179ab40";
                     post_fetch_method=custom_unpack_v2))
    register(DataDep("MIST_LSST_v2.5", "MIST v2.5 bolometric corrections for the Vera Rubin Observatory (Rubin/LSST)",
                     v2_prefix*"LSST.txz",
                     "3b465cd545d6b3f6dad798570fd796b57a19f2f9fe37db9a13de08bec26bded8";
                     post_fetch_method=custom_unpack_v2))
    register(DataDep("MIST_PanSTARRS_v2.5", "MIST v2.5 bolometric corrections for PanSTARRS",
                     v2_prefix*"PanSTARRS.txz",
                     "025d7596df5e6956e0c32425cd03a842a615bc0640932f38083fcbd56672d15e";
                     post_fetch_method=custom_unpack_v2))
    register(DataDep("MIST_RoboAO_v2.5", "MIST v2.5 bolometric corrections for RoboAO",
                     v2_prefix*"RoboAO.txz",
                     "8b1aee60f94f7589c6d0f79bb214c1b7280278a666fb5726b61e44f7303af96a";
                     post_fetch_method=custom_unpack_v2))
    register(DataDep("MIST_Roman_v2.5", "MIST v2.5 bolometric corrections for the Nancy Grace Roman Space Telescope",
                     v2_prefix*"Roman.txz",
                     "7cac7a913017181ef5c79ae3776b6a6e77955c86a3ca1630cdc735eaf5c1e7ed";
                     post_fetch_method=custom_unpack_v2))
    register(DataDep("MIST_SDSS_v2.5", "MIST v2.5 bolometric corrections for the Sloan Digital Sky Survey (SDSS)",
                     v2_prefix*"SDSSugriz.txz",
                     "10a95f87f2775c007db587ca916830b1bdaf75728cda2897a4e2bb09803a41a0";
                     post_fetch_method=custom_unpack_v2))
    register(DataDep("MIST_Spitzer_v2.5", "MIST v2.5 bolometric corrections for the Spitzer IRAC",
                     v2_prefix*"SPITZER.txz",
                     "8b7afe1b2bb99bb927ee865e0ec6c82245f6eb652e1ef1f685046786137f954c";
                     post_fetch_method=custom_unpack_v2))
    register(DataDep("MIST_SPLUS_v2.5", "MIST v2.5 bolometric corrections for the S-PLUS survey",
                     v2_prefix*"SPLUS.txz",
                     "677c22e9e1f6df5e3256c92cd53cd5c233665f923899376665b66e216c16b6e0";
                     post_fetch_method=custom_unpack_v2))
    register(DataDep("MIST_SkyMapper_v2.5", "MIST v2.5 bolometric corrections for SkyMapper",
                     v2_prefix*"SkyMapper.txz",
                     "df1f14c7f487e4717e93ee14f8c271ff79a04323275e3c907d66e88de869cc7e";
                     post_fetch_method=custom_unpack_v2))
    register(DataDep("MIST_Swift_v2.5", "MIST v2.5 bolometric corrections for the Swift Observatory",
                     v2_prefix*"Swift.txz",
                     "572ce4c58bb5b6ed5bb80666f46b2936e647a92a0f5637b6728545d206da43b0";
                     post_fetch_method=custom_unpack_v2))
    register(DataDep("MIST_UBVRIplus_v2.5", "MIST v2.5 bolometric corrections for the Bessell filters, 2MASS, Kepler, Hipparcos, Gaia, and Tycho",
                     v2_prefix*"UBVRIplus.txz",
                     "e39a5d8871a36b450ee79928d82a9abbe104fe3771c8cb3617da60a31d911176";
                     post_fetch_method=custom_unpack_v2))
    register(DataDep("MIST_UKIDSS_v2.5", "MIST v2.5 bolometric corrections for the UKIDSS survey",
                     v2_prefix*"UKIDSS.txz",
                     "95a4ebe83d014dd40fe1128e363aa585d33e104b4ad7a315f4aeefaf028ae08e";
                     post_fetch_method=custom_unpack_v2))
    register(DataDep("MIST_UVIT_v2.5", "MIST v2.5 bolometric corrections for the Ultra-Violet Imaging Telescope (UVIT)",
                     v2_prefix*"UVIT.txz",
                     "117553e716951e697ea7707061754f88eb186ea0d6ab6cd24782574acd24cb8a";
                     post_fetch_method=custom_unpack_v2))
    register(DataDep("MIST_VISTA_v2.5", "MIST v2.5 bolometric corrections for the Visible and Infrared Survey Telescope for Astronomy (VISTA)",
                     v2_prefix*"VISTA.txz",
                     "275f408cf8be60430377ee286af84df6525668dc1c9a27f21732b5fef04bd40d";
                     post_fetch_method=custom_unpack_v2))
    register(DataDep("MIST_Washington_v2.5", "MIST v2.5 bolometric corrections for Washington, Strömgren, and DDO51 filters",
                     v2_prefix*"WashDDOuvby.txz",
                     "8c491b3cf66d91a10f679367aa15d6ea3d87ba2ce54cb19bb156fc8ccad408be";
                     post_fetch_method=custom_unpack_v2))
    register(DataDep("MIST_WISE_v2.5", "MIST v2.5 bolometric corrections for the Wide-field Infrared Explorer (WISE)",
                     v2_prefix*"WISE.txz",
                     "58779b3e7b3c8f655b89e0e36b620dae58bed05f8d82a74084b0eed0f43c6d77";
                     post_fetch_method=custom_unpack_v2))

end
