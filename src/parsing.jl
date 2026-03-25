# Given an input NamedTuple or other type that supports `keys`
# we want to pick out the necessary values from the keys
# to pass along to the AbstractBCTable.
# Since keys are encoded in the types of NamedTuples,
# checking for the existence of the keys has no runtime cost.
# These methods should expect that the keys of the input
# e.g., iso.logL, may be vectors/arrays and broadcast appropriately.

@inline function _parse_teff(iso)
    iso_keys = keys(iso)
    if :logTe in iso_keys
        return exp10.(iso.logTe)
    elseif :logTeff in iso_keys
        return exp10.(iso.logTeff)
    elseif :Te in iso_keys
        return iso.Te
    elseif :T in iso_keys
        return iso.T
    elseif :Teff in iso_keys
        return iso.Teff
    else
        throw(ArgumentError("Provided `iso` argument does not contain a recognized effective temperature key, `(:logTe, :logTeff, :Te, :T, :Teff)`."))
    end 
end

@inline function _parse_logg(iso)
    iso_keys = keys(iso)
    if :logg in iso_keys
        return iso.logg
    elseif :log_g in iso_keys
        return iso.log_g
    else
        throw(ArgumentError("Provided `iso` argument does not contain a recognized surface gravity key, `(:logg, :log_g)`."))
    end
end

@inline function _parse_Mbol(iso)
    if :Mbol in keys(iso)
        return iso.Mbol
    else
        try # Try to use conversion from logL if available
            return Mbol.(_parse_logL(iso))
        catch
            throw(ArgumentError("Provided `iso` argument does not contain a recognized bolometric magnitude key, `(:Mbol,)`."))
        end
    end
end

@inline function _parse_logL(iso)
    iso_keys = keys(iso)
    if :logL in iso_keys
        return iso.logL
    elseif :log_L in iso_keys
        return iso.log_L
    # Calling up to _parse_Mbol would create a stack overflow if neither succeeds,
    # so we need to do an independent check here to prevent recursive loop.
    elseif :Mbol in iso_keys
        return logL.(iso.Mbol)
    else
        throw(ArgumentError("Provided `iso` argument does not contain a recognized log(L) key, `(:logL, :log_L)`."))
    end
end

@inline function _parse_Z(iso)
    iso_keys = keys(iso)
    if :Z in iso_keys
        return iso.Z
    else
        throw(ArgumentError("Provided `iso` argument does not contain a recognized metal mass fraction key, `(:Z,)`."))
    end
end

@inline function _parse_Mdot(iso)
    iso_keys = keys(iso)
    if :Mdot in iso_keys
        return iso.Mdot
    else
        throw(ArgumentError("One-argument `_parse_Mdot` failed as none of the supported keys `(:Mdot,)` found in $iso. Use two or three argument version to calculate mass-loss rate from other parameters."))
    end
end

@inline function _parse_Mdot(iso, model::AbstractMassLoss)
    iso_keys = keys(iso)
    if :Mdot in iso_keys
        return iso.Mdot
    else
        Z = _parse_Z(iso)
        return _parse_Mdot(iso, Z, model)
    end
end

@inline function _parse_Mdot(iso, Z, model::Bjorklund2021MassLoss)
    if :Mdot in keys(iso)
        return iso.Mdot
    else
        try # Try to calculate from other quantities
            logL = _parse_logL(iso)
            return model.(Z, logL)
        catch
            throw(ArgumentError("Provided `iso` argument does not contain a recognized Mdot key, `(:Mdot,)` and cannot be calculated from the keys of $iso for the mass-loss model $model."))
        end
    end
end

@inline function _parse_Mdot(iso, Z, model::Krticka2025MassLoss)
    if :Mdot in keys(iso)
        return iso.Mdot
    else
        try # Try to calculate from other quantities
            logL = _parse_logL(iso)
            Teff = _parse_teff(iso)
            return model.(Z, logL, Teff)
        catch
            throw(ArgumentError("Provided `iso` argument does not contain a recognized Mdot key, `(:Mdot,)` and cannot be calculated from the keys of $iso for the mass-loss model $model."))
        end
    end
end