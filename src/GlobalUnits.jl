
#functions / variables for global units option. See set_units in ThreeBodyTB.jl

global global_energy_units="eV"
global global_length_units="Å"

const eV=13.6056980659
const Ang=0.529177210903

function convert_energy(energy)
    if global_energy_units=="eV"
        return energy * eV
    else
        return energy
    end
end

function convert_dos(dos)
    if global_energy_units=="eV"
        return dos / eV
    else
        return dos
    end
end

function convert_force(force)
    if global_energy_units=="eV"
        f1 = eV
    else
        f1 = 1.0
    end
    if global_length_units=="Å"
        f2 = 1/Ang
    else
        f2 = 1.0
    end
    
    return force * f1 * f2
    
end

function convert_stress(stress)
    if global_energy_units=="eV"
        f1 = eV
    else
        f1 = 1.0
    end
    if global_length_units=="Å"
        f2 = (1/Ang)^3
    else
        f2 = 1.0
    end
    
    return stress * f1 * f2
    
end

function convert_length(length)

    if global_length_units=="Å"
        f1 = Ang
    else
        f1 = 1.0
    end
    
    return length * f1
    
end

global no_display=false

function set_no_display(d)
    no_display = d
end
