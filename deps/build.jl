cd(@__DIR__)
println("building...")
println(@__DIR__)

#unpack dats
cd("../dats/pbesol/v1.2/")
run(`tar -xf els.tar`)
run(`tar -xf binary.tar`)
cd(@__DIR__)

#unpack pseudos
cd("../pseudo/gbrv_pbesol/")
run(`tar -xf all.tar.gz`)
cd(@__DIR__)


