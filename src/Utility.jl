"""
    module Utility

Some useful functions, mostly for converting stuff and loading files and reshaping stuff.
"""
module Utility

using Printf
"""
Some useful functions
"""


"""
    function reshape_vec(x, nat; strain_mode=false)

The force and relax algorithms from outside codes take in vectors, not `crystal` 's. So we have to reshape
vectors to and from `crystals`
"""    
function reshape_vec(x, nat; strain_mode=false)
#    println("RESHAPEVEC RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRr ", strain_mode)

    T=typeof(x[1])
    
#    println("size x ", size(x))
    x_r = zeros(T, nat, 3)
    for n = 1:nat
        for j = 1:3
            x_r[n,j] = x[3*(n-1) + j]
        end
    end
    
    x_r_strain = zeros(T, 3,3)

    if length(x) == 3*nat+6 && strain_mode
        x_r_strain[1,1] = x[3*nat+1]
        x_r_strain[2,2] = x[3*nat+2]
        x_r_strain[3,3] = x[3*nat+3]

        x_r_strain[2,3] = 0.5*x[3*nat+4]
        x_r_strain[3,2] = 0.5*x[3*nat+4]

        x_r_strain[1,3] = 0.5*x[3*nat+5]
        x_r_strain[3,1] = 0.5*x[3*nat+5]

        x_r_strain[1,2] = 0.5*x[3*nat+6]
        x_r_strain[2,1] = 0.5*x[3*nat+6]
    elseif length(x) == 3*nat+6 
        x_r_strain[1,1] = x[3*nat+1]
        x_r_strain[2,2] = x[3*nat+2]
        x_r_strain[3,3] = x[3*nat+3]

        x_r_strain[2,3] = x[3*nat+4]
        x_r_strain[3,2] = x[3*nat+4]

        x_r_strain[1,3] = x[3*nat+5]
        x_r_strain[3,1] = x[3*nat+5]

        x_r_strain[1,2] = x[3*nat+6]
        x_r_strain[2,1] = x[3*nat+6]

    elseif length(x) == 3*nat+9
        x_r_strain[1,1] = x[3*nat+1]
        x_r_strain[1,2] = x[3*nat+2]
        x_r_strain[1,3] = x[3*nat+3]
        x_r_strain[2,1] = x[3*nat+4]
        x_r_strain[2,2] = x[3*nat+5]
        x_r_strain[2,3] = x[3*nat+6]
        x_r_strain[3,1] = x[3*nat+7]
        x_r_strain[3,2] = x[3*nat+8]
        x_r_strain[3,3] = x[3*nat+9]
    else
        println("I'm confusing about the length reshape_vec $nat ", length(x) )
    end

    return x_r, x_r_strain
end

"""
    function inv_reshape_vec(x, strain, nat; strain_mode=true)

The force and relax algorithms from outside codes take in vectors, not `crystal` 's. So we have to reshape
vectors to and from `crystals`
"""
function inv_reshape_vec(x, strain, nat; strain_mode=true)
    T=typeof(x[1])
    if strain_mode
        x_r = zeros(T, nat*3 + 6)
    else
        x_r = zeros(T, nat*3 + 9)
    end

    for n = 1:nat
        for j = 1:3
            x_r[3*(n-1) + j] = x[n,j] 
        end
    end
    if strain_mode
        x_r[3*nat + 1] = strain[1,1]
        x_r[3*nat + 2] = strain[2,2]
        x_r[3*nat + 3] = strain[3,3]
        x_r[3*nat + 4] = (strain[2,3]+strain[3,2])
        x_r[3*nat + 5] = (strain[1,3]+strain[3,1])
        x_r[3*nat + 6] = (strain[1,2]+strain[2,1])
    else
        x_r[3*nat+1] = strain[1,1]  
        x_r[3*nat+2] = strain[1,2]  
        x_r[3*nat+3] = strain[1,3]  
        x_r[3*nat+4] = strain[2,1]  
        x_r[3*nat+5] = strain[2,2]  
        x_r[3*nat+6] = strain[2,3]  
        x_r[3*nat+7] = strain[3,1]  
        x_r[3*nat+8] = strain[3,2]  
        x_r[3*nat+9] = strain[3,3]  
    end

    return x_r
end


function cutoff_fn(num, min_c, max_c)

    if num < 1e-4
        return 0.0
    end
    if num < min_c
        return 1.0
    elseif num > max_c
        return 0.0
    else
        t = (num - min_c)/(max_c-min_c)
        return 1.0 - 10.0 * t^3 + 15.0 *  t^4  - 6.0 * t^5
    end
end

function cutoff_fn_fast(num, min_c, max_c)

    t = (num - min_c)/(max_c-min_c)
    return max.(min.(1.0 - 10.0 * t^3 + 15.0 *  t^4  - 6.0 * t^5, 1.0), 0.0)
end

function arr2str(a::Array{Int,1})
    st=""
    for i = 1:size(a,1)
        t= @sprintf(" % 7s", a[i])
        st=st*t
    end
    return st

end

function arr2str(a::Array{Float64,1})
    st=""
    for i = 1:size(a,1)
        t= @sprintf(" % 2.16e", a[i])
        st=st*t
    end
    return st

end

#now matches all number arrays!!!
function arr2str(a::AbstractArray{<:Number,2})
    st=""
    for i = 1:size(a,1)
        t= arr2str(a[i,:])
        st=st*t*"\n"
    end
    return st

end

#now matches all number arrays!!!
function arr2str(a::AbstractArray{Any,2})
    st=""
    for i = 1:size(a,1)
        t = str_w_spaces(a[i,:])
        st=st*t*"\n"
    end
    return st

end


function str_w_spaces(a)
    st=""
    for i = 1:(size(a,1)-1)
        st= st*string(a[i])* "  "
    end
    st= st*string(a[size(a,1)])

    return st

end

function str_w_spaces(a::Set)
    st=""
    for s in a
        st= st*string(s)* "  "
    end
#    st= st*string(a[size(a,1)])

    return st

end


function parse_str_ARR_float(sp)
   if typeof(sp) == String
        return map(x->parse(Float64,x),split(sp))
    else
        return map(x->parse(Float64,x),sp)
    end
end


function parse_str_ARR_complex(sp)

    
    t = map(x->parse(Float64,x),split(replace(sp, "," => " ")))
    return t[1:2:end] + im*t[2:2:end]

end

function dict2str(a::Dict)
    st=""
    for k in keys(a)
        st *= string(k) *" => "* string(a[k]) *"\n"
    end
    return st
end

function str2tuplesdict(st::String, d=missing)
    if ismissing(d)
        d = Dict()
    end
    
    for line in split(st, "\n")
        sp = split(line, " => ")
        if length(sp) == 2
            key = eval(Meta.parse(sp[1]))
            val = eval(Meta.parse(sp[2]))        
            d[key] = val
        end
    end
    return d
end

function write_to_file(str, filename, directory="./")

    f=open("$directory/$filename", "w")
    write(f, str)
    close(f)

end



end
