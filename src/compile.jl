#Based on code from ITensors

using PackageCompiler

default_compile_dir() = joinpath(homedir(), ".julia", "sysimages")

default_compile_filename() = "sys_threebodytb.so"

default_compile_path() = joinpath(default_compile_dir(), default_compile_filename())

function compile_note(; dir=default_compile_dir(), filename=default_compile_filename())
  path = joinpath(dir, filename)
  return """
  You will be able to start Julia with a compiled version of ThreeBodyTB using:
  ```
  ~ julia --sysimage $path
  ```
  and you should see that the startup times and JIT compilation times are substantially improved when you are using ThreeBodyTB.

  In unix, you can create an alias with the Bash command:
  ```
  ~ alias julia_threebodytb="julia --sysimage $path -e 'using ThreeBodyTB' -i"
  ```
  which you can put in your `~/.bashrc`, `~/.zshrc`, etc. This also executes `using ThreeBodyTB` so that ThreeBodyTB is loaded and ready to use, you can leave off ` -e 'using ThreeBodyTB' -i` if you don't want that. Then you can start Julia with a version of ThreeBodyTB installed with the command:
  ```
  ~ julia_threebodytb
  ```

  Note that if you update ThreeBodyTB to a new version, for example with `using Pkg; Pkg.update("ThreeBodyTB")`, you will need to run the `ThreeBodyTB.compile()` command again to recompile the new version of ThreeBodyTB.
  """
end

function compile(;
  dir::AbstractString=default_compile_dir(),
  filename::AbstractString=default_compile_filename(),
)
  println("start compile.jl")
  if !isdir(dir)
    println("""The directory "$dir" doesn't exist yet, creating it now.""")
    println()
    mkdir(dir)
  end
  path = joinpath(dir, filename)
  println(
    """Creating the system image "$path" containing the compiled version of ThreeBodyTB. This may take a few minutes.""",
  )
  create_sysimage(
    :ThreeBodyTB;
    sysimage_path=path,
    precompile_execution_file=joinpath(@__DIR__, "precompile_threebodytb.jl"),
  )
  println(compile_note(; dir=dir, filename=filename))
  return path
end

@doc """
    ThreeBodyTB.compile(; dir = "$(default_compile_dir())",
                       filename = "$(default_compile_filename())")

Compile ThreeBodyTB.jl with [PackageCompiler](https://julialang.github.io/PackageCompiler.jl/dev/). This will take some time, perhaps a few minutes.

This will create a system image containing the compiled version of ThreeBodyTB located at `dir/filename`, by default `$(default_compile_path())`.

$(compile_note())
""" compile

