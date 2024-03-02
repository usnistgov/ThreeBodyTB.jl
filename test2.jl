### A Pluto.jl notebook ###
# v0.18.0

using Markdown
using InteractiveUtils

# ╔═╡ f3aae5f6-893a-11ec-33a8-0961474d4592
using LinearAlgebra

# ╔═╡ cdf555aa-3e3a-4863-810a-b3a6922637c6
H = [1.0 -0.1; -0.1 1.0]

# ╔═╡ ce92365b-142f-46e5-9571-b6e383930385
vals, vects = eigen(H)

# ╔═╡ f944d037-48a7-4597-b442-358ef5ef1dee
md"bonding state energy:  $(round(vals[1], digits=3))"

# ╔═╡ 6758a7d8-82e7-4d57-add6-915bd634c747
vects[:,1]

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.2"
manifest_format = "2.0"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
"""

# ╔═╡ Cell order:
# ╠═f3aae5f6-893a-11ec-33a8-0961474d4592
# ╠═cdf555aa-3e3a-4863-810a-b3a6922637c6
# ╠═ce92365b-142f-46e5-9571-b6e383930385
# ╠═f944d037-48a7-4597-b442-358ef5ef1dee
# ╠═6758a7d8-82e7-4d57-add6-915bd634c747
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
