### A Pluto.jl notebook ###
# v0.14.8

using Markdown
using InteractiveUtils

# ╔═╡ c8abedfc-102d-4cb9-b8da-61acf49fcb77
using StatsBase, LinearAlgebra, AbstractAlgebra

# ╔═╡ eae23ae4-dc0d-11eb-0774-ed8da89022a4
begin
	n = 10
	k = 2
end

# ╔═╡ 61a9a089-7f88-427a-a894-ca507f547910
M(m) = matrix(GF(2), m)

# ╔═╡ 9090d304-3bea-4fe8-abd9-f597420b1b3a
begin
	H = nothing
	G = nothing
	while true
		H = zeros(Int, n-k, n)
		for r in 1:n-k
			H[r,sample(1:n, 5, replace=false)] .= 1
		end
		H = M(H)
		G = transpose(nullspace(H)[2])
		if size(G)[1] == 2
			break
		end
	end
end

# ╔═╡ d69579c1-3ec6-4591-9220-8f6a821d9bd1
Matrix(H)

# ╔═╡ 49819c37-a198-45d4-8b58-652f38088495
Matrix(G)

# ╔═╡ 0550cd83-ee09-416c-a276-362f3c31e8ab
syndrome(c) = H * c'

# ╔═╡ 743ee253-9c16-404f-8205-0b887d650c5c
function decode(c)
	while true
		s = sum(getfield.(syndrome(c), :d))
		@info "Starting weight $(s)"
		if s == 0
			break
		end
		best = nothing
		for i in 1:length(c)
			c[1,i] = 1-c[i]
			s1 = sum(getfield.(syndrome(c), :d))
			if s1 < s
				best = i
				s = s1
			end
			c[1,i] = 1-c[i]
		end
		@info "Improving syndrome weight to $(s) with $(best)"
		if best == nothing
			break
		end
		c[1,best] = 1-c[best]
	end
	c
end

# ╔═╡ b94d33c1-4bd3-4726-971b-cea97831401a
begin
	c = M([0 1]) * G
	w = c + M([0 0 1 0 0 0 0 0 0 1])
	cp = decode(copy(w))
	c, w, sum(syndrome(w)), cp
end

# ╔═╡ 889e6034-a09e-447b-9e6a-43446a9fa4e3
getfield.(solve(G', cp'), :d)

# ╔═╡ Cell order:
# ╠═c8abedfc-102d-4cb9-b8da-61acf49fcb77
# ╠═eae23ae4-dc0d-11eb-0774-ed8da89022a4
# ╠═9090d304-3bea-4fe8-abd9-f597420b1b3a
# ╠═d69579c1-3ec6-4591-9220-8f6a821d9bd1
# ╠═49819c37-a198-45d4-8b58-652f38088495
# ╠═0550cd83-ee09-416c-a276-362f3c31e8ab
# ╠═61a9a089-7f88-427a-a894-ca507f547910
# ╠═743ee253-9c16-404f-8205-0b887d650c5c
# ╠═b94d33c1-4bd3-4726-971b-cea97831401a
# ╠═889e6034-a09e-447b-9e6a-43446a9fa4e3
