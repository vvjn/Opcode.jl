module Opcode

export AABBTree, intersect, intersect!, deform!, raycast, raycast!

using CxxWrap

const depsfile = joinpath(dirname(dirname(@__FILE__)), "deps", "deps.jl")
if !isfile(depsfile)
    error("$depsfile not found, Opcode did not build properly")
end
include(depsfile)
const jlcxx_path = _l_jlcxx

# Wrap the functions defined in C++
wrap_module(Opcode._l_opcode, Opcode)

const AABBTree = Opcode.ModelRef

function AABBTree(verts::Array{<:Real,2}, faces::Array{<:Integer,2})
    size(verts,1)==3 && size(faces,1)==3 || error("require verts : 3 x n and faces : 3 x m")
    AABBTree(Opcode.aabb_tree_create(verts,size(verts,2),faces,size(faces,2)))
end

# cast first n rays in (origins,dirs)
function raycast(tree::AABBTree, origins::Array{<:Real,2}, dirs::Array{<:Real,2}, n::Integer = size(origins,2);
                 mode::Int=0)
    size(origins,1)==size(dirs,1)==3 &&
    n <= min(size(origins,2), size(dirs,2)) || error("require origins : 3 x n and dirs : 3 x n")
    if mode==0
        hits = Vector{Int}()
    	dists = Vector{Float64}(); 
    	faceinds = Vector{Int}();  
    	barys = Vector{Float64}();
    	Opcode.aabb_tree_intersect(tree, origins, dirs, hits, dists, faceinds, barys, n)
    	(hits,dists,faceinds,reshape(barys,2,length(hits)))
    else
        hits = Vector{Int}(n)
        dists = Vector{Float64}(n);
        faceinds = Vector{Int}(n);
        barys = Matrix{Float64}(2,n);
        Opcode.aabb_tree_intersect_mut(tree, origins, dirs, hits, dists, faceinds, barys, n)
        (hits,dists,faceinds,reshape(barys,2,length(hits)))
    end
end

function raycast!(tree::AABBTree, origins::Array{<:Real,2}, dirs::Array{<:Real,2},
                    hits::Vector{<:Integer}, dists::Vector{<:Real}, faceinds::Vector{<:Integer}, 
		    barys::Array{<:Real,2}, n :: Integer = size(origins,2))
    size(origins,1)==size(dirs,1)==3 && size(barys,1)==2 &&
    n <= min(size(origins,2),size(dirs,2),length(hits),length(dists),length(faceinds),size(barys,2)) ||
    error("input dimensions")
    Opcode.aabb_tree_intersect_mut(tree, origins, dirs, hits, dists, faceinds, barys, n)
    tree
end

function deform!(tree::AABBTree, verts::Array{<:Real,2})
    size(verts,1)==3 || error("require verts : 3 x n")
    Opcode.aabb_tree_deform(tree, verts, size(verts,2))
    tree
end

function bary2eucl!(x::Array{<:Real,2}, bary::Array{<:Real,2}, faceinds::Vector{<:Integer},
                   verts::Array{<:Real,2}, faces::Array{<:Integer})	 	   
    n = size(bary,2)
    size(bary,1)==2 || error("bary : 2 x n")
    @simd for i = 1:n
    	b2 = bary[1,i]
	b3 = bary[2,i]
        b1 = 1 - b2 - b3
        v1,v2,v3 = (faces[1,faceinds[i]], faces[2,faceinds[i]], faces[3,faceinds[i]])
        v1x,v1y,v1z = (verts[1,v1], verts[2,v1], verts[3,v1])
        v2x,v2y,v2z = (verts[1,v2], verts[2,v2], verts[3,v2])
        v3x,v3y,v3z = (verts[1,v3], verts[2,v3], verts[3,v3])
        x[1,i] = b1*v1x + b2*v2x + b3*v3x
        x[2,i] = b1*v1y + b2*v2y + b3*v3y
        x[3,i] = b1*v1z + b2*v2z + b3*v3z
    end
    x
end
function bary2eucl!(bary::Array{T,2}, faceinds::Vector{<:Integer},
                    verts::Array{<:Real,2}, faces::Array{<:Integer}) where {T<:Real}
    x = Array{T,2}(3,n)
    bary2eucl(x,bary,faceinds,verts,faces)
end

end # module
