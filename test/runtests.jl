using Opcode
using Base.Test

using PlyIO

function imwritepnm(x::Array{T,2},fname) where {T <: Real}
    width, height = size(x)
    xmax = maximum(x)
    open(fname,"w") do stream
        write(stream, "P6\n$width $height\n255\n")
        for j in 1:height, i in 1:width, k = 1:3
            write(stream,round(UInt8,255*x[i,j]/xmax))
        end
    end
end

function readply(filename)
    mesh = load_ply(filename)
    verts = vcat(Float64.(mesh["vertex"]["x"]).', Float64.(mesh["vertex"]["y"]).', Float64.(mesh["vertex"]["z"]).')
    faces = hcat(collect(Int.(y).+1 for y in mesh["face"]["vertex_indices"])...)
    (verts,faces)
end

verts,faces = readply("chopper.ply")
tree = AABBTree(verts,faces)

function imgit(tree, verts)
    imgsz = 512;
    bmin = minimum(verts[[2,3],:],2)
    bmax = maximum(verts[[2,3],:],2)
    mindist = 1e9;
    maxdist = -1e9;

    y = 1:imgsz;
    x = 1:imgsz;

    # view from the side
    fracy = (imgsz-y)./imgsz;
    cy = (bmax[1] - bmin[1]) .* fracy .+ bmin[1];
    fracz = x./imgsz;
    cz = (bmax[2] - bmin[2]) .* fracz .+ bmin[2];

    fy = repmat(cy.', length(cz), 1)
    fz = repmat(cz, 1, length(cy))
    #[fz,fy] = meshgrid(cy,cz); christ no meshgrid??

    from = vcat(-1000*ones(size(fy[:]))', fy[:]', fz[:]');
    to = vcat(1000*ones(size(fy[:]))', fy[:]', fz[:]');

    hits,dists,faceinds,barys = raycast(tree, from, to, mode=1)
    img = reshape(dists,size(fz))
end

img = imgit(tree,verts)
imwritepnm(img, "x1.ppm")

verts .= verts .+ rand(size(verts)).*20
deform!(tree, verts)

img = imgit(tree,verts)
imwritepnm(clamp.(img,0,1), "x2.ppm")



if false
verts = rand(3,20); faces = hcat(randperm(size(verts,2)), randperm(size(verts,2)),randperm(size(verts,2)))';
a = Opcode.aabb_tree_create(verts,size(verts,2),faces,size(faces,2));
origins = rand(3,10); dirs = rand(3,10); hits = Vector{Int}(); dists = Vector{Float64}(); faceinds = Vector{Int}(); barys = Vector{Float64}(0);
Opcode.aabb_tree_intersect(a, origins, dirs, hits, dists, faceinds, barys, size(origins,2))
end
