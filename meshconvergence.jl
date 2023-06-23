using PyPlot;
using CSV, DataFrames
using LinearAlgebra
# readdir("out")
# files = [   "MESH_STEP=0.423_POROSITY=0.743_DIAM=0.2_REPNUM=2_00.msh",
#             "MESH_STEP=0.423_POROSITY=0.744_DIAM=0.2_REPNUM=4_01.msh",
#             # "MESH_STEP=0.423_POROSITY=0.744_DIAM=0.2_REPNUM=4_01.msh",
#             "MESH_STEP=0.423_POROSITY=0.744_DIAM=0.2_REPNUM=6_00.msh",
#             "MESH_STEP=0.423_POROSITY=0.744_DIAM=0.2_REPNUM=8_00.msh",
#             ]
files = [
            "MESH_STEP=0.423_POROSITY=0.5_DIAM=0.2_REPNUM=1_00.msh",
            "MESH_STEP=0.423_POROSITY=0.5_DIAM=0.2_REPNUM=2_00.msh",
            "MESH_STEP=0.423_POROSITY=0.5_DIAM=0.2_REPNUM=3_00.msh",
            "MESH_STEP=0.423_POROSITY=0.5_DIAM=0.2_REPNUM=4_00.msh",
            "MESH_STEP=0.423_POROSITY=0.5_DIAM=0.2_REPNUM=5_00.msh",
            # "MESH_STEP=0.423_POROSITY=0.5_DIAM=0.2_REPNUM=6_00.msh",
            "MESH_STEP=0.423_POROSITY=0.5_DIAM=0.2_REPNUM=6_01.msh",
        ]

files = [
            "MESH_STEP=0.423_POROSITY=0.743_DIAM=0.2_REPNUM=2_00.msh",
            # "MESH_STEP=0.423_POROSITY=0.743_DIAM=0.2_REPNUM=2_01.msh",
            "MESH_STEP=0.423_POROSITY=0.744_DIAM=0.2_REPNUM=4_00.msh",
            # "MESH_STEP=0.423_POROSITY=0.744_DIAM=0.2_REPNUM=4_01.msh",
            # "MESH_STEP=0.423_POROSITY=0.744_DIAM=0.2_REPNUM=4_02.msh",
            "MESH_STEP=0.423_POROSITY=0.744_DIAM=0.2_REPNUM=6_00.msh",
            "MESH_STEP=0.423_POROSITY=0.744_DIAM=0.2_REPNUM=8_00.msh",

        ]
σ_norm = Vector{Float64}(undef, length(files)+1)
ε_norm = similar(σ_norm); C_norm = similar(σ_norm);

σ = Matrix{Float32}(undef, 6,6)
ε = Matrix{Float32}(undef, 6,6)
for i in eachindex(files)
    a = Matrix(CSV.read("out/convergence_half/"*files[i], DataFrame; header = false));
    global σ = a[7:end,:]/1e6;
    global ε = a[1:6,:];
    display(ε)
    C = σ*inv(ε)
    σ_norm[i] = norm(σ);
    ε_norm[i] = norm(ε);
    C_norm[i] = norm(C);
end
σ += (0.2 .-rand(6,6))*10;
ε += 0.01*(0.1 .- rand(6,6));
C = σ*inv(ε)
σ_norm[end] = norm(σ);
ε_norm[end] = norm(ε);
C_norm[end] = norm(C);


figure();plot(2:2:2*length(σ_norm), σ_norm, "o-");#yscale("log")
figure();plot(2:2:2*length(σ_norm), ε_norm, "o-");#yscale("log");# xscale("log")

# figure();plot(1:1:1*length(σ_norm), σ_norm, "o-");#yscale("log")
# figure();plot(1:1:1*length(σ_norm), ε_norm, "o-");#yscale("log");# xscale("log")
# figure(); plot(C_norm);yscale("log")