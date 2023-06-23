using PyPlot;
using CSV, DataFrames
using LinearAlgebra
# readdir("out")
files = [   "MESH_STEP=0.423_POROSITY=0.743_DIAM=0.2_REPNUM=2_00.msh",
            "MESH_STEP=0.423_POROSITY=0.744_DIAM=0.2_REPNUM=4_01.msh",
            # "MESH_STEP=0.423_POROSITY=0.744_DIAM=0.2_REPNUM=4_01.msh",
            "MESH_STEP=0.423_POROSITY=0.744_DIAM=0.2_REPNUM=6_00.msh",
            "MESH_STEP=0.423_POROSITY=0.744_DIAM=0.2_REPNUM=8_00.msh",
            ]

σ_norm = similar(files, Float64)
ε_norm = similar(σ_norm); C_norm = similar(σ_norm);

a = Matrix(CSV.read("out/"*files[1], DataFrame; header = false));
σ = a[7:end,:]/1e6;

imshow(log.(abs.(σ)))#, cmap = "brg")


for i in eachindex(files)
    a = Matrix(CSV.read("out/"*files[i], DataFrame; header = false));
    σ = a[7:end,:]/1e6;
    ε = a[1:6,:];
    display(ε)
    C = σ*inv(ε)
    σ_norm[i] = norm(σ);
    ε_norm[i] = norm(ε);
    C_norm[i] = norm(C);
end

figure();plot(σ_norm);
figure();plot(ε_norm)
figure(); plot(C_norm)