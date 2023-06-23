using PyPlot;
using CSV, DataFrames
using LinearAlgebra

rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams");
rcParams["font.family"] = "Times New Roman";
rcParams["text.usetex"] = true;

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
files = [
        "MESHSIZE=0.977e-3_MSFC=3_POROSITY=0.639_DIAM=0.02_00.msh",
        "MESHSIZE=1.381e-3_MSFC=3_POROSITY=0.639_DIAM=0.02_00.msh",
        "MESHSIZE=1.953e-3_MSFC=3_POROSITY=0.639_DIAM=0.02_00.msh",
        "MESHSIZE=2.762e-3_MSFC=3_POROSITY=0.639_DIAM=0.02_00.msh",
        "MESHSIZE=3.906e-3_MSFC=3_POROSITY=0.639_DIAM=0.02_00.msh",
        "MESHSIZE=5.524e-3_MSFC=3_POROSITY=0.639_DIAM=0.02_00.msh",
        "MESHSIZE=7.812e-3_MSFC=3_POROSITY=0.639_DIAM=0.02_00.msh",
        "MESHSIZE=11.049e-3_MSFC=3_POROSITY=0.639_DIAM=0.02_00.msh",
        "MESHSIZE=15.625e-3_MSFC=3_POROSITY=0.639_DIAM=0.02_00.msh",
        "MESHSIZE=22.097e-3_MSFC=3_POROSITY=0.639_DIAM=0.02_00.msh",
        "MESHSIZE=31.25e-3_MSFC=3_POROSITY=0.639_DIAM=0.02_00.msh",
        "MESHSIZE=44.194e-3_MSFC=3_POROSITY=0.639_DIAM=0.02_00.msh",
        "MESHSIZE=62.5e-3_MSFC=3_POROSITY=0.639_DIAM=0.02_00.msh",
        "MESHSIZE=88.388e-3_MSFC=3_POROSITY=0.639_DIAM=0.02_00.msh",
        "MESHSIZE=125.0e-3_MSFC=3_POROSITY=0.639_DIAM=0.02_00.msh",
]
path = "out/Mesh_convergancy/"
σ_err = Vector{Float64}(undef, length(files)-1)
msh_size = similar(σ_err);
# ε_norm = similar(σ_norm); C_norm = similar(σ_norm);

σ = Matrix{Float32}(undef, 6,6)
# ε = Matrix{Float32}(undef, 6,6)
a = Matrix(CSV.read(path*files[1], DataFrame; header = false));
σ0 = a[7:end,:]
rex = r"(?<=MESHSIZE=)[0-9]*\.?[0-9]*e-3"
mch = match(rex, files[1]);
parse(Float64, mch.match);

println("==============================");
for i in eachindex(files[2:end])
    a = Matrix(CSV.read(path*files[1+i], DataFrame; header = false));
    global σ = a[7:end,:];
    rex = r"(?<=MESHSIZE=)[0-9]*\.?[0-9]*e-3"
    mch = match(rex, files[1+i]);
    # global ε = a[1:6,:];
    display(σ0 - σ)
    # C = σ*inv(ε)
    σ_err[i] = norm(σ0 - σ)/norm(σ0);
    msh_size[i] = parse(Float64, mch.match)
    # ε_norm[i] = norm(ε);
    # C_norm[i] = norm(C);
end
# σ += (0.2 .-rand(6,6))*10;
# ε += 0.01*(0.1 .- rand(6,6));
# C = σ*inv(ε)
# σ_norm[end] = norm(σ);
# ε_norm[end] = norm(ε);
# C_norm[end] = norm(C);


figure();plot(msh_size, σ_err, "k.-"); yscale("log"); xscale("log")
xlabel("Element size, mm"); 
ylabel("\${L_2(\\hat{C}-C)}/{L_2(\\hat{C})}\$");
grid("minor"); grid("major")
# figure();plot(2:2:2*length(σ_norm), ε_norm, "o-");#yscale("log");# xscale("log")

# figure();plot(1:1:1*length(σ_norm), σ_norm, "o-");#yscale("log")
# figure();plot(1:1:1*length(σ_norm), ε_norm, "o-");#yscale("log");# xscale("log")
# figure(); plot(C_norm);yscale("log")
