using Optim
using LinearAlgebra
using CSV, DataFrames

using PyPlot

files025 = [
    "MESH_STEP=0.05_POROSITY=0.7_DIAM=0.0182_FILLNESS=1.0_00.msh",
    "MESH_STEP=0.05_POROSITY=0.7_DIAM=0.0199_FILLNESS=0.834_00.msh",
    "MESH_STEP=0.05_POROSITY=0.7_DIAM=0.02221_FILLNESS=0.666_00.msh",
    "MESH_STEP=0.05_POROSITY=0.7_DIAM=0.0256_FILLNESS=0.5_00.msh",
    "MESH_STEP=0.05_POROSITY=0.8_DIAM=0.0148_FILLNESS=1.0_00.msh",
    "MESH_STEP=0.05_POROSITY=0.8_DIAM=0.0162_FILLNESS=0.834_00.msh",
    "MESH_STEP=0.05_POROSITY=0.8_DIAM=0.0181_FILLNESS=0.666_00.msh",
    "MESH_STEP=0.05_POROSITY=0.8_DIAM=0.02085_FILLNESS=0.5_00.msh",
    "MESH_STEP=0.05_POROSITY=0.8_DIAM=0.0255_FILLNESS=0.333_00.msh",
    "MESH_STEP=0.05_POROSITY=0.9_DIAM=0.0114_FILLNESS=0.834_00.msh",
    "MESH_STEP=0.05_POROSITY=0.9_DIAM=0.01275_FILLNESS=0.666_00.msh",
    "MESH_STEP=0.05_POROSITY=0.9_DIAM=0.0147_FILLNESS=0.5_00.msh",
    "MESH_STEP=0.05_POROSITY=0.9_DIAM=0.018_FILLNESS=0.333_00.msh",
]
files040 = [
    "MESH_STEP=0.05_POROSITY=0.6_DIAM=0.02183_FILLNESS=1.0_00.msh",
    "MESH_STEP=0.05_POROSITY=0.6_DIAM=0.02378_FILLNESS=0.833_00.msh",
    "MESH_STEP=0.05_POROSITY=0.6_DIAM=0.02645_FILLNESS=0.667_00.msh",
    "MESH_STEP=0.05_POROSITY=0.6_DIAM=0.0304_FILLNESS=0.5_00.msh",
    "MESH_STEP=0.05_POROSITY=0.7_DIAM=0.01875_FILLNESS=1.0_00.msh",
    "MESH_STEP=0.05_POROSITY=0.7_DIAM=0.02045_FILLNESS=0.833_00.msh",
    "MESH_STEP=0.05_POROSITY=0.7_DIAM=0.02278_FILLNESS=0.667_00.msh",
    "MESH_STEP=0.05_POROSITY=0.7_DIAM=0.0262_FILLNESS=0.5_00.msh",
    "MESH_STEP=0.05_POROSITY=0.8_DIAM=0.01515_FILLNESS=1.0_00.msh",
    "MESH_STEP=0.05_POROSITY=0.8_DIAM=0.01515_FILLNESS=1.0_00.msh",
    "MESH_STEP=0.05_POROSITY=0.8_DIAM=0.01655_FILLNESS=0.833_00.msh",
    "MESH_STEP=0.05_POROSITY=0.8_DIAM=0.01845_FILLNESS=0.667_00.msh", 
    "MESH_STEP=0.05_POROSITY=0.8_DIAM=0.0212_FILLNESS=0.5_00.msh",
    "MESH_STEP=0.05_POROSITY=0.8_DIAM=0.0259_FILLNESS=0.333_00.msh", 
    "MESH_STEP=0.05_POROSITY=0.9_DIAM=0.01155_FILLNESS=0.833_00.msh",
    "MESH_STEP=0.05_POROSITY=0.9_DIAM=0.0129_FILLNESS=0.667_00.msh",
    "MESH_STEP=0.05_POROSITY=0.9_DIAM=0.01485_FILLNESS=0.5_00.msh",
    "MESH_STEP=0.05_POROSITY=0.9_DIAM=0.01815_FILLNESS=0.333_00.msh",
]

path = "out/meshes/";

porosities = [
    0.9, 
    0.8, 
    0.7,
    0.6,
    ];

dirs = [
    [
        [1,1],
        [2,2],
        [3,3],
    ],
    [
        [1,2],
        [1,3],
        [2,3],
        [2,1],
        [3,1],
        [3,2],
    ],
    [
        [4,4],
        [5,5],
        [6,6],
    ],
    [
        [4,5],
        [4,6],
        [5,6],
        [5,4],
        [6,4],
        [6,5],
    ],
    [
        [1,4],
        [2,5],
        [3,6],
        [4,1],
        [5,2],
        [6,3],
    ],
    [
        [1,5],
        [1,6],
        [2,6],
        [2,4],
        [3,4],
        [3,5],
        [5,1],
        [6,1],
        [6,2],
        [4,2],
        [4,3],
        [5,3],
    ]
]
rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams");
rcParams["font.family"] = "Times New Roman";
rcParams["text.usetex"] = true;
# "serif":["Times"]})
k = 6;
WIDTH_SIZE = 4
HEIGHT_SIZE = 3
linestyles = ["solid", "dashed", "dashdot", "dotted"];
figure(figsize=(WIDTH_SIZE,HEIGHT_SIZE))
iline = 0;
for prsty in porosities
    val = [];
    diam = [];
    for file in files040
        rex = r"(?<=POROSITY=)[0-9]*\.?[0-9]*"
        mch = match(rex, file);
        if parse(Float64, mch.match) == prsty
            a = Matrix(CSV.read(path*file, DataFrame; header = false));
            σ = a[7:end,:]./1e6;
            push!(val, []);
            for (i,inds) in enumerate(dirs[k])
                push!(val[end], abs(σ[dirs[k][i][1], dirs[k][i][2]]));
            end
            
            push!(diam, 1000*parse(Float64, match(r"(?<=DIAM=)[0-9]*\.?[0-9]*", file).match));

        end
    end
    iline +=1;
    Cij_mean = [mean(v) for v in val]; Cij_max = [maximum(v) for v in val]; Cij_min = [minimum(v) for v in val];
    plot(diam, Cij_mean, linestyle=linestyles[iline], label="$prsty")
    fill_between(diam, Cij_min, Cij_max, alpha=0.2, label="_nolegend_")
    # errorbar(diam, Cij_mean, [Cij_min, Cij_max], fmt="-",label="$prsty")

    # plot(diam, val, marker="o", label="$prsty")
end

legend(["$prsty" for prsty in porosities], loc="upper left")
grid("on");
xlim([10,31]);
xticks(10:5:31)

yt_fisrt = 0;
yt_last = 35;
yt_step = 10;
ylim([yt_fisrt, yt_last]);
yticks(yt_fisrt:yt_step:yt_last)

gca().minorticks_on()
ylabel("\$C_$k\$, MPa");
xlabel("Fiber diameter \$d_f\$, \$\\mathrm{\\mu}\$m");
gca().get_legend().set_title("Porosity \$ \\phi \$:");
tight_layout();
savefig("C_$(k)_040.pdf")
x = 0:0.1:10