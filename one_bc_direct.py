import numpy as np
import ufl

from mpi4py import MPI
from petsc4py.PETSc import ScalarType
from petsc4py import PETSc

from dolfinx import mesh, fem, plot, io, la
from dolfinx.io import XDMFFile, gmshio
import dolfinx.geometry as geo
import gmsh

D_TYPE = PETSc.ScalarType

ν = 0.3
E = 2.1e10

λ = ν*E/(1+ν)/(1-2*ν)
μ = E/2/(1+ν)

ORDER = 2

dirI = 1;
dirJ = 2;


## Setting up gmsh properties
gmsh.initialize()

# Choose if Gmsh output is verbose
gmsh.option.setNumber("General.Terminal", 0)

# Set elements order to the specified one
gmsh.option.setNumber("Mesh.ElementOrder", ORDER)
# Set elements size
# gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 5) # uncomment to use for mesh refinement dependending from its surface curvature
gmsh.option.setNumber("Mesh.MeshSizeMax", 5e-2)
gmsh.option.setNumber("Mesh.MeshSizeMin", 1e-2)

# Set threads number for distrebuted meshing
# gmsh.option.setNumber("Mesh.MaxNumThreads3D", 4)

# Set mesh algorithm (default is Delaunay triangulation)
# see https://gmsh.info/doc/texinfo/gmsh.html#Choosing-the-right-unstructured-algorithm
gmsh.option.setNumber("Mesh.Algorithm3D", 3)

# gmsh.option.setNumber("Mesh.RecombinationAlgorithm",3)
# gmsh.option.setNumber("Mesh.Recombine3DAll",1)

# Set the usage of hexahedron elements 
gmsh.option.setNumber("Mesh.SubdivisionAlgorithm", 0)
## Importing RVE geometry
gmsh.open("in/MESH_STEP=0.4_POROSITY=0.637_diam=0.2.msh")

model = gmsh.model()
# model.add("main_domain")
model_name = model.getCurrent()
tags = [dimtag[1] for dimtag in model.get_entities(3)]

model.add_physical_group(dim=3, tags=tags)


# Synchronize OpenCascade representation with gmsh model
model.occ.synchronize()


# Generate the mesh
# model.mesh.generate(2)
# model.mesh.recombine()
model.mesh.generate(dim=3)

bbox = [np.Inf,
        np.Inf,
        np.Inf,
        -np.Inf,
        -np.Inf,
        -np.Inf]
for tag in tags:
    buf_bbox = model.get_bounding_box(3, tag)
    for i in range(3):
        if bbox[i] > buf_bbox[i]:
            bbox[i] = buf_bbox[i]
    for j in range(3,6):
        if bbox[j] < buf_bbox[j]:
            bbox[j] = buf_bbox[j]


# Create a DOLFINx mesh (same mesh on each rank)
msh, cell_markers, facet_markers = gmshio.model_to_mesh(model, MPI.COMM_SELF,0)
msh.name = "Box"
cell_markers.name = f"{msh.name}_cells"
facet_markers.name = f"{msh.name}_facets"

# Finalize gmsh to be able to use it again
gmsh.finalize()

def epsilon(u):
    return ufl.sym(ufl.grad(u)) # Equivalent to 0.5*(ufl.nabla_grad(u) + ufl.nabla_grad(u).T)
def sigma(u):
    return λ * ufl.nabla_div(u) * ufl.Identity(len(u)) + 2*μ*epsilon(u)


V = fem.VectorFunctionSpace(msh, ("CG", ORDER))
u = ufl.TrialFunction(V)
v = ufl.TestFunction(V)
f = fem.Constant(msh, ScalarType((0., 0., 0.)))
a = ufl.inner(sigma(u), epsilon(v)) * ufl.dx(metadata={'quadrature_degree': ORDER})
L = ufl.dot(f, v) * ufl.dx(metadata={'quadrature_degree': ORDER}) #+ ufl.dot(T, v) * ds)

eps = np.linalg.norm(np.array(bbox[0:3]) + np.array(bbox[3:]));

def left(x):
    return np.isclose(x[0], bbox[0], atol = eps);

def right(x):
    return np.isclose(x[0], bbox[3], atol = eps);

def bottom(x):
    return np.isclose(x[2], bbox[2], atol = eps);

def top(x):
    return np.isclose(x[2], bbox[5], atol = eps);

def front(x):
    return np.isclose(x[1], bbox[1], atol = eps);

def back(x):
    return np.isclose(x[1], bbox[4], atol = eps);

fdim = msh.topology.dim - 1

# find all facets on top, bottom and left boundary
left_facets = mesh.locate_entities_boundary(msh, fdim, left);
right_facets = mesh.locate_entities_boundary(msh, fdim, right);
bottom_facets = mesh.locate_entities_boundary(msh, fdim, bottom);
top_facets = mesh.locate_entities_boundary(msh, fdim, top);
front_facets = mesh.locate_entities_boundary(msh, fdim, front);
back_facets = mesh.locate_entities_boundary(msh, fdim, back);

marked_facets = np.hstack([left_facets, 
                           right_facets, 
                           bottom_facets,
                           top_facets,
                           front_facets,
                           back_facets,
                          ]);

markers = np.hstack([np.full_like(left_facets, 1),
                     np.full_like(right_facets, 2),
                     np.full_like(bottom_facets, 3),
                     np.full_like(top_facets, 4),
                     np.full_like(front_facets, 5),
                     np.full_like(back_facets, 6),
                    ]);

facets_order = np.argsort(marked_facets);

facets_tags = mesh.meshtags(msh, 
                            fdim, 
                            marked_facets[facets_order],
                            markers[facets_order]);

ds = ufl.Measure('ds', domain=msh, subdomain_data=facets_tags);

unit_disp = np.mean(np.array(bbox[3:]) - np.array(bbox[:3]));

def KUBC(x, i, j, ud):
    values = np.zeros(x.shape);
    
    values[i,:] += 0.5*ud*(x[j])/(bbox[j+3] - bbox[j]);
    values[j,:] += 0.5*ud*(x[i])/(bbox[i+3] - bbox[i]);

    return values;

# apply 2nd, 3rd and 4th constraints
facets = np.hstack([facets_tags.find(1),
                    facets_tags.find(2),
                    facets_tags.find(4),
                    facets_tags.find(5),
                    facets_tags.find(6),
                   ]);

ub_ = fem.Function(V);

full_bc = lambda x: KUBC(x, dirI, dirJ, unit_disp);

ub_.interpolate(full_bc);

nonbottom_dofs = fem.locate_dofs_topological(V,
                                         facets_tags.dim,
                                         marked_facets);
bc_ = fem.dirichletbc(ub_, nonbottom_dofs);

problem = fem.petsc.LinearProblem(a, L, bcs=[ bc_], petsc_options={"ksp_type": "preonly", "pc_type": "lu"})

uh = problem.solve()