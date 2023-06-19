import numpy as np;
import ufl;

from mpi4py import MPI;
from petsc4py import PETSc

from dolfinx import mesh, fem, plot, io, la;
from dolfinx.io import XDMFFile, gmshio;
import dolfinx.geometry as geo;
import gmsh;

import sys;

D_TYPE = PETSc.ScalarType

ν = 0.43
E = 1.8e9

λ = ν*E/(1+ν)/(1-2*ν)
μ = E/2/(1+ν)

def epsilon(u):
    return ufl.sym(ufl.grad(u)) # Equivalent to 0.5*(ufl.nabla_grad(u) + ufl.nabla_grad(u).T)
def sigma(u):
    return λ * ufl.nabla_div(u) * ufl.Identity(len(u)) + 2*μ*epsilon(u)

def build_nullspace(V: fem.VectorFunctionSpace):
    """Build  PETSc nullspace for 3D elasticity"""
    
    # Create vector that span the nullspace
    bs = V.dofmap.index_map_bs;
    length0 = V.dofmap.index_map.size_local;
    length1 = length0 + V.dofmap.index_map.num_ghosts;
    basis = [np.zeros(bs * length1, dtype = D_TYPE) for i in range(6)];
    
    # Get dof indices for each subspace (x, y and z dofs)
    dofs = [V.sub(i).dofmap.list.array for i in range(3)];
    
    # Set the three translational rigid body modes
    for i in range(3):
        basis[i][dofs[i]] = 1.0;
    
    # Set the three rotational rigid body modes
    x = V.tabulate_dof_coordinates();
    dofs_block = V.dofmap.list.array;
    x0, x1, x2 = x[dofs_block, 0], x[dofs_block, 1], x[dofs_block, 2];
    
    basis[3][dofs[0]] = -x1;
    basis[3][dofs[1]] = x0;
    basis[4][dofs[0]] = x2;
    basis[4][dofs[2]] = -x0;
    basis[5][dofs[2]] = x1;
    basis[5][dofs[1]] = -x2;
    
    # Create PETSc Vec objects (excluding ghosts) and normalise
    basis_petsc = [PETSc.Vec().createWithArray(x[:bs*length0], bsize=3, comm=V.mesh.comm) for x in basis]
    la.orthonormalize(basis_petsc);
    assert la.is_orthonormal(basis_petsc);
    
    #Create and return a PETSc nullspace
    return PETSc.NullSpace().create(vectors=basis_petsc);



def macro_strain(i):
    """returns the macroscopic strain for the 3 elementary load cases"""
    ϵ = np.zeros((6,), dtype = np.float64)
    ϵ[i] = 1;
    ϵ[3:] /= 2;
    return np.array([[ϵ[0], ϵ[5], ϵ[4]],
                     [ϵ[5], ϵ[1], ϵ[3]],
                     [ϵ[4], ϵ[3], ϵ[2]]]);
# ε = fem.Constant
def stress2Voigt(s):
    return ufl.as_vector([s[0,0], 
                          s[1,1], 
                          s[2,2],
                          s[1,2],
                          s[0,2],
                          s[0,1]])

def strein2Voigt(s):
    return ufl.as_vector([s[0,0], 
                          s[1,1], 
                          s[2,2],
                          2*s[1,2],
                          2*s[0,2],
                          2*s[0,1]])

def indexVoigt(i, j):
    if i==j:
        return i
    elif (i==0):
        if (j==1):
            return 5;
        else:
            return 4;
    else:
        return 3;

def main():
    


    ORDER = 2

    
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
    gmsh.open("in/"+sys.argv[1]);

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
    msh, cell_markers, facet_markers = gmshio.model_to_mesh(model, MPI.COMM_WORLD,0)
    msh.name = "Box"
    cell_markers.name = f"{msh.name}_cells"
    facet_markers.name = f"{msh.name}_facets"

    # Finalize gmsh to be able to use it again
    gmsh.finalize()
    # with io.XDMFFile(msh.comm, "out/imported_mesh.xdmf", "w") as file:
    #     file.write_mesh(msh)
    #     file.write_meshtags(cell_markers)
    #     msh.topology.create_connectivity(msh.topology.dim - 1, msh.topology.dim)
    #     file.write_meshtags(facet_markers)
        


    V = fem.VectorFunctionSpace(msh, ("CG", ORDER))
    u = ufl.TrialFunction(V)
    v = ufl.TestFunction(V)
    f = fem.Constant(msh, PETSc.ScalarType((0., 0., 0.)))
    a = fem.form(ufl.inner(sigma(u), epsilon(v)) * ufl.dx(metadata={"quadrature_degree": ORDER}))
    L = fem.form(ufl.dot(f, v) * ufl.dx(metadata={"quadrature_degree": ORDER})) #+ ufl.dot(T, v) * ds
    
    eps = np.linalg.norm(np.array(bbox[0:3]) + np.array(bbox[3:]));
    
    unit_disp =np.mean(np.array(bbox[3:]) - np.array(bbox[:3]));
    
    fdim = msh.topology.dim - 1
    
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

    def KUBC(x, i, j, ud):
        values = np.zeros(x.shape);

        values[i,:] += 0.5*ud*(x[j])/(bbox[j+3] - bbox[j]);
        values[j,:] += 0.5*ud*(x[i])/(bbox[i+3] - bbox[i]);

        return values;
    
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
    
    
    m_σ = np.zeros((6,6), dtype = np.float64);
    m_ε = np.zeros((6,6), dtype = np.float64);
    
    dx = ufl.Measure('dx', domain=msh);
    volume = fem.assemble_scalar(fem.form(fem.Constant(msh, PETSc.ScalarType(1.0)) * dx()));
    
    # set solver options
    opts = PETSc.Options();
    opts["ksp_type"] = "cg";
    opts["ksp_rtol"] = 1.0e-7;
    opts["pc_type"] = "gamg"; # geometric algebraic multigrid preconditioner

    # Use Chebyshev smothing for multigrid
    opts["mg_levels_ksp_type"] = "chebyshev";
    opts["mg_levels_pc_type"] = "jacobi";

    # Improve estimation of eigenvalues for Chebyshev smoothing
    opts["mg_levels_esteig_ksp_type"] = "cg";
    opts["mg_levels_ksp_chebyshev_esteig_steps"] = 20;

    # Create PETSc Krylov solver and turn convergence monitoring on
    solver = PETSc.KSP().create(msh.comm);
    solver.setFromOptions();
    solver.setMonitor(lambda _, its, rnorm: print(f"Iteration: {its}, rel. residual: {rnorm}"));
    ns = build_nullspace(V);
    

    for i in range(3):
        for j in range(i,3):

            ub_ = fem.Function(V);
            full_bc = lambda x: KUBC(x, i,j, unit_disp);
            ub_.interpolate(full_bc);
            nonbottom_dofs = fem.locate_dofs_topological(V,
                                                     facets_tags.dim,
                                                     marked_facets);
            bc_ = fem.dirichletbc(ub_, nonbottom_dofs);

            A = fem.petsc.assemble_matrix(a, bcs=[bc_]);
            A.assemble()
            A.setNearNullSpace(ns);
            
            b = fem.petsc.assemble_vector(L);
            fem.petsc.apply_lifting(b, [a], bcs=[[bc_]]);
            b.ghostUpdate(addv=PETSc.InsertMode.ADD, mode=PETSc.ScatterMode.REVERSE);
            fem.petsc.set_bc(b, [bc_])
            
            uh = fem.Function(V);
            # Set matrix operator
            solver.setOperators(A)
            # Set a monitor, solve linear system and display the solver
            # configuration
            
            solver.solve(b, uh.vector);
            solver.view();

            # Scatter forward the the solution ector to update ghost values
            uh.x.scatter_forward()

            for (k, case) in enumerate(["xx", "yy", "zz", "yz", "xz", "xy"]):
                ϵ_i = fem.assemble_scalar(fem.form(strein2Voigt(epsilon(uh))[k]*dx)) / volume;
                σ_i = fem.assemble_scalar(fem.form(stress2Voigt(sigma(uh))[k]*dx)) / volume;
                m_σ[indexVoigt(i, j), k] = σ_i;
                m_ε[indexVoigt(i, j), k] = ϵ_i;
                print("ε{} = {}; σ{} = {} ".format(case, ϵ_i, case, σ_i));
    
    
    np.savetxt("out/"+sys.argv[1], np.vstack((m_ε,m_σ)) , delimiter = ", ")
    # np.savetxt("out/", m_σ, delimiter = ", ");    
    
    
    
    

if __name__=="__main__":
    main();