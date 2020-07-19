from dolfin import *
import matplotlib.pyplot as plt
import time
start_time = time.time()

mesh = Mesh('figures/case_2_2d.xml')
sub_domains = MeshFunction("size_t", mesh, "figures/case_2_2d_physical_region.xml")
boundaries = MeshFunction("size_t", mesh, "figures/case_2_2d_facet_region.xml")

dx = dx(domain = mesh,subdomain_data=sub_domains)
ds = ds(domain = mesh,subdomain_data=boundaries)

#plot(sub_domains)
#plt.savefig("sub_domains.pdf")


mark = { "v_domain": 3,
         "p_domain": 2,
         "v_inflow": 9,
         "v_outflow": 10,
         "right_boundary": 5,
         "bottom_boundary" : 6,
         "left_boundary": 7,
         "top_boundary": 8,
         "interface_domain": 4,
         "wall_domain": 1,
         }

print(sub_domains.size())
#plot(sub_domains)
#plt.savefig("sub_domains.pdf")

# Define function spaces
# Define function spaces
V = VectorElement("Lagrange", mesh.ufl_cell(), 2)
Q = FiniteElement("Lagrange", mesh.ufl_cell(), 1)
TH = V * Q
W = FunctionSpace(mesh, TH)

# VISCOUS FLOW -- DIRICHLET BOUNDARIES --
noslip = Constant((0, 0))
#p0 = Constant((0.0))
#bc0_v = DirichletBC(W.sub(0), noslip, boundaries, mark["v_left_boundary"])
#bc1_v = DirichletBC(W.sub(0), noslip, boundaries, mark["v_right_boundary"])
#bc2_v = DirichletBC(W.sub(0), noslip, boundaries, mark["interface"])

# Inflow boundary condition for velocity
#inflow = Expression(("sin(pi*x[0])", "cos(pi*x[1])"), degree=2)
inflow = Expression(("0", "(6-x[0])*(7.5-x[0])"), degree=2)
#inflow = Expression(("-300", "-300"), degree=2)
bc3_v = DirichletBC(W.sub(0), inflow, boundaries, mark["v_inflow"])

#bc4_v = DirichletBC(W.sub(1), p0, boundaries, mark["v_left_boundary"])
#bc5_v = DirichletBC(W.sub(1), p0, boundaries, mark["v_right_boundary"])

# Boundary condition for pressure at inflow
# x0 = 0
#outflow = Expression(("sin(2*pi*x[0])"), degree=1)
#bc6_v = DirichletBC(W.sub(1), outflow, boundaries, mark["v_outflow"])

# POROUS FLOW -- Dirichlet BOUNDARIES --
# Inflow boundary condition for velocity
#inflow = Expression(("x[1]*(1-x[1])", "x[1]"), degree=2)
#inflow = Expression(("sin(pi*x[0])", "cos(pi*x[1])"), degree=2)
bc0_p = DirichletBC(W.sub(0), noslip, boundaries, mark["top_boundary"])

# Collect boundary conditions
bc1_p = DirichletBC(W.sub(0), noslip, boundaries, mark["right_boundary"])
bc2_p = DirichletBC(W.sub(0), noslip, boundaries, mark["left_boundary"])
#bc8 = DirichletBC(W.sub(0), noslip, boundaries, mark["v_left_boundary"])

# Boundary condition for pressure at outflow
# x0 = 0
#outflow = Expression(("-2+2*x[0]"), degree=1)
#bc3_p = DirichletBC(W.sub(1), outflow, boundaries, mark["p_bottom_boundary"])
bc3_p = DirichletBC(W.sub(0), noslip, boundaries, mark["bottom_boundary"])
p0 = Constant((0.0))
#bc6_p = DirichletBC(W.sub(1), p0, boundaries, mark["p_top1_boundary"])
#bc7_p = DirichletBC(W.sub(1), p0, boundaries, mark["p_top2_boundary"])
#bc0_i = DirichletBC(W.sub(1), p0, boundaries, mark["interface"])

bcs = [bc3_v,bc0_p, bc1_p, bc2_p, bc3_p]


'''
# viscous fluid
a_s = inner(grad(u), grad(v)) + div(v) * p + q * div(u)  # * si  # mangler sum K ledded
b_s1 = div(v) * p  # * si
b_s2 = div(u) * p  # * si
f_s = Expression(("2*mu -2 ", "0.0"), mu=4.22 * 10**(-3), degree=2)

# porous fluid
K = 1.36 * 10**(-18)
mu = 4.22 * 10**(-3)
a_d = ((mu / K) * inner(u, v) + inner(grad(p), v) + q * div(u))  # * di
b_d1 = div(v) * p  # * di  # litt feil, delvis-int
b_d2 = div(u) * p  # * di  # litt feil, delvis-int
'''

n = FacetNormal(mesh)
K_p = 1.36 * 10**(-8)
K_interface = 1.36 * 10**(-8)
K_wall = 1.36 * 10**(-18)
mu = 4.22 * 10**(-3)
#f = Expression(("pi*pi*sin(pi*x[1]) - 2*pi*cos(2*pi*x[0])", "pi*pi*cos(pi*x[0])"), degree=2)
f_p = Expression(("0", "mu/K_p*(6-x[0])*(7.5-x[0])"),mu = mu,K_p = K_p, degree=2)
f_v = Constant((0, 0))
p
#g_s = Expression(("sin(pi*x[0])", "cos(pi*x[1])"), degree=2)
#g_d = Expression(("sin(2*pi*x[0])"), degree=1)
while K_interface <= 1.36 * 10**(-1) :
    # Define variational problem
    (u, p) = TrialFunctions(W)
    (v, q) = TestFunctions(W)
    a = (inner(grad(u), grad(v)) + div(v) * p + q * div(u)) * dx(mark["v_domain"]) \
    + ((mu/ K_p) * inner(u, v) + inner(grad(p), v) + q * div(u)) * dx(mark["p_domain"]) \
    + ((mu/ K_wall) * inner(u, v) + inner(grad(p), v) + q * div(u)) * dx(mark["wall_domain"]) \
    + ((mu/ K_interface) * inner(u, v) + inner(grad(p), v) + q * div(u)) * dx(mark["interface_domain"])
    L = inner(f_p, v) * ds(mark["p_domain"]) + inner(f_v, v) * ds(mark["v_domain"])

    # Compute solution
    w = Function(W)
    #assemble_system(a, L, bcs, keep_diagonal=True)
    #assemble(a, keep_diagonal=True)
    solve(a == L, w, bcs, solver_parameters={'linear_solver':'mumps'})

    # Split the mixed solution using deepcopy
    # (needed for further computation on coefficient vector)
    (u, p) = w.split(True)

    '''
    ufile_pvd = File("outputs/parav/velocity_K_%.3e.pvd"%(K_interface))
    ufile_pvd << u
    pfile_pvd = File("outputs/parav/pressure_K_%.3e.pvd"%(K_interface))
    pfile_pvd << p
    '''

    plt.figure()
    plot(p)
    #plt.savefig("outputs/pressure_K_2_%.3e.pdf"%(K_p))
    plot(u)
    plt.title("outputs/velocity_2_K_%.3e"%(K_interface))
    plt.savefig("outputs/velocity_2_K_%.3e.png"%(K_interface), dpi=1200)
    K_interface *=10
    K_p *=10

    print("P-average (p_domain):",assemble(p*dx(mark["p_domain"])))
    print("v-average (v_domain):",assemble(p*dx(mark["v_domain"])))



print("--- %s seconds ---" % (time.time() - start_time))
