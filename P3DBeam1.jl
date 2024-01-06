include("Elaslin.jl")

# Empty structure
EST = Elaslin.EmptyEstructura()

# Optimization parameters
EST.Problema = "P3DBeam1"   # Name of the problem
EST.beta = 3.0              # Reliability index
EST.vtol = 0.0001           # Volume tolerance for drawings
EST.Vol = 1.0               # Maximum volume
EST.TolGap = 0.0            # Gap tolerance in optimization (zero if not used)
EST.dim = 3                 # Space dimension

# Mesh grid
xi = 0.0; xf = 1.0
yi = 0.0; yf = 4.0
zi = 0.0; zf = 1.0

nx = 2
ny = 8
nz = 2
lv = 100 # Connectivity level

# Mesh generation
Elaslin.Malla3D!(EST, xi, xf, yi, yf, zi, zf, nx, ny, nz, lv)

# Material properties
EST.PropMat = [100.0] # Young modulus

# Support and load matrices
# EST.McndC = [Entry, Node, Load state, TypeX, TypeY, TypeZ]
# EST.McndV = [ValX, ValY, ValZ]
# Type = 0: Force in the node
# Type = 1: Support in the node

# Supports at nodes with coordinate Y=0
Aux = [1:EST.Nnodo...]
Aux = Aux[EST.Mnodo[Aux,2].==0.0,:]

NN = length(Aux)
EST.Ncond = 0
EST.Ncond += NN
EST.McndC = vcat(EST.McndC, hcat([1:NN...], Aux, zeros(Int,NN,1), ones(Int,NN,3)))
EST.McndV = vcat(EST.McndV, zeros(NN,3))

# Vertical force 1
Aux = [1:EST.Nnodo...]
Aux = Aux[EST.Mnodo[Aux,1].==0.5]
Aux = Aux[EST.Mnodo[Aux,2].==4.0]
Aux = Aux[EST.Mnodo[Aux,3].==0.0]

EST.Ncond += 1
EST.McndC = vcat(EST.McndC, [EST.Ncond Aux[1] 1 0 0 0])
EST.McndV = vcat(EST.McndV, [0.0  0.0 -1.0])

# Vertical force 2
EST.Ncond += 1
EST.McndC = vcat(EST.McndC, [EST.Ncond Aux[1] 2 0 0 0])
EST.McndV = vcat(EST.McndV, [0.3  0.0  0.0])

# Vertical force 3
Aux = [1:EST.Nnodo...]
Aux = Aux[EST.Mnodo[Aux,1].==0.5]
Aux = Aux[EST.Mnodo[Aux,2].==4.0]
Aux = Aux[EST.Mnodo[Aux,3].==1.0]

EST.Ncond += 1
EST.McndC = vcat(EST.McndC, [EST.Ncond Aux[1] 3 0 0 0])
EST.McndV = vcat(EST.McndV, [0.0  0.0 -1.0])

# Vertical force 4
EST.Ncond += 1
EST.McndC = vcat(EST.McndC, [EST.Ncond Aux[1] 4 0 0 0])
EST.McndV = vcat(EST.McndV, [0.3  0.0  0.0])

# Matrix R
EST.RCarg = [
     0.0 0.0 0.0 0.0
     0.0 0.0 0.0 0.0
     0.0 0.0 0.0 0.0
     0.0 0.0 0.0 0.0
    ]

# Optimization routines
Elaslin.Elaslin!(EST)
Elaslin.PCG!(EST)
Elaslin.Write(EST)

print("Variables:   "); println(length(EST.xopt))
print("Iterations:  "); println(EST.iter)
print("Time:        "); println(EST.time)
