module Elaslin

using Mosek, SparseArrays, LinearAlgebra, Printf, FinancialToolbox

mutable struct Estructura
    Problema
    beta
    vtol
    alpha
    Vol
    TolGap
    PropMat
    dim
    Nnodo
    Nelem
    Ncond
    Mnodo
    Melem
    McndC
    McndV
    NGrLi
    GrLib
    Belem
    NumFu
    PCarg
    RCarg
    xopt
    copt
    iter
    time
end

function EmptyEstructura()
    Estructura("",0.0,0.0,0.0,0.0,0.0,[],0,0,0,0,[],[],[],[],0,[],[],0,[],[],[],0.0,0,0.0)
end

function Elaslin!(EST)
    dim = EST.dim
    Nnodo = EST.Nnodo
    Nelem = EST.Nelem
    Ncond = EST.Ncond
    PE = 1
    XX = 1
    PM = 2
    N1 = 4; N2 = 5
    NC = 2; ES = 3
    Belem = zeros(2*dim,Nelem)
    for ele = 1:Nelem
        Be = zeros(2*dim,1)
        Nod1 = EST.Melem[ele,N1]
        Nod2 = EST.Melem[ele,N2]
        V1 = EST.Mnodo[Nod1,XX:XX+dim-1]
        V2 = EST.Mnodo[Nod2,XX:XX+dim-1]
        dV = V2-V1
        L = norm(dV)
        dV = dV/L
        Ee = EST.PropMat[EST.Melem[ele,PM],PE]
        EAL = sqrt(Ee)/L
        Be[1:dim] = EAL*dV
        Be[dim+1:end] = -Be[1:dim]
        Belem[:, ele] = Be
    end
    GrLib = zeros(Int,dim*Nnodo)
    NumFu = maximum(EST.McndC[:,ES])
    PCarg = zeros(dim*Nnodo,NumFu)
    for cond = 1:Ncond
        Nod = EST.McndC[cond,NC]
        Est = EST.McndC[cond,ES]
        for GL = 1:dim
            tipo = EST.McndC[cond,3+GL]
            if (tipo==1)
                GrLib[dim*(Nod-1)+GL] = 1
            else
                PCarg[dim*(Nod-1)+GL,Est] += EST.McndV[cond,GL]
            end
        end
    end
    NGrLi = 0
    for gr = 1:dim*Nnodo
        if (GrLib[gr] == 1)
            GrLib[gr] = 0
        else
            NGrLi += 1
            GrLib[gr] = NGrLi
        end
    end
    PCarg = PCarg[GrLib.>0,:]
    EST.NGrLi = NGrLi
    EST.GrLib = GrLib
    EST.Belem = Belem
    EST.NumFu = NumFu
    EST.PCarg = PCarg
    return nothing
end

function PCG!(EST)
    EST.alpha = 2.0*log(normcdf(-EST.beta))
    dim = EST.dim
    Nelem = EST.Nelem
    Melem = EST.Melem
    Belem = EST.Belem
    NGrLi = EST.NGrLi
    GrLib = EST.GrLib
    NumFu = EST.NumFu
    PCarg = EST.PCarg
    RCarg = EST.RCarg
    RR = RCarg'*RCarg
    Vol = EST.Vol
    alpha = EST.alpha
    TolGap = EST.TolGap
    numvar = 3 + NumFu*NumFu + 3*NumFu
    numcon = 2*NGrLi*NumFu
    numcon += NumFu*(NumFu+1)÷2
    numcon += NumFu*NumFu
    numcon += NumFu*(NumFu+1)÷2
    numcon += NumFu*(NumFu+1)÷2
    numcon += NumFu*NumFu
    numcon += NumFu*(NumFu+1)÷2
    numcon += NumFu*(NumFu-1)÷2
    numcon += NumFu
    numcon += NumFu
    numcon += 2
    bkc = [MSK_BK_FX for i = 1:numcon]
    blc = vcat(PCarg[:], zeros(NGrLi*NumFu), zeros(4*(NumFu*(NumFu+1)÷2)), zeros((NumFu*(NumFu-1)÷2)), zeros(2*NumFu*NumFu), zeros(2*NumFu), 0.0, Vol)
    buc = vcat(PCarg[:], zeros(NGrLi*NumFu), zeros(4*(NumFu*(NumFu+1)÷2)), zeros((NumFu*(NumFu-1)÷2)), zeros(2*NumFu*NumFu), zeros(2*NumFu), 0.0, Vol)
    Ai = zeros(Int,9*NumFu + 3*NumFu*NumFu + 3)
    Aj = zeros(Int,9*NumFu + 3*NumFu*NumFu + 3)
    Av = zeros(9*NumFu + 3*NumFu*NumFu + 3)
    Bi = zeros(Int,2*NumFu + 4*NumFu*NumFu + Nelem*(1 + 4*dim*NumFu + 2*NumFu + 4*NumFu*NumFu))
    Bj = zeros(Int,2*NumFu + 4*NumFu*NumFu + Nelem*(1 + 4*dim*NumFu + 2*NumFu + 4*NumFu*NumFu))
    Bk = zeros(Int,2*NumFu + 4*NumFu*NumFu + Nelem*(1 + 4*dim*NumFu + 2*NumFu + 4*NumFu*NumFu))
    Bl = zeros(Int,2*NumFu + 4*NumFu*NumFu + Nelem*(1 + 4*dim*NumFu + 2*NumFu + 4*NumFu*NumFu))
    Bv = zeros(2*NumFu + 4*NumFu*NumFu + Nelem*(1 + 4*dim*NumFu + 2*NumFu + 4*NumFu*NumFu))
    PosA, PosB = 0, 0
    N1 = 4; N2 = 5
    for k = 1:Nelem
        nl = 0
        for NN = N1:N2
            nl += 1
            nod = Melem[k,NN]
            for gr = 1:dim
                GrG = GrLib[dim*(nod-1)+gr]
                if (GrG > 0)
                    for i = 1:NumFu
                        PosB += 1
                        Bi[PosB] = (i-1)*NGrLi + GrG
                        Bj[PosB] = k
                        Bk[PosB] = 1 + i
                        Bl[PosB] = 1
                        Bv[PosB] = 0.5*Belem[dim*(nl-1)+gr,k]
                        PosB += 1
                        Bi[PosB] = NGrLi*NumFu + (i-1)*NGrLi + GrG
                        Bj[PosB] = k
                        Bk[PosB] = 1 + NumFu + i
                        Bl[PosB] = 1
                        Bv[PosB] = 0.5*Belem[dim*(nl-1)+gr,k]
                    end
                end
            end
        end
    end
    Eqi1 = 2*NGrLi*NumFu
    Eqi2 = Eqi1 + NumFu*(NumFu+1)÷2
    Eqi3 = Eqi2 + NumFu*NumFu
    Eqi4 = Eqi3 + NumFu*(NumFu+1)÷2
    Eqi5 = Eqi4 + NumFu*(NumFu+1)÷2
    Eqi6 = Eqi5 + NumFu*NumFu
    Eqi7 = Eqi6 + NumFu*(NumFu+1)÷2
    Eqi8 = Eqi7 + NumFu*(NumFu-1)÷2
    Eqi9 = Eqi8 + NumFu
    Eqi0 = Eqi9 + NumFu
    for i = 1:NumFu
        j = i
        PosA += 1
        Ai[PosA] = Eqi1 + (i-1)*i÷2 + j
        Aj[PosA] = 2
        Av[PosA] = -1.0
        PosB += 1
        Bi[PosB] = Eqi1 + (i-1)*i÷2 + j
        Bj[PosB] = Nelem + 1
        Bk[PosB] = i
        Bl[PosB] = j
        Bv[PosB] = 1.0
        for j = 1:i-1
            PosB += 1
            Bi[PosB] = Eqi1 + (i-1)*i÷2 + j
            Bj[PosB] = Nelem + 1
            Bk[PosB] = i
            Bl[PosB] = j
            Bv[PosB] = 0.5
        end
        for j = 1:NumFu
            PosA += 1
            Ai[PosA] = Eqi2 + (i-1)*NumFu + j
            Aj[PosA] = 2
            Av[PosA] = -RCarg[j,i]
            PosB += 1
            Bi[PosB] = Eqi2 + (i-1)*NumFu + j
            Bj[PosB] = Nelem + 1
            Bk[PosB] = NumFu + i
            Bl[PosB] = j
            Bv[PosB] = 0.5
        end
        j = i
        PosA += 1
        Ai[PosA] = Eqi3 + (i-1)*i÷2 + j
        Aj[PosA] = 3
        Av[PosA] = -1.0
        PosA += 1
        Ai[PosA] = Eqi3 + (i-1)*i÷2 + j
        Aj[PosA] = 2
        Av[PosA] = -RR[i,j]
        PosB += 1
        Bi[PosB] = Eqi3 + (i-1)*i÷2 + j
        Bj[PosB] = Nelem + 1
        Bk[PosB] = NumFu + i
        Bl[PosB] = NumFu + j
        Bv[PosB] = 1.0
        for j = 1:i-1
            PosA += 1
            Ai[PosA] = Eqi3 + (i-1)*i÷2 + j
            Aj[PosA] = 2
            Av[PosA] = -RR[i,j]
            PosB += 1
            Bi[PosB] = Eqi3 + (i-1)*i÷2 + j
            Bj[PosB] = Nelem + 1
            Bk[PosB] = NumFu + i
            Bl[PosB] = NumFu + j
            Bv[PosB] = 0.5
        end
        j = i
        PosA += 1
        Ai[PosA] = Eqi4 + (i-1)*i÷2 + j
        Aj[PosA] = 2
        Av[PosA] = -1.0
        PosB += 1
        Bi[PosB] = Eqi4 + (i-1)*i÷2 + j
        Bj[PosB] = Nelem + 2
        Bk[PosB] = i
        Bl[PosB] = j
        Bv[PosB] = 1.0
        for j = 1:i-1
            PosB += 1
            Bi[PosB] = Eqi4 + (i-1)*i÷2 + j
            Bj[PosB] = Nelem + 2
            Bk[PosB] = i
            Bl[PosB] = j
            Bv[PosB] = 0.5
        end
        for j = 1:NumFu
            PosA += 1
            Ai[PosA] = Eqi5 + (i-1)*NumFu + j
            Aj[PosA] = 3 + (j-1)*NumFu + i
            Av[PosA] = -1.0
            PosB += 1
            Bi[PosB] = Eqi5 + (i-1)*NumFu + j
            Bj[PosB] = Nelem + 2
            Bk[PosB] = NumFu + i
            Bl[PosB] = j
            Bv[PosB] = 0.5
        end
        j = i
        PosA += 1
        Ai[PosA] = Eqi6 + (i-1)*i÷2 + j
        Aj[PosA] = 3 + (i-1)*NumFu + j
        Av[PosA] = -1.0
        PosB += 1
        Bi[PosB] = Eqi6 + (i-1)*i÷2 + j
        Bj[PosB] = Nelem + 2
        Bk[PosB] = NumFu + i
        Bl[PosB] = NumFu + j
        Bv[PosB] = 1.0
        for j = 1:i-1
            PosB += 1
            Bi[PosB] = Eqi6 + (i-1)*i÷2 + j
            Bj[PosB] = Nelem + 2
            Bk[PosB] = NumFu + i
            Bl[PosB] = NumFu + j
            Bv[PosB] = 0.5
        end
    end
    for k = 1:Nelem
        for i = 1:NumFu
            j = i
            PosB += 1
            Bi[PosB] = Eqi1 + (i-1)*i÷2 + j
            Bj[PosB] = k
            Bk[PosB] = 1 + i
            Bl[PosB] = 1 + j
            Bv[PosB] = 1.0
            PosB += 1
            Bi[PosB] = Eqi2 + (i-1)*NumFu + j
            Bj[PosB] = k
            Bk[PosB] = 1 + NumFu + i
            Bl[PosB] = 1 + j
            Bv[PosB] = 0.5
            PosB += 1
            Bi[PosB] = Eqi3 + (i-1)*i÷2 + j
            Bj[PosB] = k
            Bk[PosB] = 1 + NumFu + i
            Bl[PosB] = 1 + NumFu + j
            Bv[PosB] = 1.0
            PosB += 1
            Bi[PosB] = Eqi4 + (i-1)*i÷2 + j
            Bj[PosB] = k
            Bk[PosB] = 1 + i
            Bl[PosB] = 1 + j
            Bv[PosB] = 1.0
            PosB += 1
            Bi[PosB] = Eqi5 + (i-1)*NumFu + j
            Bj[PosB] = k
            Bk[PosB] = 1 + NumFu + i
            Bl[PosB] = 1 + j
            Bv[PosB] = 0.5
            PosB += 1
            Bi[PosB] = Eqi6 + (i-1)*i÷2 + j
            Bj[PosB] = k
            Bk[PosB] = 1 + NumFu + i
            Bl[PosB] = 1 + NumFu + j
            Bv[PosB] = 1.0
            for j = 1:i-1
                PosB += 1
                Bi[PosB] = Eqi1 + (i-1)*i÷2 + j
                Bj[PosB] = k
                Bk[PosB] = 1 + i
                Bl[PosB] = 1 + j
                Bv[PosB] = 0.5
                PosB += 1
                Bi[PosB] = Eqi2 + (i-1)*NumFu + j
                Bj[PosB] = k
                Bk[PosB] = 1 + NumFu + i
                Bl[PosB] = 1 + j
                Bv[PosB] = 0.5
                PosB += 1
                Bi[PosB] = Eqi3 + (i-1)*i÷2 + j
                Bj[PosB] = k
                Bk[PosB] = 1 + NumFu + i
                Bl[PosB] = 1 + NumFu + j
                Bv[PosB] = 0.5
                PosB += 1
                Bi[PosB] = Eqi4 + (i-1)*i÷2 + j
                Bj[PosB] = k
                Bk[PosB] = 1 + i
                Bl[PosB] = 1 + j
                Bv[PosB] = 0.5
                PosB += 1
                Bi[PosB] = Eqi5 + (i-1)*NumFu + j
                Bj[PosB] = k
                Bk[PosB] = 1 + NumFu + i
                Bl[PosB] = 1 + j
                Bv[PosB] = 0.5
                PosB += 1
                Bi[PosB] = Eqi6 + (i-1)*i÷2 + j
                Bj[PosB] = k
                Bk[PosB] = 1 + NumFu + i
                Bl[PosB] = 1 + NumFu + j
                Bv[PosB] = 0.5
            end
            for j = i+1:NumFu
                PosB += 1
                Bi[PosB] = Eqi2 + (i-1)*NumFu + j
                Bj[PosB] = k
                Bk[PosB] = 1 + NumFu + i
                Bl[PosB] = 1 + j
                Bv[PosB] = 0.5
                PosB += 1
                Bi[PosB] = Eqi5 + (i-1)*NumFu + j
                Bj[PosB] = k
                Bk[PosB] = 1 + NumFu + i
                Bl[PosB] = 1 + j
                Bv[PosB] = 0.5
            end
        end
    end
    for i = 1:NumFu
        j = i
        PosA += 1
        Ai[PosA] = Eqi8 + i
        Aj[PosA] = 3 + NumFu*NumFu + 3*(i-1) + 1
        Av[PosA] = -1.0
        PosA += 1
        Ai[PosA] = Eqi8 + i
        Aj[PosA] = 3 + (i-1)*NumFu + j
        Av[PosA] = 1.0
        for j = 1:i-1
            PosA += 1
            Ai[PosA] = Eqi7 + (i-2)*(i-1)÷2 + j
            Aj[PosA] = 3 + (j-1)*NumFu + i
            Av[PosA] = 1.0
        end
    end
    for i = 1:NumFu
        PosA += 1
        Ai[PosA] = Eqi9 + i
        Aj[PosA] = 2
        Av[PosA] = -1.0
        PosA += 1
        Ai[PosA] = Eqi9 + i
        Aj[PosA] = 3 + NumFu*NumFu + 3*(i-1) + 2
        Av[PosA] = 1.0
    end
    for i = 1:NumFu
        PosA += 1
        Ai[PosA] = Eqi0 + 1
        Aj[PosA] = 3 + NumFu*NumFu + 3*(i-1) + 3
        Av[PosA] = 1.0
    end
    PosA += 1
    Ai[PosA] = Eqi0 + 1
    Aj[PosA] = 1
    Av[PosA] = 1.0
    PosA += 1
    Ai[PosA] = Eqi0 + 1
    Aj[PosA] = 2
    Av[PosA] = alpha
    PosA += 1
    Ai[PosA] = Eqi0 + 1
    Aj[PosA] = 3
    Av[PosA] = -1.0
    for k = 1:Nelem
        PosB += 1
        Bi[PosB] = Eqi0 + 2
        Bj[PosB] = k
        Bk[PosB] = 1
        Bl[PosB] = 1
        Bv[PosB] = 1.0
    end
    A = sparse(Ai, Aj, Av, numcon, numvar)
    Bi = Bi[1:PosB]
    Bj = Bj[1:PosB]
    Bk = Bk[1:PosB]
    Bl = Bl[1:PosB]
    Bv = Bv[1:PosB]
    maketask() do task
        appendvars(task, numvar)
        putvarboundsliceconst(task, 1, numvar+1, MSK_BK_FR, -Inf, Inf)
        appendcons(task, numcon)
        putconboundslice(task, 1, numcon+1, bkc, blc, buc)
        putacolslice(task, 1, numvar+1, A)
        afei = getnumafe(task) + 1
        appendafes(task, 3*NumFu)
        putafefentrylist(task,
            [i               for i = 1:3*NumFu],
            [i+3+NumFu*NumFu for i = 1:3*NumFu],
            ones(3*NumFu))
        for i = 1:NumFu
            expdomain = appendprimalexpconedomain(task)
            appendaccseq(task,
                expdomain,
                afei,
                nothing)
            afei += 3
        end
        dimbarvar = zeros(Int, Nelem + 2) .+ (1 + 2*NumFu)
        dimbarvar[end-1:end] = [2*NumFu, 2*NumFu] 
        appendbarvars(task, dimbarvar)
        putbarablocktriplet(task, Bi, Bj, Bk, Bl, Bv)
        putcj(task, 1, 1.0)
        putobjsense(task, MSK_OBJECTIVE_SENSE_MINIMIZE)
        if (TolGap > 0.0)
            putdouparam(task, MSK_DPAR_INTPNT_CO_TOL_REL_GAP, TolGap)
        end
        optimize(task)
        solutionsummary(task, MSK_STREAM_MSG)
        prosta = getprosta(task, MSK_SOL_ITR)
        solsta = getsolsta(task, MSK_SOL_ITR)
        ct = getxx(task, MSK_SOL_ITR)
        c = ct[1]
        xx = zeros(Nelem)
        for j = 1:Nelem
            barx = getbarxj(task, MSK_SOL_ITR, j)
            xx[j] = barx[1]
        end
        EST.xopt = xx
        EST.copt = c
        EST.time = getdouinf(task, MSK_DINF_INTPNT_TIME)
        EST.iter = getintinf(task, MSK_IINF_INTPNT_ITER)
        if solsta == MSK_SOL_STA_OPTIMAL
        elseif solsta == MSK_SOL_STA_DUAL_INFEAS_CER
            println("Primal or dual infeasibility.")
        elseif solsta == MSK_SOL_STA_PRIM_INFEAS_CER
            println("Primal or dual infeasibility.")
        elseif solsta == MSK_SOL_STA_UNKNOWN
            println("Unknown solution status")
        else
            println("Other solution status")
        end
    end
    return nothing
end

function Malla2D!(EST, xi, xf, yi, yf, nx, ny, lv)
    dx = (xf-xi)/nx
    dy = (yf-yi)/ny
    XX = 1; YY = 2
    Tol = 1e-5
    Mnodo = zeros((nx+1)*(ny+1),2)
    Nnodo = 0
    for x = xi:dx:xf
        for y = yi:dy:yf
            Nnodo += 1
            Mnodo[Nnodo,:] = [x y]
        end
    end
    Melem = zeros(Int,0,5)
    Nelem = 0
    for nod1 = 1:Nnodo-1
        x1 = Mnodo[nod1,XX]
        y1 = Mnodo[nod1,YY]
        for nod2 = nod1+1:Nnodo
            x2 = Mnodo[nod2,XX]
            y2 = Mnodo[nod2,YY]
            ix = round(Int,(x2-x1)/dx)
            iy = round(Int,(y2-y1)/dy)
            if (abs(ix) <= lv) && (abs(iy) <= lv)
                v1 = [x2-x1, y2-y1]
                L1 = norm(v1)
                v1 = v1/L1
                solapa = false
                for nod3 = 1:Nnodo
                    x3 = Mnodo[nod3,XX]
                    y3 = Mnodo[nod3,YY]
                    v2 = [x3-x1, y3-y1]
                    L2 = norm(v2)
                    if (L2 > Tol)
                        v2 = v2/L2
                    end
                    if (norm(v2-v1) < Tol)
                        if (L2 > Tol) && (L2 < L1 - Tol)
                            solapa = true
                        end
                    end
                end
                if ~solapa
                    Nelem += 1
                    Melem = vcat(Melem, [Nelem 1 Nelem nod1 nod2])
                end
            end
        end
    end
    EST.Nnodo = Nnodo
    EST.Mnodo = Mnodo
    EST.Nelem = Nelem
    EST.Melem = Melem
    return nothing
end

function Malla3D!(EST, xi, xf, yi, yf, zi, zf, nx, ny, nz, lv)
    dx = (xf-xi)/nx
    dy = (yf-yi)/ny
    dz = (zf-zi)/nz
    XX = 1; YY = 2; ZZ = 3
    Tol = 1e-5
    Mnodo = zeros((nx+1)*(ny+1)*(nz+1),3)
    Nnodo = 0
    for x = xi:dx:xf
        for y = yi:dy:yf
            for z = zi:dz:zf
                Nnodo += 1
                Mnodo[Nnodo,:] = [x y z]
            end
        end
    end
    Melem = zeros(Int,0,5)
    Nelem = 0
    for nod1 = 1:Nnodo-1
        x1 = Mnodo[nod1,XX]
        y1 = Mnodo[nod1,YY]
        z1 = Mnodo[nod1,ZZ]
        for nod2 = nod1+1:Nnodo
            x2 = Mnodo[nod2,XX]
            y2 = Mnodo[nod2,YY]
            z2 = Mnodo[nod2,ZZ]
            ix = round(Int,(x2-x1)/dx)
            iy = round(Int,(y2-y1)/dy)
            iz = round(Int,(z2-z1)/dz)
            if (abs(ix) <= lv) && (abs(iy) <= lv) && (abs(iz) <= lv)
                v1 = [x2-x1, y2-y1, z2-z1]
                L1 = norm(v1)
                v1 = v1/L1
                solapa = false
                for nod3 = 1:Nnodo
                    x3 = Mnodo[nod3,XX]
                    y3 = Mnodo[nod3,YY]
                    z3 = Mnodo[nod3,ZZ]
                    v2 = [x3-x1, y3-y1, z3-z1]
                    L2 = norm(v2)
                    if (L2 > Tol)
                        v2 = v2/L2
                    end
                    if (norm(v2-v1) < Tol)
                        if (L2 > Tol) && (L2 < L1 - Tol)
                            solapa = true
                        end
                    end
                end
                if ~solapa
                    Nelem += 1
                    Melem = vcat(Melem, [Nelem 1 Nelem nod1 nod2])
                end
            end
        end
    end
    EST.Nnodo = Nnodo
    EST.Mnodo = Mnodo
    EST.Nelem = Nelem
    EST.Melem = Melem
    return nothing
end

function Write(EST)
    dim = EST.dim
    archivo = open(EST.Problema*".m", "w")
    @printf(archivo, "vtol = %14.6e;\n\n", EST.vtol)
    @printf(archivo, "PropGeo = [\n")
    for i = 1:EST.Nelem
        @printf(archivo, "%14.6e\n", EST.xopt[i])
    end
    @printf(archivo, "];\n\n")
    @printf(archivo, "Mnodo = [\n")
    for i = 1:EST.Nnodo
        for j = 1:dim
            @printf(archivo, "%14.6e", EST.Mnodo[i,j])
        end
        @printf(archivo, "\n")
    end
    @printf(archivo, "];\n\n")
    @printf(archivo, "Melem = [\n")
    for i = 1:EST.Nelem
        for j = 1:5
            @printf(archivo, "%8i", EST.Melem[i,j])
        end
        @printf(archivo, "\n")
    end
    @printf(archivo, "];\n")
    close(archivo)
    return nothing
end

end
