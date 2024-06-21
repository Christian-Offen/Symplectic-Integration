# fixpoint iterations
function FixIter(fun,z,minStepSize)
    while true
        zold = z
        z = fun(z)
        norm(zold - z) > minStepSize || break
    end
    return z
end


# definition of several numerical integrators

function Euler(h,z)   # not symplectic
    return z+h*XH(z)
end

function SympEuler(h,z)  # symplectic
    R(deltaQ) = h*Hp(z[1]+deltaQ,z[2])
    deltaQ=FixIter(R,0.,1e-12)
    q = z[1]+deltaQ
    p = z[2]-h*Hq(q,z[2])
    return [q;p]
end


function MidPoint(h,z)   # symplectic
    R(delta) = h*XH(z + 1/2*delta)
    delta=FixIter(R,zeros(size(z)),1e-12)
    return z + delta
end

# repreated application of integrator
function Integrator(numFlow,z,steps)
    Z = zeros(2,steps+1)
    Z[:,1] = z
    for k = 1:size(Z,2)-1
        Z[:,k + 1] = numFlow(Z[:,k])
    end
    return Z
end