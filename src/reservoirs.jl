function pumpfun(t)
    if -20 < t <= 0
        f = (t+20) / 20
    elseif 0 < t < 10
        f = -(t-10) / 10
    else
        f = 0.0
    end
    return f
end

function model_A(du,u,p,t)
    function pumpfun(t)
        if -20 < t <= 0
            f = (t+20) / 20
        elseif 0 < t < 10
            f = -(t-10) / 10
        else
            f = 0.0
        end
        return f
    end

    k13,k23,k24,k32,k34,k42,k43,kh1,k2h,k3h,k4h = p
    k13 *= pumpfun(t)

    du[1] = kh1 * u[5] - k13 * (u[1] - u[3])
    du[2] = k42 * u[4] + k32 * u[3] - k23 * u[2] - k24 * u[2] - k2h * u[2]
    du[3] = k13 * (u[1] - u[3]) + k23 * u[2] + k43 * u[4] - k32 * u[3] - k34 * u[3] - k3h * u[3]
    du[4] = k24 * u[2] + k34 * u[3] - k42 * u[4] - k43 * u[4] - k4h * u[4]
    du[5] = k2h * u[2] + k3h * u[3] + k4h * u[4] - kh1 * u[5]
end
