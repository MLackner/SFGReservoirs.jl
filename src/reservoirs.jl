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

    du[1] = p[5,1] * u[5] - p[1,3] * pumpfun(t) * (u[1] - u[3])
    du[2] = p[4,2] * u[4] + p[3,2] * u[3] - p[2,3] * u[2] - p[2,4] * u[2] - p[2,5] * u[2]
    du[3] = p[1,3] * pumpfun(t) * (u[1] - u[3]) + p[2,3] * u[2] + p[4,3] * u[4] - p[3,2] * u[3] - p[3,4] * u[3] - p[3,5] * u[3]
    du[4] = p[2,4] * u[2] + p[3,4] * u[3] - p[4,2] * u[4] - p[4,3] * u[4] - p[4,5] * u[4]
    du[5] = p[2,5] * u[2] + p[3,5] * u[3] + p[4,5] * u[4] - p[5,1] * u[5]
end
