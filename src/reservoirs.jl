function model(du,u,p,t)
    # indices of du, u, p:
    # 1:    heat bath
    # 2:    ground state
    # 3:    pumped mode
    # 4+:   additional modes

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

    # Number of Reservoirs
    N = length(du)

    # HEAT BATH
    # <-- pumped mode [3], additional modes[4+]
    # --> ground state [2]
    du[1] = -u[1] * p[1,2]
    for i = 3:N
        du[1] += u[i] * p[i,1]
    end

    # GROUND STATE
    # <-- heat bath [1]
    # --> pumped mode [3]
    du[2] = u[1] * p[1,2] - (u[2] - u[3]) * p[2,3] * pumpfun(t)

    # PUMPED MODE
    # <-- ground state [2], additional modes[4+]
    # --> heat bath [1], additional modes[4+] (in next block)
    # from ground state / to heat bath:
    du[3] = (u[2] - u[3]) * p[2,3] * pumpfun(t)

    # PUMPED MODE +
    # ADDITIONAL MODES
    for i = 4:N
        du[i] = 0.0
    end
    for n = 3:N
        du[n] -= u[n] * p[n,1] # to heat bath
        for i = 3:N
            n == i && continue
            du[n] += u[i] * p[i,n] # from other mode
            du[n] -= u[n] * p[n,i] # to other mode
        end
    end
end
