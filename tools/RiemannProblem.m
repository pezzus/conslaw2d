function [Q] = RiemannProblem(x,t,rl,ul,pl,rr,ur,pr)
    gamma = 1.4;
    cl = sqrt(gamma*pl/rl);
    cr = sqrt(gamma*pr/rr);
    [um,pm,cm1,cm2] = RI_solve(cl,pl,ul,cr,pr,ur,gamma);
    rm1 = gamma*pm/cm1^2;
    rm2 = gamma*pm/cm2^2;
    x1 = (ul-cl)*t;
    x2 = (um-cm1)*t;
    x3 = um*t;
    x4 = (rm2*um-rr*ur)/(rm2-rr)*t;
    if ( x < x1 )
        Q = [rl, ul, pl];
    else
        if ( x < x2 )
            Q = [ -(rl-rm1)/(x2-x1)*(x-x1)+rl, -(ul-um)/(x2-x1)*(x-x1)+ul, -(pl-pm)/(x2-x1)*(x-x1)+pl ];
        else
            if ( x < x3 )
                Q = [rm1, um, pm];
            else
                if ( x < x4 )
                    Q = [rm2, um, pm];
                else
                    Q = [rr, ur, pr];
                end
            end
        end
    end
end
