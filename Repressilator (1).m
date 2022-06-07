%%%% -*- Mode: Matlab -*-

%%%% repressilator.m --
%%%% The actual repressilator ODE system packaged as a function
%%%% (to be passed, properly parameterized) to the ODE integrator.

%%%% Marco Antoniotti (c) 2010-2019
%%%% Computational Biology Course
%%%% DISCo Universita` degli Studi di Milano Bicocca

function dy = repressilator(t, y, a, a0, beta, n)
    dy = zeros(6, 1);
    lacl = y(1);
    Lacl = y(2);
    tetR = y(3);
    TetR = y(4);
    cl   = y(5);
    Cl   = y(6);
    
    dy(1) = -lacl + (a / (1 + TetR^n)) + a0;
    dy(2) = -beta * (Lacl - lacl);
    
    dy(3) = -tetR + (a / (1 + Cl^n)) + a0;
    dy(4) = -beta * (TetR - tetR);
    
    dy(5) = -cl   + (a / (1 + Lacl^n)) + a0;
    dy(6) = -beta * (Cl - cl);
end

%%%% end of file -- reprissilator.m --
