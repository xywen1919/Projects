%%%% -*- Mode: Matlab -*-

%%%% reprsim.m --
%%%% Code to run the simple repressilator system.
%%%% The actual ODE system in packaged as a function in the
%%%%'repressilator' M-file.

%%%% Marco Antoniotti (c) 2010-2019
%%%% Computational Biology Course
%%%% DISCo Universita` degli Studi di Milano Bicocca

a    = 220;
a0   = 1;
beta = 2;
n    = 1.8;

repr = @(t, y) repressilator(t, y, a, a0, beta, n); % A 'closure'.

y0 = [10; 1; 1; 1; 1; 1];

[ts, ys] = ode45(repr, [0 100], y0);

plot(ts, ys(:, 2), '-', ts, ys(:, 4), '-', ts, ys(:, 6));
xlabel('Time');
ylabel('Protein concentration');
title('Protein concentrations of repressilator system');

figure
plot(ts, ys(:, 1), '-', ts, ys(:, 3), '-', ts, ys(:, 5));
xlabel('Time');
ylabel('mRNA concentration');
title('mRNA concentrations of repressilator system');

figure
plotyy(ts, ys(:, 1), ts, ys(:, 2));
xlabel('Time');
ylabel('lacl and Lacl measurements');
title('lacl and Lacl delays in repressilator system');



%%%% end of file -- reprsim.m --
