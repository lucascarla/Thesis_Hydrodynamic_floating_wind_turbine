
%% EXPANSION RATIO
x=[0.1 0.2 0.3 0.4 0.5];
y=[0.9962 0.9945 0.99 0.9846 0.9814];
h=plot(x,y,'Linewidth',1)
h.LineStyle='-';
h.Marker='o';
h.Color='black';
xlabel('Expansion Ratio')
ylabel('Computed A/Target A')

figure

y=[536.76 441.56 398.02 376.63 352.26];

h=plot(x,y,'Linewidth',1)
h.LineStyle='-';
h.Marker='o';
h.Color='red';
xlabel('Expansion Ratio')
ylabel('Execution Time [s]')

%% RELAXATION LENGTH

x=[1 2 3 4];
y=[0.9434 0.9924 0.9469 0.995];

h=plot(x,y,'Linewidth',1)
set(gca, 'XTick',1:4, 'XTickLabel',{'1 & 1' '1 & 2' '2 & 1' '2 & 2'})
h.LineStyle='-';
h.Marker='o';
h.Color='black';
xlabel('Generation (first digit) and Absorption (second digit) zone length compared to wavelength')
ylabel('Computed A/Target A')

figure

y=[316.68 351.65 356.94 379.46];

h=plot(x,y,'Linewidth',1)
set(gca, 'XTick',1:4, 'XTickLabel',{'1 & 1' '1 & 2' '2 & 1' '2 & 2'})
h.LineStyle='-';
h.Marker='o';
h.Color='red';
xlabel('Generation (first digit) and Absorption (second digit) zone length compared to wavelength')
ylabel('Execution Time [s]')

%% FREE SURFACE LEVEL REFINEMENT
x=[0 1 2 3 4];
y=[0.7618 0.9064 0.9526 0.9653 0.9763];
h=plot(x,y,'Linewidth',1)
set(gca, 'XTick',0:4)
h.LineStyle='-';
h.Marker='o';
h.Color='black';
xlabel('Free Surface Level Refinement')
ylabel('Computed A/Target A')

figure

y=[67.05 73.65 97.32 143.81 319.45];
h=plot(x,y,'Linewidth',1)
set(gca, 'XTick',0:4)
h.LineStyle='-';
h.Marker='o';
h.Color='red';
xlabel('free surface level refinement')
ylabel('Execution Time [s]')


%% Maximum residual in the Poisson iterative solution for pressure

x=[1 2 3 4 5];
y=[0.9922 0.9925 0.9927 0.9927 0.9929];

h=plot(x,y,'Linewidth',1)
set(gca, 'XTick',1:5, 'XTickLabel',{'1e-4 & 1e-5' '1e-5 & 1e-6' '1e-6 & 1e-7' '1e-7 & 1e-8' '1e-5 & 1e-7'})
h.LineStyle='-';
h.Marker='o';
h.Color='black';
xlabel('p\_rgh initial & p\_rgh final max residuals')
ylabel('Computed A/Target A')

figure

y=[230 254 260 278.5 245];
h=plot(x,y,'Linewidth',1)
set(gca, 'XTick',1:5, 'XTickLabel',{'1e-4 & 1e-5' '1e-5 & 1e-6' '1e-6 & 1e-7' '1e-7 & 1e-8' '1e-5 & 1e-7'})
h.LineStyle='-';
h.Marker='o';
h.Color='red';
xlabel('p\_rgh initial & p\_rgh final max residuals')
ylabel('Execution Time [s]')


