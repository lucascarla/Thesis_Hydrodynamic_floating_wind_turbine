%% Moody model file for 4 catenary chains

% script for case setup %
dimensionNumber = 3;
waterLevel = 0;         % [m]       z-coordinate of mean water level
waterDensity = 1025.0;     % [kg/m??]   Density of water
airDensity = 1.225;         % [kg/m??]   Density of air
gravity = 1;

time.start  = 0;
time.end = 600;
time.cfl = 0.9;
time.scheme = 'RK3';

print.dt = 0.02; % output interval.

%% extra quadpoints used for increased ground contact performance.
% it makes the difference between stable and non-stable results. 
numLib.qPointsAdded = 10;

%-------------------------------------------------------------------------%
%--------------------------- Ground model input  -------------------------%
%-------------------------------------------------------------------------%
ground.type = 'springDampGround';
ground.level = -200;
ground.dampingCoeff = 1.0;
ground.frictionCoeff = 0.1;
ground.vc = 0.01;
ground.stiffness = 300.0e6;

%-------------------------------------------------------------------------%
%---------------------------- Type definition ----------------------------%
%-------------------------------------------------------------------------%
cableType1.diameter = 0.0766;
cableType1.gamma0 = 113.35;
cableType1.CDn = 1.6;
cableType1.CDt = 0.05;
cableType1.CMn = 0.8;
cableType1.CMt = 0.25;
cableType1.materialModel.type = 'biLinear';
cableType1.materialModel.EA =  753600000;
        
%-------------------------------------------------------------------------%
%------------------------------- Geometry --------------------------------%

% note that unconnected  vertices are ignored. 
% vertice 5-8 are externally controlled. No need to specify exact value.

vertexLocations = {
                       1    [418.800   725.383    -200.000   ];
                       2    [-837.60   0.0        -200.000];
                       3    [418.800  -725.383    -200.000   ];
                       4    [42.434    35.393     -14.000];
                       5    [-18.87    0.000      -14.000   ];
                       6    [42.434   -35.393     -14.000];                       
                  };
                   
                         
%----- Object definitions -----%
cable1.typeNumber = 1;
cable1.startVertex = 1; %
cable1.endVertex = 4; %
cable1.length = 835.35; % 
cable1.IC.type = 'CatenaryStatic';
cable1.N = 20; % 
% 
% Copy remaining info from cable1. short hand
cable2=cable1;
cable2.startVertex = 2; 
cable2.endVertex = 5; 
%
cable3=cable1;
cable3.startVertex = 3; 
cable3.endVertex = 6; 
%-------------------------------------------------------------------------%
%                           Boundary conditions                           %
%-------------------------------------------------------------------------%

% Four anchors defined by vertexLocations
bc1.vertexNumber = 1;
bc1.type = 'dirichlet';
bc1.mode = 'fixed';
%
bc2=bc1;
bc2.vertexNumber = 2;
% 
bc3=bc1;
bc3.vertexNumber = 3;

bc5.vertexNumber = 4;
bc5.type = 'dirichlet';
bc5.mode = 'externalPoint';

bc6=bc5;
bc6.vertexNumber = 5;

bc7=bc5;
bc7.vertexNumber = 6;

% --- API conectivitity --- %
API.bcNames = {'bc5','bc6','bc7'};
API.reboot= 'no';
API.syncOutput = 1; 
API.staggerTimeFraction= 0.5;
API.output = 'Mooring/results';


% To remove initial fluctuations due to very stiff ground
% staticSolver.relax = 1;
% staticSolver.relaxFactor = 0.9;

% ===== END OF FILE ===== %
