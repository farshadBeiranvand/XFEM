
clc
clear all
format shortG
% Define plate properties
h=1; % Thickness of plate
w=120/2; % Half-width of plate in the x-direction
L=120/2; % Half-length of plate 

%%% Define crack geometry
crackedplate=1; % Indicates if the plate is cracked (1) or not (0)
lc=60;  % Length of crack
tc=0; % Crack angle to X-axis (degrees)   -90<tc<=90
x1c=30.001; % x-coordinate of the 1st point of crack (start from left to right)
y1c=60.001; % y-coordinate of the 1st point of crack
neltip=2 ; % Number of crack tips: 1 for edge crack, 2 for center crack

%%%% Define material properties
isort=1; % 1 for isotropic material, 2 for orthotropic material
E=1e7; % Young's modulus
v=.3; % Poisson's ratio

E1=73e9; % Young's modulus in direction 1
E2=0.1*E1; % Young's modulus in direction 2
v12=.3; % Poisson's ratio
v21=.1*v12; 
G12=0.6*E2; % Shear modulus
G13=0.5*E2;
G23=G13;

%% Define boundary and force conditions
BCtype=['s' 's' 's' 's'];    % Boundary condition type of edges [Right Down Left Up] % 'c'=clamped 's'=simply supported 'f'=free

% Foundation properties
KWbar=0;   % 0 for no foundation
KWfound=0;
GPfound=0; % 0 for no foundation

nLambda=1;  % Number of eigenvalues

% Forces (positive magnitude is compressive)
Nxxlo=0; Lfx=0;    % Load factor in the x-direction [2,1.5,1,0.5,0=uniform]
Nyylo=-1; Lfy=0;    % Load factor in the y-direction
Nxylo=0;
% Non-local forces
[Nxxnl, Nyynl, Nxynl] = nonlocal_forces(Lfx, Lfy, Nxxlo, Nyylo);

% Kappa: shear correction factor, Miu: small dimension parameter, used in stiffness matrices
kapa=5/6; miu=0;
[As, Ds, Ainp, C, NG2, NG3] = NGA(kapa, h, isort, E, v, E1, E2, v12, v21, G12, G13, G23, Nxxnl, Nxynl, Nyynl);

%% Mesh generation
Q=8; % 4 for Q4 mesh, 8 for Q8
NYE=15; % Number of elements in y (MUST BE ODD)
NXE=15; % Number of elements in x (MUST BE ODD)
dhx = 2*w/NXE; % Element size in the x-direction
dhy = 2*L/NYE; % Element size in the y-direction
[nnd, connec, nel, geom] = getMesh(NXE, NYE, dhx, dhy); % nel: number of elements, nnd: number of nodes

[x1c, y1c, x2c, y2c, tipnodes, heavynodes, tipEl, heavyEl, heavytipEl, neighbortipEl, normEl, nodvalue] = typeELnod...
    (geom, nnd, nel, x1c, y1c, lc, tc, neltip, connec, Q, crackedplate);

checkcracktips(x1c, y1c, x2c, y2c, geom);
% plot2D(NXE, NYE, heavyEl, neighbortipEl, heavytipEl, tipEl, x1c, x2c, y1c, y2c, connec, geom);

% In-plane analysis
[gdofIN] = GlobalDof(nnd, nodvalue, 1);
KIN = sparse(gdofIN, gdofIN); % Initialize global stiffness matrix for in-plane analysis
Bforgauss = {};
globxyGp = [];
elementDofIN = {};
for p = 1:nel
    [xgp, ygp, wgp, Kinp, B, Bforgauss, elementDof, nDof, globxyGp] = StiffnessIN(Q, p, tipnodes, heavynodes, normEl, heavyEl,...
        heavytipEl, neighbortipEl, tipEl, geom, connec, nnd, nodvalue, x1c, y1c, x2c, y2c, tc, dhx, dhy, h, Ainp, kapa, Bforgauss, globxyGp);
    elementDofIN{p} = elementDof;
    
    KIN(elementDof, elementDof) = KIN(elementDof, elementDof) + Kinp;
end
KK = KIN;

% Solve KIN*F=d
force = loads2(Q, gdofIN, geom, w, L, NYE, NXE, Nxxlo, Nxylo, Nyylo, Lfx, Lfy);
displacements = solveEq(Q, NXE, NYE, KIN, gdofIN, force);
newgeom = newgeom(displacements, nodvalue, heavynodes, tipnodes, geom, nnd);
plot2D(NXE, NYE, heavyEl, neighbortipEl, heavytipEl, tipEl, x1c, x2c, y1c, y2c, connec, newgeom);

% Buckling analysis
[gdof] = GlobalDof(nnd, nodvalue, 2);
K = sparse(gdof, gdof); % Initialize global stiffness matrix for buckling analysis
KG = sparse(gdof, gdof); % Initialize geometric stiffness matrix
sigma = {};
sigmanGP = [];
for p = 1:nel
    elementDof_inplane = elementDofIN{p};
    [Kbend, Kgeo, Kshear, elementDof, sigma, sigmanGP] = StiffnessOUT(Q, p, nel, nnd, normEl, heavyEl, heavytipEl, neighbortipEl, tipEl, geom,...
        connec, nodvalue, x1c, y1c, x2c, y2c, tc, dhx, dhy, h, kapa, miu, Ds, Ainp, As, Bforgauss, displacements, heavynodes, tipnodes, elementDof_inplane, NG2, NG3, KWfound, GPfound, sigma, sigmanGP);
    
    K(elementDof, elementDof) = K(elementDof, elementDof) + Kbend;
    KG(elementDof, elementDof) = KG(elementDof, elementDof) + Kgeo;
    K(elementDof, elementDof) = K(elementDof, elementDof) + Kshear;
end
kkk = K;
ggg = KG;

% Solve buckling eigenproblem
[prescribedDof, activeDof, fixedNodeW] = BC(BCtype, gdof, geom, nodvalue, nnd, heavynodes, tipnodes);
kfactors = solveEigen(K, KG, prescribedDof, w, L, E, E1, E2, v12, v21, h, v, nLambda, gdof, activeDof, connec, geom, nnd, isort, nodvalue, heavynodes, tipnodes, tc, globxyGp, sigmanGP);
