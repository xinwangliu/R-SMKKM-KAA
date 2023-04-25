clear
clc
warning off;

addpath(genpath('./'));
dataName = 'flower17';
%dataName = 'MF-first3view';
load(['./datasets/',dataName,'_Kmatrix'],'KH','Y');
numclass = length(unique(Y));
Y(Y<1) = numclass;
KH = kcenter(KH);
KH = knorm(KH);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
options.seuildiffsigma=1e-5;        % stopping criterion for weight variation
options.goldensearch_deltmax=1e-1; % initial precision of golden section search
options.numericalprecision=1e-16;   % numerical precision weights below this valueare set to zero
options.firstbasevariable='first'; % tie breaking method for choosing the base variable in the reduced gradient method
options.nbitermax=500;             % maximal number of iteration
options.seuil=0;                   % forcing to zero weights lower than this
options.seuilitermax=10;           % value, for iterations lower than this one
options.miniter=0;                 % minimal number of iterations
options.threshold = 1e-4;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
qnorm = 2;

tic;
tauset12 = 2.^-1;
%tauset12 = 2.^-0;
[H_normalized,Sigma,obj] = R_SMKKM_KAA(KH,numclass,tau,options);
[res_mean(:,1),res_std(:,1)] = myNMIACCV2(H_normalized,Y,numclass);
timecost = toc;   
