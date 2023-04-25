function [Hstar,Sigma,obj] = R_SMKKM_KAA(KH,numclass,lambda,option)

numker = size(KH,3);
Sigma = ones(numker,1)/numker;

%--------------------------------------------------------------------------------
% Options used in subroutines
%--------------------------------------------------------------------------------
if ~isfield(option,'goldensearch_deltmax')
    option.goldensearch_deltmax=5e-2;
end
if ~isfield(option,'goldensearchmax')
    optiongoldensearchmax=1e-8;
end
if ~isfield(option,'firstbasevariable')
    option.firstbasevariable='first';
end

nloop = 1;
loop = 1;
goldensearch_deltmaxinit = option.goldensearch_deltmax;
%-----------------------------------------
% Initializing H0
%------------------------------------------
Kmatrix = sumKbeta(KH,Sigma.^2);
[Hstar]= mykernelkmeans(Kmatrix,numclass);
H0 = Hstar;

[Hstar,obj1] = mySolvingHIter(Kmatrix,Hstar,H0,lambda);
obj(nloop) = obj1(end);

[grad] = Grad_R_SMKKM_KAA(KH,Hstar,Sigma);
Sigmaold  = Sigma;

%------------------------------------------------------------------------------%
% Update Main loop
%------------------------------------------------------------------------------%

while loop
    nloop = nloop+1;
    %-----------------------------------------
    % Update weigths Sigma
    %-----------------------------------------
    [Sigma,Hstar,obj(nloop)] = update_R_SMKKM_KAA(KH,Sigmaold,grad,obj(nloop-1),Hstar,H0,lambda,option);
  

    %-----------------------------------------------------------
    % Enhance accuracy of line search if necessary
    %-----------------------------------------------------------
    if max(abs(Sigma-Sigmaold))<option.numericalprecision &&...
            option.goldensearch_deltmax > optiongoldensearchmax
        option.goldensearch_deltmax=option.goldensearch_deltmax/10;
    elseif option.goldensearch_deltmax~=goldensearch_deltmaxinit
        option.goldensearch_deltmax*10;
    end
    
    [grad] = Grad_R_SMKKM_KAA(KH,Hstar,Sigma);
    %----------------------------------------------------
    % check variation of Sigma conditions
    %----------------------------------------------------
    if  max(abs(Sigma-Sigmaold))<option.seuildiffsigma
        loop = 0;
        fprintf(1,'variation convergence criteria reached \n');
    end
    %-----------------------------------------------------
    % Updating Variables
    %----------------------------------------------------
    Sigmaold  = Sigma;
end