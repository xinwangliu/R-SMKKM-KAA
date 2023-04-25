function [cost,Hstar] = cost_R_SMKKM_KAA(KH,StepSigma,DirSigma,Sigma,Hstar,H0,lambda,numclass)

global nbcall
nbcall=nbcall+1;

Sigma = Sigma+ StepSigma * DirSigma;

Kmatrix = sumKbeta(KH,(Sigma.*Sigma));
[Hstar,cost]= mykernelkmeans(Kmatrix+lambda*(H0*H0'),numclass);