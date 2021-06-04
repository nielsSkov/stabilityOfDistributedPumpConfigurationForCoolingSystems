% clear all, close all; clc

%% Load Model Parameters

run('modelParameters')

%% Hydrolic Circuit

Lambda = zeros(4,4,4);
Psi    = zeros(4,4,4);
Gamma  = zeros(4,4,4);
S      = zeros(4,4,4);
%
for i = 1:4
	
	[ G2, G3, G4 ] = deal( zeros(4,4) );
	
	G1          = ones(4,4);
	G2(2:4,2:4) = ones(3,3);
	G3(3:4,3:4) = ones(2,2);
	G4(4,4)     = 1;
	
	Lambda(:,:,i) = ( r(i) + a(i) )/b(i)*diag([(i==1) (i==2) (i==3) (i==4)]);
	
	Psi(:,:,i)    = R_c/b(i)*ones(4,4);
	
	Gamma(:,:,i) =  2*R(1)/b(i)*G1*(i>0) + ...
									2*R(2)/b(i)*G2*(i>1) + ...
									2*R(3)/b(i)*G3*(i>2) + ...
									2*R(4)/b(i)*G4*(i>3);
	
	S(:,:,i) = Lambda(:,:,i) + Psi(:,:,i) + Gamma(:,:,i);
end