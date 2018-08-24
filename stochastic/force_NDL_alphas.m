function [MVF,MVL,Fc,Lc,X]=force_NDL_alphas(X,PW,D,nmt)
%			doublet_NDL_alphas : computes net centering force & distance in a sphere
% INPUT
% 		X  : vector of positions FROM -1 to 0 !!!!!!!!!!!!!!
%			PW : power law for distance
%			D  : dimension of the space
%			nmt : number of radial elements (microtubules)
% OUTPUT
%			Fc : mean projected PW-force (on Ox)
%     Lc : mean PW-distance
% DEFINITIONS
%			Mean : average over all angles
% 		p-force : distance^p projected on axis Ox
% 		p-distance : distance^p
%			Sphere : D-Ball of radius 1
%
% Serge Dmitrieff,
% Institut Jacques Monod
% www.biophysics.fr
%
%% Variable intiation
nf=numel(X);
np=numel(PW);
Fc=zeros(nf,np);
Lc=zeros(nf,np);
L0=3.0;
L=ones(nmt,1)*L0;
Rad2=1.0;
colR=zeros(nmt,D);
R=zeros(1,D);
%% Integration
for n=1:nf
	% Using randomly chosen orientation vectors
	Uns=normalize_rows_ND(randsphere(nmt,D,1),D);
	% We have to take positions one at a time
	R(1)=X(n);
	colR(:,1)=X(n);
	% Avoids variable reinitialization
	L(:,1)=-sum(colR(:,:).*Uns(:,:),2)+sqrt(Rad2-sum(R.^2)+sum(colR(:,:).*Uns(:,:),2).^2);
	% that was easy
	for i=1:np
		pw=PW(i);
		Fc(n,i)=sum((L(:,1).^pw).*Uns(:,1),1);
		Lc(n,i)=sum((L(:,1).^pw));
	end
end

Fc(:,:)=Fc/nmt;
Lc(:,:)=Lc/nmt;

end
