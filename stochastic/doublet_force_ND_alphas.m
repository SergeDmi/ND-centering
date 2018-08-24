function [Fc,Lc]=doublet_force_ND_alphas(R1,X,PW,D,nmt)
%			doublet_force_ND_alphas : computes net centering force & distance in a doublet
% INPUT
%			R1 : radius of cells in doublet : 0.5 (fully separated) to 1 (single sphere)
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
%			Doublet : two D-balls of radius r1, truncated and symmetric at x=0
% 					ball centers at r1-1 and 1-r1 respectively
%
% Serge Dmitrieff,
% Institut Jacques Monod
% www.biophysics.fr
%
%% Variable intiation
% Preparation
nf=numel(X);
np=numel(PW);
L0=3.0;
% L is the length of radial elements ; will be corrected later
L=ones(nmt,1);
% mean force & distance
Fc=zeros(nf,np);
Lc=zeros(nf,np);
% Convenience declaration
% End points of radial elements
ends=zeros(nmt,D);
% This happens to be convenient as well
M1=[-1.0+R1 0 0];
M2=-M1;
colM2=ones(nmt,1)*M2;
col1=ones(nmt,1);
colR=zeros(nmt,D);
colRmM1=zeros(nmt,D);
colRmM2=zeros(nmt,D);

%% Integration
for n=1:nf
	% We take random radial elements : much easier in D dimensions...
	Uns=normalize_rows_ND(randsphere(nmt,D,1),D);
	% Unfortunately here we have to take one position after another
	R(1)=X(n);
	% Computing columns of coordinates
	colR(:,1)=X(n);
	colRmM1(:,1)=X(n)-M1(1);
	% Filling L in (avoiding re-allocation), assuming all radial elememnts end in the first sphere
	L(:,1)=-sum(Uns(:,:).*colRmM1(:,:),2)+sqrt(sum(Uns(:,:).*colRmM1(:,:),2).^2+R1*R1-sum(colRmM1.^2,2));
	% computing end positions
	ends(:,:)=colR(:,:)+Uns(:,:).*L(:,1);
	% Finding out which MT go beyond the center of the doublet
	secs=logical(ends(:,1)>0);
	% cutting them short !
	L(secs,1)=(0-colR(secs,1))./Uns(secs,1);
	ends(secs,:)=colR(secs,:)+Uns(secs,:).*L(secs,1);
	% Now is a bit more tricky : we compute the length from mid-plane to the second sphere of the doublet
	colRmM2(secs,:)=ends(secs,:)-colM2(secs,:);
	L(secs,1)=L(secs,1)-sum(Uns(secs,:).*colRmM2(secs,:),2)+sqrt(sum(Uns(secs,:).*colRmM2(secs,:),2).^2+R1*R1-sum(colRmM2(secs,:).^2,2));

	% Now computing the quantities. Wasn't it easy ?
	for i=1:np
		pw=PW(i);
		Fc(n,i)=sum((L(:,1).^pw).*Uns(:,1),1);
		Lc(n,i)=sum((L(:,1).^pw));
	end
end

Fc(:,:)=Fc/nmt;
Lc(:,:)=Lc/nmt;
end
