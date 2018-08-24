function [Fc,Lc]=force_ND_truncated(H,X,PW,d,nmt)
%			force_ND_truncated : computes net centering force & distance in a truncated sphere
% INPUT
%			H : height of trunction ; from 0 (disc) to 1 : sphere
% 		X  : vector of positions FROM -1 to 0 !!!!!!!!!!!!!!
%			PW : power law for distance
%			d  : dimension of the space
%			nmt : number of radial elements (microtubules)
% OUTPUT
%			Fc : mean projected PW-force (on Ox)
%     Lc : mean PW-distance
% DEFINITIONS
%			Mean : average over all angles
% 		p-force : distance^p projected on axis Ox
% 		p-distance : distance^p
%			Truncated sphere : D-Ball of radius 1 truncated along y (second dimension)
%
% Serge Dmitrieff,
% Institut Jacques Monod
% www.biophysics.fr

if d<2
	error('Cannot work in d<2')
end
%% Variable intiation
nf=numel(X);
np=numel(H);
Fc=zeros(nf,np);
Lc=zeros(nf,np);
LL=ones(nmt,1);
L=ones(nmt,1);
colR=zeros(nmt,d);
R=zeros(1,d);
out=false(nmt,1);
%% Integration
for n=1:nf
	% We chose random orientation vectors
	Uns=normalize_rows_ND(randsphere(nmt,d,1),d);
	% We have to iterate over position x
	R(1)=X(n);
	colR(:,1)=X(n);

	% We compute the size of radial elements in the non-truncated sphere
	LL(:,1)=(-sum(colR(:,:).*Uns(:,:),2)+sqrt(Rad2-sum(R.^2)+sum(colR(:,:).*Uns(:,:),2).^2));

	% Then we loop over truncation heights
	% Ok there is a non-loop way to do this, I know !
	for i=1:np
		h=H(i);
		L(:,1)=LL(:,1);
		out(:)=logical(abs(colR(:,2)+Uns(:,2).*L(:,1))>h);
		L(out,1)=(h./abs(Uns(out,2)));
		Fc(n,i)=sum((L(:,1).^pw).*Uns(:,1),1);
		Lc(n,i)=sum((L(:,1).^pw));
	end

end

Fc(:,:)=Fc/nmt;
Lc(:,:)=Lc/nmt;

end
