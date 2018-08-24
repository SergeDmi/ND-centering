function [Fc,Lc,Uns,L]=cone_force_ND_alphas(R1,X,PW,D,nmt)
%			cone_force_ND_alphas : computes net centering force & distance in a cone
% INPUT
%			R1 : final radius of cone
% 		X  : vector of positions FROM 0 to 1 !!!!!!!!!!!
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
%			cone : D-cone of length 1 and final radius R1
%
% Serge Dmitrieff,
% Institut Jacques Monod
% www.biophysics.fr


%% Variable intiation
nf=numel(X);
np=numel(PW);
Fc=zeros(nf,np);
Lc=zeros(nf,np);
L=ones(nmt,1);
colR=zeros(nmt,D);
ends=zeros(nmt,D);
secs=false(nmt,1);
becs=false(nmt,1);

%% Integration
for n=1:nf
	% Orientation vectors are randomly chosen
	Uns=normalize_rows_ND(randsphere(nmt,D,1),D);

	%L=ones(np,N)*L0;
	x=X(n);
	%R(1)=X(n);
	colR(:,1)=X(n);
	% firs we assume all length as long as to reach the corner at the base
	L(:,1)=R1/sin(atan(R1/(1-x)));
	ends(:,:)=colR(:,:)+Uns(:,:).*L(:,1);
	% Finding out which MT go beyond the cone base
	secs(:)=ends(:,1)>1;
	% cutting them short !
	L(secs,1)=(1-colR(secs,1))./Uns(secs,1);
	% The other ones : making sure they don't stick out of the cone !
	becs(:)=~secs;
	L(becs,1)=R1*x./(sqrt(1-Uns(becs,1).^2)-R1*Uns(becs,1));

	for i=1:np
		pw=PW(i);
		Fc(n,i)=sum((L(:,1).^pw).*Uns(:,1),1);
		Lc(n,i)=sum((L(:,1).^pw));
	end
end

Fc(:,:)=Fc/nmt;
Lc(:,:)=Lc/nmt;

end
