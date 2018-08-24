function [Fc,Ac]=integ_doublet_ND_Newton(R1,PW,D,X,nmt)
%			integ_doublet_ND_Newton : computes net centering force & distance in a doublet
% INPUT
%			R1 : radius of cells in doublet : 0.5 (fully separated) to 1 (single sphere)
% 		X  : vector of positions FROM -1 to 0 !!!!!!!!!!!!!!
%			PW : power law for distance
%			D  : vector of dimension of the space
%			nmt : number of integration elements per sphere ( = pi / d \theta )
% OUTPUT
%			Fc : mean projected PW-force (on Ox)
%     Ac : Surface visible from X
% DEFINITIONS
%			Mean : average over points on the surface visible from X
% 		p-force : distance^p projected on axis Ox
% 		p-distance : distance^p
%			Doublet : two D-balls of radius r1, truncated and symmetric at x=0
% 					ball centers at r1-1 and 1-r1 respectively
%
% Serge Dmitrieff,
% Institut Jacques Monod
% www.biophysics.fr
%

if max(X)>0
	error('Does not support X>0')
end

%% Variable intiation
nf=numel(X);
np=numel(R1);
dt=pi/nmt;
th=(0:dt:pi)';
nmt=numel(th);
TH=th*ones(1,nf);
Fc=zeros(np,nf);
Ac=zeros(np,nf);
X=ones(nmt,1)*X;
Kept1=false(nmt,nf);
Kept2=false(nmt,nf);
L1=zeros(nmt,nf);
L2=zeros(nmt,nf);
geom=sinTH.^(D-2);

% And sometimes, one can be cheap
cosTH=cos(TH);
sinTH=sin(TH);
FF=zeros(nmt,nf);
AA=zeros(nmt,nf);

% Integration
for i=1:np
  r1=R1(i);

	% Slightly cheaper...
	FF(:)=0;
	AA(:)=0;
  L1(:,:)=sqrt(r1*r1+(r1-1-X).^2-2*r1*(r1-1-X).*cosTH);
  L2(:,:)=sqrt(r1*r1+(X-1+r1).^2-2*r1*(X-1+r1).*cosTH);

	% Keeping only visible points on both spheres...
  Kept1(:,:)=cosTH<((1-r1)/r1);
  Kept2(:,:)=(cosTH>((r1-1)/r1)).*(cosTH>((r1-1)+X-L2(:,:).*X./sqrt(X.^2+2*r1-1))/r1);

	% Now doing the grunt work
	FF(Kept1)=(L1(Kept1).^(PW-1)).*geom(Kept1).*(r1-1+r1*cosTH(Kept1)-X(Kept1));
	FF(Kept2)=FF(Kept2)+(L2(Kept2).^(PW-1)).*geom(Kept2).*(1-r1+r1*cosTH(Kept2)-X(Kept2));
	Fc(i,:)=sum(FF,1);

	AA(Kept1)=geom(Kept1);
	AA(Kept2)=AA(Kept2)+geom(Kept2);
  Ac(i,:)=sum(AA,1);
end


% Oh, wait, we're computing an average !
Fc=(Fc./Ac)';


end
