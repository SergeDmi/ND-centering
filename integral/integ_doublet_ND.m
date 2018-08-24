function [Fc,Lc]=integ_doublet_ND(R1,PW,D,X,nmt)
%			doublet_force_ND_alphas : computes net centering force & distance in a doublet
% INPUT
%			R1 : radius of cells in doublet : 0.5 (fully separated) to 1 (single sphere)
% 		X  : vector of positions FROM -1 to 0 !!!!!!!!!!!!!!
%			PW : power law for distance
%			D  : dimension of the space
%			nmt : number of integration elements ( = pi / d \theta )
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
%ysical constants

if D<2
	error('Cannot work in d<2')
end


%% Variable intiation
nf=numel(X);
np=numel(R1);
dt=pi/nmt;
th=(0:dt:pi)';
nmt=numel(th);
TH=th*ones(1,nf);
Fc=zeros(np,nf);
Lc=zeros(np,nf);
MVF=zeros(1,np);
MVL=zeros(1,np);

X=ones(nmt,1)*X;
L=zeros(nmt,nf);
Out=false(nmt,nf);

cosTH=cos(TH);
sinTH=sin(TH);

geom=sinTH.^(D-2);
Yend=zeros(nmt,nf);
for i=1:np

    r1=R1(i);
    L(:,:)=-(X+1-r1).*cosTH+sqrt(r1*r1-(sinTH.*(X+1-r1)).^2);

		% Finding all MT going to the second sphere and recomputing the length
    Out(:,:)=(X+L.*cosTH)>0;
    L(Out)=-X(Out)./cosTH(Out);
    Yend(Out)=L(Out).*sinTH(Out);
    L(Out)=L(Out)-(r1-1).*cosTH(Out)-Yend(Out).*sinTH(Out)+sqrt(r1*r1-((r1-1)^2).*(sinTH(Out).^2)-(Yend(Out).*cosTH(Out)).^2+2.0*(r1-1)*Yend(Out).*cosTH(Out).*sinTH(Out));

		% Yes it's done
    Fc(i,:)=sum((L(:,:).^PW).*cosTH.*geom,1);
    Lc(i,:)=sum((L(:,:).^PW).*geom,1);
end

% Normalizing by the sum of geom
vol=sqrt(pi)*gamma((D-1)/2.0)/gamma(D/2.0);
Fc=Fc'*dt/vol;
Lc=Lc'*dt/vol;

end
