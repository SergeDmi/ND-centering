function [Fc,Lc]=integ_force_ND(PW,D,X,nmt)
%			integ_force_ND : computes net centering force & distance in a sphere
% INPUT
%			R1 : radius of cells in doublet : 0.5 (fully separated) to 1 (single sphere)
% 		X  : vector of positions FROM -1 to 0 !!!!!!!!!!!!!!
%			PW : power law for distance (Vector)
%			D  : vector of dimension of the space
%			nmt : number of integration elements  ( = pi / d \theta )
% OUTPUT
%			Fc : mean projected PW-force (on Ox)
%     Lc : mean PW-distance
% DEFINITIONS
%			Mean : average over orientations from X
% 		p-force : distance^p projected on axis Ox
% 		p-distance : distance^p
%			sphere : D-balls of radius r1
%
% Serge Dmitrieff,
% Institut Jacques Monod
% www.biophysics.fr
%
%% Variable intiation
nf=numel(X);
dt=pi/nmt;
th=(0:dt:pi)';
nmt=numel(th);
TH=th*ones(1,nf);
np=numel(PW);
Fc=zeros(np,nf);
Lc=zeros(np,nf);

X=ones(nmt,1)*X;
cosTH=cos(TH);
L=-X.*cosTH+sqrt(1-(X.*sin(TH)).^2);
geom=sin(TH).^(D-2);


for i=1:np
    pw=PW(i);
    Fc(i,:)=sum((L(:,:).^pw).*cosTH.*geom,1);
    Lc(i,:)=sum((L(:,:).^pw).*geom,1);
end

% We normalize by the sum of geom (analytical result)
vol=sqrt(pi)*gamma((D-1)/2.0)/gamma(D/2.0);
Fc=Fc'*dt/vol;
Lc=Lc'*dt/vol;

end
