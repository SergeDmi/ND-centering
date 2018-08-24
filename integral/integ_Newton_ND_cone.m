function [MVF,MVL,Fc,Lc]=integ_Newton_ND_cone(R1,PW,D,X,nmt)
%			integ_force_ND_Newton : computes net centering force & distance in a cone
%  Not that a cone is convex and the whole space is visible
% INPUT
%			R1 : radius of cells in doublet : 0.5 (fully separated) to 1 (single sphere)
% 		X  : vector of positions FROM -1 to 0 !!!!!!!!!!!!!!
%			PW : power law for distance (Vector)
%			D  : vector of dimension of the space
%			nmt : number of integration elements per sphere ( = pi / d \theta )
% OUTPUT
%			Fc : mean projected PW-force (on Ox)
%     Lc : mean PW-distance
% DEFINITIONS
%			Mean : average over points on the surface visible from X
% 		p-force : distance^p projected on axis Ox
% 		p-distance : distance^p
%			cone : D-cone of length 1 and radius r1
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
geom=sin(TH).^(D-2);

% First we assume all ends to be at the edge of the base of the cone
L=R1./sin(atan(R1./(1-X)));
endsX=X+cos(TH).*L;
secs=endsX>1;
% All the ends beyond the base
L(secs)=(1-X(secs))./(cos(TH(secs)));
becs=~secs;
% All the rest
L(becs)=R1*X(becs)./(sqrt(1-cos(TH(becs)).^2)-R1*cos(TH(becs)));



for i=1:np
    pw=PW(i);
    Fc(i,:)=sum((L(:,:).^pw).*cos(TH).*geom,1);
    Lc(i,:)=sum((L(:,:).^pw).*geom,1);
end

% We normalize by sum of geom
vol=sqrt(pi)*gamma((D-1)/2.0)/gamma(D/2.0);
Fc=Fc'*dt/vol;
Lc=Lc'*dt/vol;

end
