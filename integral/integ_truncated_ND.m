function [Fc,Lc]=integ_truncated_ND(H,pw,D,X,nmt)
%			integ_truncated_ND : computes net centering force & distance in a truncated sphere
%  Not that a cone is convex and the whole space is visible
% INPUT
%			R1 : radius of cells in doublet : 0.5 (fully separated) to 1 (single sphere)
% 		X  : vector of positions FROM -1 to 0 !!!!!!!!!!!!!!
%			PW : power law for distance (Vector)
%			D  : vector of dimension of the space
%			nmt : number of integration elements ( = pi / d \theta )
% OUTPUT
%			Fc : mean projected PW-force (on Ox)
%     Lc : mean PW-distance
% DEFINITIONS
%			Mean : average over all orientations from X
% 		p-force : distance^p projected on axis Ox
% 		p-distance : distance^p
%			truncated sphere : D-ball truncated on its second dimension
%
% Serge Dmitrieff,
% Institut Jacques Monod
% www.biophysics.fr

if D<2
	error('Cannot work in d<2')
end


%% Variable intiation
nf=numel(X);
np=numel(H);
dt=pi/nmt;
th=(0:dt:pi)';
nmt=numel(th);
TH=th*ones(1,nf);

Fc=zeros(np,nf);
Lc=zeros(np,nf);
X=ones(nmt,1)*X;
L=zeros(nmt,nf);
Out=false(nmt,nf);

cosTH=cos(TH);
sinTH=sin(TH);
geom=sinTH.^(D-2);

for i=1:np
		% Loop over truncation heights
    h=H(i);
    L(:,:)=-X.*cosTH+sqrt(1-(sinTH.*X).^2);
    Out(:,:)=(L.*sinTH)>h;
    L(Out)=h./sinTH(Out);


    Fc(i,:)=sum((L(:,:).^pw).*cosTH.*geom,1);
    Lc(i,:)=sum((L(:,:).^pw).*geom,1);
end

% Normalization by sum of geom (analytical result)
vol=sqrt(pi)*gamma((D-1)/2.0)/gamma(D/2.0);
Fc=Fc'*dt/vol;
Lc=Lc'*dt/vol;

end
