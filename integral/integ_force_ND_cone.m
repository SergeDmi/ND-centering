function [Fc,Lc]=integ_force_ND_cone(R1,PW,D,X,nmt)
%			cone_force_ND_alphas : computes net centering force & distance in a cone
% INPUT
%			R1 : final radius of cone
% 		X  : vector of positions FROM 0 to 1 !!!!!!!!!!!
%			PW : power law for distance
%			D  : dimension of the space
%			nmt : number of radial elements (pi / d \theta)
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
dt=pi/nmt;
th=(0:dt:pi)';
nmt=numel(th);
TH=th*ones(1,nf);
np=numel(PW);
Fc=zeros(np,nf);
Lc=zeros(np,nf);

% Bits of tricks to compute all in one baduk... in one go.
X=ones(nmt,1)*X;
geom=sin(TH).^(D-2);

% First we assume all elements as long as to reach the edge of the base
L=R1./sin(atan(R1./(1-X)));
endsX=X+cos(TH).*L;
% Finding out which MT go beyond the base of the cone
secs=endsX>1;
L(secs)=(1-X(secs))./(cos(TH(secs)));
becs=~secs;
% And the rest !
L(becs)=R1*X(becs)./(sqrt(1-cos(TH(becs)).^2)-R1*cos(TH(becs)));

% Fast !
for i=1:np
    pw=PW(i);
    Fc(i,:)=sum((L(:,:).^pw).*cos(TH).*geom,1);
    Lc(i,:)=sum((L(:,:).^pw).*geom,1);
end

% Normalizing by sum of geom
vol=sqrt(pi)*gamma((D-1)/2.0)/gamma(D/2.0);
Fc=Fc'*dt/vol;
Lc=Lc'*dt/vol;

end
