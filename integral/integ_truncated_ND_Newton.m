function [Fc,Ac]=integ_truncated_ND_Newton(H,PW,D,X,ds)
%			integ_truncated_ND_Newton : computes net centering force & distance in a truncated sphere
%  Not that a cone is convex and the whole space is visible
% INPUT
%			R1 : radius of cells in doublet : 0.5 (fully separated) to 1 (single sphere)
% 		X  : vector of positions FROM -1 to 0 !!!!!!!!!!!!!!
%			PW : power law for distance (Vector)
%			D  : vector of dimension of the space
%			nmt : number of integration elements per sphere ( = pi / d \theta )
%						ok, it's actually more complicated that that, the number depends on
%						the truncation height
% OUTPUT
%			Fc : mean projected PW-force (on Ox)
%     Lc : mean PW-distance
% DEFINITIONS
%			Mean : average over points on the surface visible from X
% 		p-force : distance^p projected on axis Ox
% 		p-distance : distance^p
%			truncated sphere : D-ball truncated on its second dimension
%
% Serge Dmitrieff,
% Institut Jacques Monod
% www.biophysics.fr

if max(X)>0
	error('Does not support X>0')
end

%% Variable intiation
nf=numel(X);
np=numel(H);
Fc=zeros(np,nf);
Ac=zeros(np,nf);
row1=ones(1,nf);

%% 3,2,1, let's jam !
% Ok this is actually a bit tricky
% We need to take regularly spaced points on the truncated sphere
% We do that and we try to keep a constant distance even at corners...
for i=1:np
  	% For truncation height h
    h=H(i);
		% First we deal with the geometry, drawing the points on the truncated 2-ball
    th1=asin(h);
    th2=pi-th1;
    angs1=0:ds:th1;
    xend=cos(angs1(end));
    xx=(xend-ds):(-ds):(-xend);
    nx=numel(xx);
    angs2=(acos(xx(end))+ds):ds:pi;

		% Now we get the angles from the intermediate points (on the truncation)
    sin_mid=abs(h./xx)./sqrt(1+(h./xx).^2);
    cos_mid=sign(xx)./sqrt(1+(h./xx).^2);

		% we gather all the points
    XP=([cos(angs1) xx cos(angs2)]')*row1;
    YP=([sin(angs1) h*ones(1,nx) sin(angs2)]')*row1;

		% We compute the sines for the projection
    sinP=([sin(angs1) sin_mid sin(angs2)]')*row1;
    % Not needed : cosP=([sin(angs1) cos_mid sin(angs2)]')*row1;

		% The number of points obviously depends on truncation height
    nmt=size(XP,1);
    XX=ones(nmt,1)*X;

		% Ok now it's actually super straightforward
		% Because we have the coordinates of the points !
    L=sqrt(((XP-XX).^2)+YP.^2);
    Fc(i,:)=sum((L.^(PW-1)).*(XP-XX).*sinP.^(D-2),1);
    Ac(i,:)=sum(sinP.^(D-2),1);

end

% Normalization by the number of points and geometry
Fc=(Fc./Ac)';

end
