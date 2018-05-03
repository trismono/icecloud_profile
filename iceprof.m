% A function to derive profile of cirrus clouds
% trismono_candra.krisna@uni-leipzig.de
% version 26 Feb. 2018
% ============================================ %
% to run, you need input parameters
% tc = total cloud optical thickness 
% rt = particle effective radius at cloud top (µm)
% rb = particle effective radius at cloud base (µm)
% zt = cloud altitude at cloud top (m)
% zb = cloud altitude at cloud base (m)
% n_layer = number of layer in the cloud (20 or 30 is an ideal number)
% k = parameter influencing the curvature of the profile (typically 2-4)
% 
% Example :
% [dre,dzl,dtau] = prof(tc,rb,rt,zb,zt,n_layer,k)
%                    or
% [dre,dzl,dtau] = prof(2,40,10,10000,12000,20,3)
% this will derive profile of particle effective radius as a function
% of optical thickness and geometrical altitude
% ============================================= %

function [prof,tau] = iceprof(tc,rb,rt,iwcb,iwct,zb,zt,n_layer,k)

% basic parameters
n_layer = n_layer+1         ; % n-layer
h = abs(zt-zb)              ; % geometrical thickness 
z = linspace(0,h,n_layer)   ; % altitude of each cloud layers

a0 = rt+rb          ;
a1 = rt^k           ;
a2 = rt^k-rb^k      ;
re_dummy = a0 - (a1 - a2.*z./h).^(1/k) ;

% linear decrease as the mirror
re_linear = linspace(rb,rt,n_layer)    ;

% line equation
p1 = [rb,0]         ;
p2 = [rt,h]         ;
m = (p2(2)-p1(2))/(p2(1)-p1(1)) ;
c = p1(2)-m*p1(1)   ;
% sudut rotasi
alpha = atand(m)    ;

% Transformasi / percerminan terhadap garis y = mx+c
for i = 1:numel(re_dummy)
    var(:,i) = [cosd(2*alpha) sind(2*alpha);sind(2*alpha) -cosd(2*alpha)]*...
        [re_dummy(i);z(i)-c]+[0;c]  ;
end

% output parameters
re = (var(1,:)')                    ; % particle effective radius
zl = (var(2,:)+zb)'./1000           ; % altitude of effective radius (km)
iwc = linspace(iwcb,iwct,n_layer)'  ; % ice water content (g/m^3)
prof = flipud([zl,iwc,re])          ; % sorting in order from cloud top

% tau profile
tau = (linspace(tc,0,n_layer))'     ; % homogeneous thin optical thickness
tau = flipud(tau)                   ; % tau starts from cloud base

end