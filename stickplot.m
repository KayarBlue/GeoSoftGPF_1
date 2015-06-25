function [h] = stickplot(t,u,v,offset);
%function [h] = stickplot(t,u,v,offset);
%
% Makes a stickplot of the field u(t),v(t).
%
% Sticks will begin at each time interval in T, U is the
% velocity in the abscissa-direction, V is the velocity in the
% ordinate-direction
%
% If you desire to have a vertical offset, include a value for
% "OFFSET" (default is 0).
%
% The output is the handle for the stickplot (if you want to change
% its colour or linewidth or something.
%
% e.g. 
%
%  T = [0:2:100];
%  U = 10 * sin(pi * T / 50);
%  V = 10 * cos(pi * T / 50);
%
%  h = stickplot(T,U,V);
%  g = stickplot(T,-U,-V,20);
%
%  set(h,'color','red','linestyle','--');
%  set(g,'color','black','linewidth',2);
%  axis equal

%                                                         rsm || 11Nov03
%--
% Comments added 21Jul04, rsm
%--
 
t = t(:);
t = [t NaN*ones(size(t))]';
t = t(:);

u = u(:);
u = [u NaN*ones(size(u))]';
u = u(:);

v = v(:);
v = [v NaN*ones(size(v))]';
v = v(:);

t = [t t + u]';
t = t(:);
v = [zeros(size(v)) v]';
v = v(:);
if(exist('offset')), v = v + offset; end,
h = line(t,v); 

