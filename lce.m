% This m-file generates the standard colormap used for centuries by students
% of Elliot McVeigh...
%
% DBE 4/5/01

function lce_color=lce(m);

if nargin==0;
  m=64;
end

s=linspace(0,1,m);

for i=1:length(s)
	if s(i)<0.5
      red(i)=2*s(i);
      green(i)=0;
      blue(i)=1-red(i);
	else
      red(i)=1;
      green(i)=2*s(i)-1;
      blue(i)=0;
	end
end

lce_color=[red; green; blue;]';