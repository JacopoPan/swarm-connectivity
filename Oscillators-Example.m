#! /usr/bin/octave
#Ubuntu's interpreter /usr/bin/octave
#OSX's interpreter /usr/local/bin/octave

addpath("~/Dropbox/papers/2016/saso-carlo/");


x = 0.15;

a2 = 1.3;

for iter=1:100
	dx = (-1.0*a2*sign(x)*(1 - x))-x;
	
	x = x + dx
	
	input("");
endfor

input("");

x = 1;
y = 1;
z = 0.5;

sig = 10;
rho = 28;
bet = 8/3;

step = 100

for iter=1:1000
	dx = sig*(y - x);
	dy = x*(rho - z) - y;
	dz = x*y - bet*z;
	
	x = x + dx/step
	y = y + dy/step
	z = z + dz/step
	
	input("");
endfor
