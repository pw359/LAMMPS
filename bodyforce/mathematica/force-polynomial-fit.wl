(* ::Package:: *)

(* ::Text:: *)
(*We want to fit a 3rd-order polynomial to a function f(x) such that function values and derivatives are continuous at x = x0+- \[Delta]/2. We first consider the case x0 = 0 and then translate the solution for the more general case.*)


ClearAll[p, \[Delta], x,f,ptmp, pfit,a0,a1,a2,a3]
p[x_]= a0 + a1*x + a2*x^2 + a3*x^3;
ptmp =  p[x] /. Solve[{   p[\[Delta]/2] == f[\[Delta]/2], p'[\[Delta]/2] == f'[\[Delta]/2], p[-\[Delta]/2] == f[-\[Delta]/2] , p'[-\[Delta]/2] == f'[-\[Delta]/2]}, {a0,a1,a2,a3}]  [[1]]
pfit[y_] := ptmp /. x-> y
pfit[y] // TraditionalForm
g[z_] := If[ Abs[z] <= \[Delta]/2, pfit[z ], f[z]]


(* ::Text:: *)
(*Let's see if it works for x=0*)


ClearAll[f]
\[Delta] = 0.4;
f[x_] := x/Sqrt[x*x]+ x;
Plot[pfit[z],{z,-0.7,0.7}]
Plot[ {pfit[z], f[z]},{z,-0.7,0.7}]
Plot[g[z],{z,\[Delta]/2 - 0.1, \[Delta]/2 + 0.1}]
Plot[D[g[zz],zz] /.zz->z,{z,-1,1}, PlotRange->Full]


(* ::Text:: *)
(*Now let's use the previous result to fit a polynomial centred around an arbitrary point x0. The function polyFit below performs the polynomial fit in the region of interest and returns the original function values outside that region.*)


polyFit[y_, y0_, \[Delta]_,f_] := Module[ {x = y-y0,  fl = f[y0-\[Delta]/2], fr=f[y0+\[Delta]/2], dfl=f'[y0-\[Delta]/2], dfr=f'[y0+\[Delta]/2]},  
If [Abs[x] <= \[Delta]/2, \!\(TraditionalForm\`
\*FractionBox[\(1\), \(8\)]\ \((\[Delta]\ dfl - \[Delta]\ dfr + 4\ fl + 4\ fr)\) - 
\*FractionBox[\(
\*SuperscriptBox[\(x\), \(3\)]\ \((\(-\[Delta]\)\ dfl - \[Delta]\ dfr - 2\ fl + 2\ fr)\)\), 
SuperscriptBox[\(\[Delta]\), \(3\)]] - 
\*FractionBox[\(
\*SuperscriptBox[\(x\), \(2\)]\ \((dfl - dfr)\)\), \(2\ \[Delta]\)] - 
\*FractionBox[\(x\ \((\[Delta]\ dfl + \[Delta]\ dfr + 6\ fl - 6\ fr)\)\), \(4\ \[Delta]\)]\),f[y]]]


(* ::Text:: *)
(*Let's see if it works for x = x0*)


Plot[polyFit[x,0, \[Delta], f], {x, -1,1}]


x0 = 1;
g[x_] := f[x - x0];
Plot[polyFit[x,x0, 3*\[Delta], g], {x, x0-1,x0+1}]
