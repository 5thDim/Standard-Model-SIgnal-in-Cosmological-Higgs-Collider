(* Reproduce the P and soft-loop reductions from the byte-identical original notebook.
   Runs in a fresh Wolfram 15 kernel. Outputs go to /tmp via WSL UNC paths.
   This is a derivation/audit script; the optimized notebook does not import it.
   Gauss reduction is applied only to parameter derivatives with a surviving
   untouched equal upper/lower parameter. No regulator derivative is dropped.
*)
tmp="\\\\wsl.localhost\\Ubuntu\\tmp\\";
root=DirectoryName[DirectoryName[$InputFileName]];
nb=Get[FileNameJoin[{root,"checks","Gauge Boson Signal.before-optimization.nb"}]];
inputs=Cases[nb,Cell[BoxData[b_],"Input",___]:>ToExpression[b,StandardForm,HoldComplete],Infinity];
Scan[ReleaseHold,inputs[[{32,33}]]];
old={P2LSimp[nu],P2PMSimp[nu]};
canon[expr_]:=expr/.{HoldPattern[Derivative[aa_List,bb_List,0][HypergeometricPFQRegularized][a_List,b_List,1]]:>Module[{o=Ordering[Transpose[{b,bb}]]},Derivative[aa,bb[[o]],0][HypergeometricPFQRegularized][a,b[[o]],1]],HoldPattern[HypergeometricPFQRegularized[a_List,b_List,1]]:>HypergeometricPFQRegularized[a,Sort[b],1]};
p=canon[old];
masters=DeleteDuplicates[Cases[p,(Derivative[_,_,0][HypergeometricPFQRegularized][_,_,1])|_HypergeometricPFQRegularized|_HypergeometricPFQ,Infinity]];
gaussDerivative[oa_,ob_,a_,b_]:=Module[{zz=Array[z,5],ap,bp,ia,ib,pairs,gamma,der},
 ap=a+Take[zz,3] Map[If[#==0,0,1]&,oa];bp=b+Take[zz,-2] Map[If[#==0,0,1]&,ob];
 pairs=Select[Tuples[{Range[3],Range[2]}],SameQ[ap[[#[[1]]]],bp[[#[[2]]]]]&];
 If[pairs==={},Return[Derivative[oa,ob,0][HypergeometricPFQRegularized][a,b,1]]];
 {ia,ib}=First[pairs];
 gamma=With[{u=Delete[ap,ia],v=First[Delete[bp,ib]]},Gamma[v] Gamma[v-Total[u]]/(Times@@Gamma[v-u] Times@@Gamma[bp])];
 der=Fold[D[#1,{#2[[1]],#2[[2]]}]&,gamma,Transpose[{zz,Join[oa,ob]}]]/.Thread[zz->0];
 der];
rules=Map[#->If[FreeQ[#,Derivative],TimeConstrained[FunctionExpand[#,nu>0],3,#],#/.HoldPattern[Derivative[oa_List,ob_List,0][HypergeometricPFQRegularized][a_List,b_List,1]]:>gaussDerivative[oa,ob,a,b]]&,masters];
p=p/.rules;
remaining=DeleteDuplicates[Cases[p,(Derivative[_,_,0][HypergeometricPFQRegularized][_,_,1])|_HypergeometricPFQRegularized|_HypergeometricPFQ,Infinity]];
Put[remaining,tmp<>"gauge_final_masters.wl"];
ss=Array[h,Length[remaining]]; pp=Expand[p/.Thread[remaining->ss]];
cc=Table[Prepend[Table[Coefficient[e,v],{v,ss}],e/.Thread[ss->0]],{e,pp}];
Put[{cc,remaining},tmp<>"gauge_p_before_normalize.wl"];

tmp="\\\\wsl.localhost\\Ubuntu\\tmp\\";
{cc,mm}=Get[tmp<>"gauge_p_before_normalize.wl"];
log={};out[x_]:=(AppendTo[log,ToString[x,InputForm]];Export[tmp<>"gauge_normalize_report.txt",StringRiffle[log,"\n"],"Text"]);
rf[z_,n_Integer]:=If[n>=0,Product[z+j,{j,0,n-1}],1/Product[z+j,{j,n,-1}]];
pgshift[k_,z_,n_Integer]:=PolyGamma[k,z]+(-1)^k k! If[n>=0,Sum[1/(z+j)^(k+1),{j,0,n-1}],-Sum[1/(z+j)^(k+1),{j,n,-1}]];
ee=cc/.nu->-I u;
ee=ee/.{Gamma[a_]:>With[{c=a/.u->0},If[IntegerQ[2c],With[{n=Floor[c]},Gamma[a-n]rf[a-n,n]],Gamma[a]]],PolyGamma[k_Integer,a_]:>With[{c=a/.u->0},If[IntegerQ[2c],With[{n=Floor[c]},pgshift[k,a-n,n]],PolyGamma[k,a]]]};
out[{"Gamma bases",DeleteDuplicates[Cases[ee,_Gamma,Infinity]]}];
out[{"PolyGamma bases",DeleteDuplicates[Cases[ee,_PolyGamma,Infinity]]}];
ee=ee/.{Gamma[-2u]->-Pi/(2u Sin[2Pi u]Gamma[2u]),Gamma[-u]->-Pi/(u Sin[Pi u]Gamma[u]),Gamma[1/2-u]->Pi/(Cos[Pi u]Gamma[1/2+u]),PolyGamma[0,-2u]->PolyGamma[0,2u]+Pi Cot[2Pi u]+1/(2u),PolyGamma[0,-u]->PolyGamma[0,u]+Pi Cot[Pi u]+1/u,PolyGamma[0,1/2-u]->PolyGamma[0,1/2+u]-Pi Tan[Pi u],PolyGamma[1,-2u]->Pi^2 Csc[2Pi u]^2-PolyGamma[1,2u]+1/(4u^2),PolyGamma[1,-u]->Pi^2 Csc[Pi u]^2-PolyGamma[1,u]+1/u^2,PolyGamma[1,1/2-u]->Pi^2 Sec[Pi u]^2-PolyGamma[1,1/2+u]};
ee=ee/.{Gamma[u]->2^(1-2u)Sqrt[Pi]Gamma[2u]/Gamma[1/2+u],PolyGamma[0,u]->2PolyGamma[0,2u]-2Log[2]-PolyGamma[0,1/2+u],PolyGamma[1,u]->4PolyGamma[1,2u]-PolyGamma[1,1/2+u]};
ee=ee/.{Gamma[2u]->gg,Gamma[1/2+u]->hh,PolyGamma[0,2u]->ps,PolyGamma[0,1/2+u]->qh,PolyGamma[1,2u]->ps1,PolyGamma[1,1/2+u]->qh1};
ee=TrigToExp[ee];
ee=ee/.Power[E,a_]:>w^Simplify[a/(-I Pi u)];
ee=Map[Factor,ee,{2}];
out[{"leaves",LeafCount/@ee}];
Put[{ee,mm},tmp<>"gauge_p_factored.wl"];
out[ee];

tmp="\\\\wsl.localhost\\Ubuntu\\tmp\\"; 
log={};out[x_]:=(AppendTo[log,ToString[x,InputForm]];Export[tmp<>"gauge_soft_report.txt",StringRiffle[log,"\n"],"Text"]);
Scan[ReleaseHold,Flatten[inputs[[{4,5,6,12,13,21,45,58,60,63,64,67}]]]];
raw=Block[{f0},f0[a_,b_]:=fSym[a,b];Expand[{bubblePolarizationSum[nu],IL[nu,t1],IP[nu,t1],SLL[nu,t1,t3,ph],SHH[nu,t1,t3,ph],SLH[nu,t1,t3,ph]}]];
seeds=DeleteDuplicates[Cases[raw,_fSym,Infinity]];
out[{"seeds",seeds}];
shift[j_,k_,x_]:=FunctionExpand[Pochhammer[2 x-3/2,j+k] Pochhammer[3/2-x,-j] Pochhammer[3/2-x,-k]/(Pochhammer[3-2 x,-j-k] Pochhammer[x,j] Pochhammer[x,k])];
rules=Map[#->shift[Simplify[#[[1]]+I nu],Simplify[#[[2]]+I nu],-I nu]&,seeds];
ratios=Map[FullSimplify[#,nu>0 && Element[{t1,t3,ph},Reals],TimeConstraint->8]&,raw/.rules];
out[{"ratio leaves",LeafCount/@ratios}];
out[{"triangle",FullSimplify[ratios[[{2,3}]]-ratios[[1]]/3,nu>0]}];
out[{"completeness",FullSimplify[{ratios[[5]]-ratios[[4]]-ratios[[1]]/3,ratios[[6]]+2ratios[[4]]-2ratios[[1]]/3},nu>0,TimeConstraint->10]}];
proj[e_]:=Module[{a,aa1,aa2,aa3,bb,dd},a[u_,v_]:=((e/.{t1->u,t3->v,ph->0})+(e/.{t1->u,t3->v,ph->Pi})+2(e/.{t1->u,t3->v,ph->Pi/2}))/4;aa1=a[Pi/2,Pi/2];aa2=a[Pi/2,0]-aa1;aa3=a[0,0]-aa1-2aa2;bb=((e/.{t1->Pi/4,t3->Pi/4,ph->0})-(e/.{t1->Pi/4,t3->Pi/4,ph->Pi}))/2;dd=((e/.{t1->Pi/2,t3->Pi/2,ph->0})-(e/.{t1->Pi/2,t3->Pi/2,ph->Pi/2}))/2;Factor/@{aa1,aa2,aa3,bb,dd}];
mat=proj/@ratios[[4;;6]];
out[{"matrix",mat}];
Put[{ratios,mat},tmp<>"gauge_soft_reduced.wl"];
out[{"old compact numerical",Table[With[{v=vv},N[({SLLSimp[v,2/5,7/10,9/10],SHHSimp[v,2/5,7/10,9/10],SLHSimp[v,2/5,7/10,9/10]}-f0[-I v,-I v](ratios[[4;;6]]/.{nu->v,t1->2/5,t3->7/10,ph->9/10})),40]],{vv,{1/10,3/5,1}}]}];

tmp="\\\\wsl.localhost\\Ubuntu\\tmp\\";
{cc,mm}=Get[tmp<>"gauge_p_factored.wl"];
cc=Map[Factor[#/.{Log[4]->2Log[2],Log[16]->4Log[2]}]&,cc,{2}];
cc[[All,1]]=Map[Collect[Expand[#],{ps,qh1},Factor]&,cc[[All,1]]];
Put[{cc,mm},tmp<>"gauge_p_display.wl"];
Export[tmp<>"gauge_p_display_report.txt",ToString[{LeafCount/@cc,cc},InputForm],"Text"];
Export[tmp<>"gauge_reduced_sources.json",<|"pMatrix"->ToString[cc,InputForm],"pMasters"->ToString[mm,InputForm],"soft"->ToString[Get[tmp<>"gauge_soft_reduced.wl"],InputForm]|>,"RawJSON"];
Quit[];
