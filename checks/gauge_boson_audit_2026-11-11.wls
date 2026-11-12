(* Read-only audit of the 2026-11-11 gauge notebook and manuscript formulas.
Run in a fresh Wolfram 15 kernel:
  WolframKernel -script checks/gauge_boson_audit_2026-11-11.wls
This reads the notebook and P2 script; it does not save or evaluate the archive's Get cell.
Inputs are selected by index against the recorded notebook SHA256 to avoid silent drift.
All notebook expressions are held before evaluation; BeginPackage cells are skipped so
this standalone audit consistently evaluates definitions in Global`.
*)
$HistoryLength=0;
root=DirectoryName[DirectoryName[$InputFileName]];
path=FileNameJoin[{root,"Gauge Boson Signal.nb"}];
If[FileHash[path,"SHA256","HexString"]!="c232ddaf2a2ef7749f555458962db10832bfe04cdb7df2aa75270b0b313aa9e6",Print["Notebook changed; update cell selection before running this audit."];Quit[2]];
nb=Get[path];
inputs=Flatten[List@@Cases[nb,Cell[BoxData[b_],"Input",___]:>ToExpression[b,StandardForm,HoldComplete],Infinity]];
Scan[ReleaseHold,inputs[[Join[Range[2,15],{19,20},Range[34,49],{207}]]]];
out[s_,v_]:=Print[s," = ",ToString[v,InputForm,PageWidth->Infinity]];
out["Notebook identities",{loopIdentities,loopReconstruction,a3Identity}];
(* Independently transcribed from paper equations 115-117. *)
paperT=gg/hh^2 ((1-w^-2)(ps^2-qh1-(4u ps+2)/(1+2u))-4I Pi/(1+w^2)(ps-2u/(1+2u))+Pi^2(1+4(w^2-1)/(1+w^2)^2))-I((2+u)w-u/w)/(1+2u)+I(-w master[1]/(4u^2)+2u master[2]/w+master[4]/w);
paperD=4gg/hh^2/((1-2u)(3+2u))((1-w^-2)((1+4u)/(1+2u)+2u ps)-4I Pi u/(1+w^2))+2I/((1-2u)(3+2u))((8+10u+3u^2)w-(6+8u+3u^2)/w)-2I u(1+2u)/((1-2u)w)master[3];
out["Paper finite P coefficient residual",Factor/@Together[gg {(1+2u)/(1-2u)paperT+paperD,paperT}-pTemplates]];
paperB[v_]:=3(v-I)(v-2I)^4 Sinh[Pi v]^2/(256Pi^(7/2)v^2(2v-I))Gamma[-2-I v]^2 Gamma[5/2+I v]^2 Gamma[-3/2-2I v]Gamma[1+2I v];
out["Paper bubble numerical residual",Table[{v,N[paperB[v]/(commonLoop[v]p1Gauge[v]^2(rBubble/.nu->v))-1,35]},{v,{1/100,1/10,3/5,1,2,5}}]];
out["Paper bubble exact ratio",FullSimplify[paperB[nu]/(commonLoop[nu]p1Gauge[nu]^2 rBubble),nu>0,TimeConstraint->30]];
out["Raw vs reduced soft loops",Table[With[{vv=N[v,55]}, {v, N[({bubblePolarizationSum[vv],IL[vv,2/5],IP[vv,2/5],SLL[vv,2/5,7/10,9/10],SHH[vv,2/5,7/10,9/10],SLH[vv,2/5,7/10,9/10]}/(gSoft[vv](softRatios/.{nu->vv,t1->2/5,t3->7/10,ph->9/10})))-1,35]}],{v,{1/100,3/5,2}}]];
out["P samples",Table[{v,hardValues[v,30]},{v,{1/100,1/10,3/5,1,2}}]];
out["Regulator convergence",Table[With[{vv=N[v,60],target=hardValues[v,35]}, {v,Table[{e,N[({P2L[vv,N[e,60]],P2PM[vv,N[e,60]]}-target),30]},{e,{1/1000,1/10000,1/100000}}]}],{v,{1/10,3/5,2}}]];

(* Reconstruct all soft sums directly from S_ij(q)=delta_ij+(beta-1)q_i q_j/q^2.
   qx,qy azimuthal moments and z powers are integrated independently of notebook J definitions. *)
sr[j_,k_]:=sr[j,k]=Factor[FunctionExpand[Pochhammer[2x-3/2,j+k]Pochhammer[3/2-x,-j]Pochhammer[3/2-x,-k]/(Pochhammer[3-2x,-j-k]Pochhammer[x,j]Pochhammer[x,k])]];
jr[r_Integer,m_Integer,n_Integer]:=jr[r,m,n]=Factor[Sum[(-1)^(r-j) Multinomial[j,k,r-j-k]sr[-m/2+r/2-k,n/2-j],{j,0,r},{k,0,r-j}]/2^r];
moment[p_,a_,b_]:=Module[{av,terms},
 av=Expand[Total[Map[Function[row,With[{ij=row[[1]],cc=row[[2]]},If[OddQ[ij[[1]]]||OddQ[ij[[2]]],0,cc qm^Total[ij] z^ij[[3]] (1-z^2)^((ij[[1]]+ij[[2]])/2) Factorial2[ij[[1]]-1]Factorial2[ij[[2]]-1]/Factorial2[ij[[1]]+ij[[2]]]]]],CoefficientRules[Expand[p],{qx,qy,qz}]]]];
 terms=CoefficientRules[av,{qm,z}];Factor[Total[Map[#[[2]]jr[#[[1,2]],#[[1,1]]-2a,2b]&,terms]]]];
softV={qx,qy,qz};ww={qx,qy,qz+1};hardA={Sin[t1],0,Cos[t1]};hardB={Sin[t3]Cos[ph],Sin[t3]Sin[ph],Cos[t3]};
c=4x/(1-2x);
ib=Factor[3+2c+c^2 moment[(softV.ww)^2,1,1]];
it=Factor[1+c(moment[(hardA.softV)^2,1,0]+moment[(hardA.ww)^2,0,1])+c^2 moment[(hardA.softV)(softV.ww)(ww.hardA),1,1]];
ill=Factor[(hardA.hardB)^2+c(moment[(hardA.hardB)(hardA.softV)(hardB.softV),1,0]+moment[(hardA.hardB)(hardA.ww)(hardB.ww),0,1])+c^2 moment[(hardA.ww)(ww.hardB)(hardB.softV)(softV.hardA),1,1]];
out["Independent bubble residual",Factor[(ib/.x->-I nu)-rBubble]];
out["Independent triangle isotropy residual",FullSimplify[it-ib/3,Element[t1,Reals],TimeConstraint->20]];
out["Independent LL box residual",FullSimplify[(ill/.x->-I nu)-softRatios[[4]],nu>0&&Element[{t1,t3,ph},Reals],TimeConstraint->40]];
out["Independent all-group numerical residuals",Table[With[{vv=N[val,50],ar=Thread[{t1,t3,ph}->angles]},N[((({ib,it,(ib-it)/2,ill,ib-2it+ill,2(it-ill)}/.x->-I nu)/softRatios-1)/.Join[{nu->vv},ar]),30]],{val,{1/100,3/5,2}},{angles,{{2/5,7/10,9/10},{Pi/3,Pi/4,Pi/6},{0,Pi/2,0}}}]];
(* General J recurrence versus notebook definitions, at all shifts actually used. *)
out["Independent J recurrence numeric",Table[With[{m=shift[[1]],n=shift[[2]]},N[(jr[r,m,n]-({J0,J1,J2,J3,J4}[[r+1]][m,n,x]/f0[x,x]))/.x->N[-3I/5,50],30]],{r,0,4},{shift,{{0,0},{0,2},{1,2},{2,2}}}]];

out["Box paper five polynomial residuals",gaugeAngularCheck];out["Box paper Gamma normalization ratio",gaugeNCheck];
timeTerm[a_,n_]:=Gamma[2I n]/(I a)^(2I n)+I a Gamma[1+2I n]/(I a)^(1+2I n)-Gamma[2+2I n]/(4(I a)^(2+2I n));
out["Independent P1 ratio",FullSimplify[Sum[I a timeTerm[a,nu],{a,{-1,1}}]/p1Gauge[nu],nu>0,TimeConstraint->20]];
src[t_]=t^(-3/2+I nu)(1+I a t)Exp[-I a t];
out["Longitudinal integration-by-parts residual for a=+-1",Table[FullSimplify[(D[src[t],t]+2src[t]/t-t^(-5/2+I nu)Exp[-I a t]((1/2+I nu)(1+I a t)+t^2))/.a->skSign,t>0],{skSign,{-1,1}}]];
out["Fresh C samples",Table[{v,N[coefficientValues[v,30],16]},{v,{1/10,3/5,1,2}}]];
out["Fresh P extrapolation checks",Table[With[{vv=N[v,60],pv=hardValues[v,35]},Module[{es={1/1000,1/2000,1/4000,1/8000,1/16000},vals,ws},vals=Table[N[{P2L[vv,N[e,60]],P2PM[vv,N[e,60]]},40],{e,es}];ws=Table[Product[If[k==j,1,-es[[k]]/(es[[j]]-es[[k]])],{k,5}],{j,5}];{v,Max[Abs[(ws.vals-pv)/pv]]}]],{v,{1/10,3/5,2}}]];
out["Machine precision regulator loss",Table[{e,N[{P2L[0.6,e],P2PM[0.6,e]}]},{e,{0.001,0.00001,0.00000001}}]];


out["Loading finite-part script",AbsoluteTime[]];
Get[FileNameJoin[{root,"P2_reg0_analytic.wl"}]];
out["Script finite P at nu=3/5 (L,T)",TimeConstrained[N[{P2LReg0[3/5],P2PMReg0[3/5]},25],120,"timeout"]];
out["Pole coefficients at nu=3/5",N[{P2PMPoleCoefficients,P2LPoleCoefficients}/.nu0->3/5,25]];
Quit[];
