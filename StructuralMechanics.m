(* ::Package:: *)

BeginPackage["StructuralMechanics`"]
(*V0 .4*)

(*Foundamental functions*)
energyParser::usage = "energyParser[EXPR, FUN, VAR, TAR] returns EXPR whose FUN[VAR,t] are replaced with TAR."
varD1D::usage = "varD1D[EXPR, U, Z] returns variation of EXPR with respect to U[Z,t]."
(*Basic utility*)
hamiltonP::usage = "hamiltonP[EXPR, FUN, VAR] returns EOM and BC's from energy EXPR, which is in the form of {T, U, W} with unknowns as FUN[VAR,t]"
lagrangeP::usage = "lagrangeP[EXPR, FUN, SFUN, VAR, RNG] returns EOM from energy EXPR, which is in the form of {T, U, W} with unknowns as FUN[VAR,t]
that will be replaced by shape functions. Set of shape functions is a 2D list. Generalized coordinates will be generated automatically."
(*Extended functions*)
eqnRdc::usage = "eqnRdc[EXPR, IDX] returns reduced EOM by retaining BC's specified by IDX, which is 2-D list"
freqRdc::usage = "freqRdc[EXPR, FUN, VAR] returns reduced EOM by removing temporal dependence, whose unknowns are FUN[VAR]"
evpSolve::usage = "evpSolve[EXPR, FUN, VAR] returns the solution to eigenvalue problem EXPR, whose unknowns are FUN[VAR]"
rrSolve::usage = "rrSolve[EQN] returns Rayleigh-Ritz solution to EQN obtained from lagrangeP"
(*Auxiliary*)
coeSmp::usage = "coeSmp[EXPR, U, LIST] returns simplified EXPR with coefficients of U replaced by those from LIST."
coeAppd::usage = "coeAppd[EXPR, LISTN, APPD] returns simplified EXPR according to LISTN with modification introduced by APPD."

Begin["`Private`"]

energyParser[EXPR_,FUN_,VAR_,TAR_]:=Module[{expr=Expand[EXPR],fun=FUN,var=VAR,tar=TAR,const=4,lst,cvt},
(*Convert shorthand expression to full form*)
cvt=Flatten@Table[ToExpression@StringJoin[ToString[fun[[k]]],Table[ToString[var],{i}],Table["t",{j}]]->D[tar[[k]],{var,i},{ToExpression["t"],j}],{i,0,const},{j,0,const},{k,1,Length[fun]}];
expr/.cvt];

varD1D[EXPR_,FUN_,VAR_]:=Module[{expr=EXPR,fun=FUN,var=VAR,t=ToExpression["t"],res,df,tmp,tmpVar,prt},
(*Variantional operation on Integrals based on Gateaux derivative*)
res=Table[
(*Construct variational symbol*)
df=ToExpression["\[Delta]"<>ToString[f]];
(*Gateaux derivative*)
tmp=D[expr/.{f[a__]->(f[a]+g\[Epsilon] df[a]),Derivative[a__][f][b__]->Derivative[a][f][b]+g\[Epsilon] Derivative[a][df][b]},g\[Epsilon]]/.{g\[Epsilon]->0};
If[var===t,
(*If integral variable is time*)
prt={
Coefficient[tmp/.{Derivative[0,1][df][_,t]->tmpVar},tmpVar],
Coefficient[tmp/.{Derivative[1,1][df][_,t]->tmpVar},tmpVar],
Coefficient[tmp/.{Derivative[2,1][df][_,t]->tmpVar},tmpVar]};
{df,-D[prt[[1]],t],-D[prt[[2]],t],-D[prt[[3]],t]},
(*If integral variable is space*)
prt={
Coefficient[tmp/.{df[var,t]->tmpVar},tmpVar],
Coefficient[tmp/.{Derivative[1,_][df][var,_]->tmpVar},tmpVar],
Coefficient[tmp/.{Derivative[2,_][df][var,_]->tmpVar},tmpVar]};
{df,prt[[1]]-D[prt[[2]],var]+D[prt[[3]],var,var],prt[[2]]-D[prt[[3]],var],prt[[3]]}]
,{f,fun}];
(*Result*)
res];

hamiltonP[EXPR_,FUN_,VAR_]:=
Module[{expr=EXPR,fun=FUN,var=VAR,tuw,\[CapitalLambda]=ToExpression[{"\[Lambda]l","\[Lambda]r"}],t=ToExpression["t"],eb,eqT,eqU,eqTb,eqUb,bcl,bcr,res},
(*Derive Hamiltonian EOM of a beam with BC's on both ends*)
(*Convert energy expression to normal form*)
tuw=energyParser[expr,fun,var,Table[f[var,t],{f,fun}]];
(*Seperate domain terms and boundary terms*)
eb=CoefficientArrays[tuw,\[CapitalLambda]]//Normal;
eb=If[Length[eb]==1,{eb[[1]],Table[{0,0},{Dimensions[eb][[2]]}]},eb];
(*Do variation on kinetic energy and the rest, domain terms*)
eqT=varD1D[eb[[1,1]],fun,t];
eqU=varD1D[eb[[1,2]]-eb[[1,3]],fun,var]; (*External work is included*)
(*Do variation on kinetic energy and the rest, boundary terms*)
eqTb=varD1D[eb[[2,1]],fun,t];
eqUb=varD1D[eb[[2,2]]-eb[[2,3]],fun,var]; (*External work is included*)
(*Assemble the EOM*)
res=Table[
bcl={
{fun[[idx]][var,t],-eqU[[idx,3]]-eqT[[idx,3]]+eqTb[[idx,2,1]]+eqUb[[idx,2,1]]},
{D[fun[[idx]][var,t],var],-eqU[[idx,4]]-eqT[[idx,4]]+eqTb[[idx,3,1]]+eqUb[[idx,3,1]]}};
bcr={
{fun[[idx]][var,t],eqU[[idx,3]]-eqT[[idx,3]]+eqTb[[idx,2,2]]+eqUb[[idx,2,2]]},
{D[fun[[idx]][var,t],var],eqU[[idx,4]]-eqT[[idx,4]]+eqTb[[idx,3,2]]+eqUb[[idx,3,2]]}};
{eqT[[idx,1]],eqU[[idx,2]]-eqT[[idx,2]],bcl,bcr}
,{idx,1,Length[fun]}];
res
];

lagrangeP[EXPR_,FUN_,SFUN_,VAR_,RNG_]:=
Module[{expr=EXPR,fun=FUN,sfun=SFUN,var=VAR,rng=RNG,\[CapitalLambda]=ToExpression[{"\[Lambda]l","\[Lambda]r"}],t=ToExpression["t"],len=Length@Flatten@SFUN,
qfun,cnt,tar,lentmp,tuw,tuws,ener,dq,KK,MM,res,f0,fq,fdq},
(*Derive Lagrangian EOM of a beam with BC's on both ends*)
(*Convert energy expression to discretized form*)
qfun=Table[ToExpression["q"<>ToString[idx]<>"[t]"],{idx,1,len}];
cnt=0;
tar=Table[lentmp=Length[sfun[[idx]]];cnt=cnt+lentmp;qfun[[cnt-lentmp+1;;cnt]].sfun[[idx]],{idx,1,Length@sfun}];
tuw=energyParser[expr,fun,var,tar];
(*Seperate domain terms and boundary terms*)
tuws=CoefficientArrays[tuw[[1]]-tuw[[2]]+tuw[[3]],\[CapitalLambda]]//Normal;
tuws=If[Length[tuws]==1,{tuws[[1]],{0,0}},tuws];
(*Derive the total energy, including boundary terms*)
ener=Expand[Integrate[tuws[[1]],Flatten@{var,rng}]+(tuws[[2,1]]/.{var->rng[[1]]})+(tuws[[2,2]]/.{var->rng[[2]]})];
(*Derive the K and M matrix*)
dq=D[qfun,t];
KK=Table[D[ener,{qfun[[i]],1},{qfun[[j]],1}],{i,1,len},{j,1,len}];
MM=Table[D[ener,{dq[[i]],1},{dq[[j]],1}],{i,1,len},{j,1,len}];
(*Derive the generalized forces*)
res=Expand[ener-dq.MM.dq/2-qfun.KK.qfun/2];
fq=Table[D[res,{qfun[[i]],1}],{i,1,len}];
fdq=Table[D[res,{dq[[i]],1}],{i,1,len}];
f0=Simplify[res-fq.qfun-fdq.dq];
(*Generate the output*)
{MM,-KK,f0,fq,fdq}
];

eqnRdc[EXPR_,IDX_]:=Module[{expr=EXPR,idx=IDX,tmp},
(*Simplify EOM by specifying the BC's*)
Table[
tmp=Switch[Length[idx[[jdx]]],
(*For 2nd-order equations*)
2,
{{expr[[jdx,3,1,idx[[jdx,1]]]]},
{expr[[jdx,4,1,idx[[jdx,2]]]]}},
(*For 4th-order equations*)
4,
{{expr[[jdx,3,1,idx[[jdx,1]]]],expr[[jdx,3,2,idx[[jdx,2]]]]},
{expr[[jdx,4,1,idx[[jdx,3]]]],expr[[jdx,4,2,idx[[jdx,4]]]]}}];
(*Variable, EOM, Left BC, Right BC*)
{expr[[jdx,1]],expr[[jdx,2]],tmp[[1]],tmp[[2]]},{jdx,1,Length[expr]}]];

freqRdc[EXPR_,FUN_,VAR_]:=Module[{expr=EXPR,fun=FUN,var=VAR,p=ToExpression["p"],t=ToExpression["t"],cvt,fsub,coe,res},
(*Seperate temporal and spatial dependence and get the L & M operators*)
(*Construct conversion from [var,t] to [var]*)
cvt=Flatten@Table[{Derivative[i_,j_][f][var,t]->D[f[var],{var,i}](I p)^j,f[var,t]->f[var]},{f,fun}];
(*Function for extracting the operators*)
fsub=Switch[Length[#],2,{-#[[1]],#[[2]]},1,{#[[1]],0},0,{0,0}]&;
res=Table[
{expr[[idx,1]],
fsub[CoefficientList[expr[[idx,2]]/.cvt,p^2]],
Table[fsub[CoefficientList[expr[[idx,3,jdx]]/.cvt,p^2]],{jdx,1,Length[expr[[idx,3]]]}],
Table[fsub[CoefficientList[expr[[idx,4,jdx]]/.cvt,p^2]],{jdx,1,Length[expr[[idx,4]]]}]},
{idx,1,Length[expr]}]];

evpSolve[EXPR_,FUN_,VAR_]:=Module[{expr=EXPR,fun=FUN,var=VAR,p=ToExpression["p"],l=ToExpression["l"],len=Length[EXPR],eqn,bcl,bcr,sol,coe,cvt0,cvt,det},
(*Solve eigenvalue problem defined by L & M operators*)
(*Construct the EOM from L & M*)
eqn=Table[expr[[idx,2,1]]-p^2 expr[[idx,2,2]]==0,{idx,1,len}];
(*Construct the BC's, left and right respectively*)
bcl=Flatten@Table[expr[[idx,3,All,1]]-p^2 expr[[idx,3,All,2]],{idx,1,len}];
bcr=Flatten@Table[expr[[idx,4,All,1]]-p^2 expr[[idx,4,All,2]],{idx,1,len}];
(*Solve the EOM without BC's*)
sol=DSolve[eqn,Table[f[var],{f,fun}],var];
(*Generate list of constant coefficients*)
coe=Table[ToExpression["C["<>ToString[idx]<>"]"],{idx,1,Length[Sum[sol[[1,sdx,2]],{sdx,1,len}]]}];
(*Generate substitution of derivatives*)
cvt0=Table[sol[[1,idx,1]]->Coefficient[sol[[1,idx,2]],coe],{idx,1,len}];
cvt=Flatten[Table[D[cvt0,{var,idx}],{idx,0,4}],1];
(*Generate the determinant*)
det=Flatten[{bcl/.cvt/.{var->0},bcr/.cvt/.{var->l}},1];
(*Print middle results and solve for eigenvalue*)
{sol,det,Solve[Det[det]==0,p]}
];

rrSolve[EQN_]:=Module[{eqn=EQN,res},
(*Apply Rayleigh-Ritz/Galerkin method to the beam problem described by L & M operators, based on lagrangeP*)
(*Solve EVP*)
res=Eigensystem[Inverse[eqn[[1]]].eqn[[2]]];
(*Construct output*)
{Table[Sqrt[res[[1,idx]]],{idx,1,Length[res[[1]]]}],res[[2]]}
];

coeSmp[EXPR_,U_,LIST_]:=Module[{expr=Expand[EXPR],u=U,list=LIST,n=Length[U],clr,lst,sub,tmpVar},
clr=Table[u[[i]]->0,{i,1,n}];
sub=Table[{u[[i]],list[[i]],Coefficient[expr/.{u[[i]]->tmpVar},tmpVar]/.clr},{i,1,n}];
sub=Append[sub,{1,ToExpression[ToString[list[[1]]]<>"0"],Simplify[expr-Sum[sub[[i,1]]sub[[i,3]],{i,1,n}]]}];
{Sum[sub[[i,1]]sub[[i,2]]Boole[!PossibleZeroQ[sub[[i,3]]]],{i,1,n+1}],sub,Table[sub[[i,2]]->sub[[i,3]],{i,1,n+1}]}];

coeAppd[EXPR_,LISTN_,APPD_]:=EXPR/.Table[lst->ToExpression[ToString[lst]<>APPD],{lst,LISTN}];

End[]
EndPackage[]



