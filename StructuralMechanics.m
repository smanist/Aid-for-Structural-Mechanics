(* ::Package:: *)

BeginPackage["StructuralMechanics`"]

energyParser::usage = "energyParser[EXPR, FUN, VAR, TAR] returns EXPR whose FUN[VAR,t] are replaced with TAR."
varD1D::usage = "varD1D[EXPR, U, Z] returns variation of EXPR with respect to U[Z,t]."
hamiltonP::usage = "hamiltonP[TUW, FUN, VAR] returns EOM and BC's from energy TUW, whose unknowns are FUN[VAR,t]"
lagrangeP::usage = "lagrangeP[TUW, FUN, VAR, RNG, QFUN] returns EOM from energy TUW, which is integrated over RNG with unknowns as FUN[VAR,t] and generalized coordinates as QFUN."
coeSmp::usage = "coeSmp[EXPR, U, LIST] returns simplified EXPR with coefficients of U replaced by those from LIST."
coeAppd::uasge = "coeAppd[EXPR, LISTN, APPD] returns simplified EXPR according to LISTN with modification introduced by APPD."

Begin["`Private`"]

energyParser[EXPR_,FUN_,VAR_,TAR_]:=Module[{expr=Expand[EXPR],fun=FUN,var=VAR,tar=TAR,const=4,lst,cvt},
cvt=Flatten@Table[ToExpression@StringJoin[ToString[fun[[k]]],Table[ToString[var],{i}],Table["t",{j}]]->D[tar[[k]],{var,i},{ToExpression["t"],j}],{i,0,const},{j,0,const},{k,1,Length[fun]}];
expr/.cvt];

varD1D[EXPR_,U_,Z_]:=Module[{expr=EXPR,u=U,z=Z,res,df,tmp,tmpVar,prt},
res=Table[
df=ToExpression["\[Delta]"<>ToString[f]];
tmp=D[expr/.{f[a__]->(f[a]+g\[Epsilon] df[a]),Derivative[a__][f][b__]->Derivative[a][f][b]+g\[Epsilon] Derivative[a][df][b]},g\[Epsilon]]/.{g\[Epsilon]->0};
prt={
Coefficient[tmp/.{df[z,ToExpression["t"]]->tmpVar,Derivative[0,1][df][_,ToExpression["t"]]->tmpVar},tmpVar],
Coefficient[tmp/.{Derivative[1,_][df][z,_]->tmpVar,Derivative[1,1][df][_,ToExpression["t"]]->tmpVar},tmpVar],
Coefficient[tmp/.{Derivative[2,_][df][z,_]->tmpVar,Derivative[2,1][df][_,ToExpression["t"]]->tmpVar},tmpVar]};
If[z===ToExpression["t"],
Flatten@{df,-D[prt,ToExpression["t"]]},
{df,prt[[1]]-D[prt[[2]],z]+D[prt[[3]],z,z],prt[[2]]-D[prt[[3]],z],prt[[3]]}]
,{f,u}];
{Sum[res[[i,1]]res[[i,2]],{i,1,Length[u]}],res}];

hamiltonP[TUW_,FUN_,VAR_]:=Module[{e=TUW,u=FUN,z=VAR,eb,eqT,eqTb,eqUW,eqUWb,res},
eb=CoefficientArrays[e,ToExpression["\[Lambda]"]]//Normal;
eb=If[Length[eb]==1,{eb[[1]],Table[{0},{Dimensions[eb][[2]]}]},eb];
eqT=varD1D[eb[[1,1]],u,ToExpression["t"]][[2]];
eqTb=varD1D[eb[[2,1,1]],u,ToExpression["t"]][[2]];
eqUW=varD1D[eb[[1,2]]-eb[[1,3]],u,z][[2]];
eqUWb=varD1D[eb[[2,2,1]]-eb[[2,3,1]],u,z][[2]];
res=eqUW;
res[[All,2;;]]=res[[All,2;;]]-eqT[[All,2;;]];
res[[All,3;;]]=res[[All,3;;]]+eqTb[[All,2;;3]]+eqUWb[[All,2;;3]];
res
];

lagrangeP[TUW_,FUN_,VAR_,RNG_,QFUN_]:=Module[{tuw=TUW,fun=FUN,var=VAR,rng=RNG,qfun=QFUN,len=Length[QFUN],expr,ener,dq,KK,MM,res},
expr=CoefficientList[tuw[[1]]-tuw[[2]]+tuw[[3]],ToExpression["\[Lambda]"]];
expr=If[Length[expr]==1,{expr[[1]],0},expr];
ener=Expand@(Integrate[expr[[1]],Flatten@{var,rng}]+expr[[2]]/.{var->rng[[2]]});
dq=D[qfun,ToExpression["t"]];
KK=Table[D[ener,{qfun[[i]],1},{qfun[[j]],1}],{i,1,len},{j,1,len}];
MM=Table[D[ener,{dq[[i]],1},{dq[[j]],1}],{i,1,len},{j,1,len}];
res=Expand[ener-dq.MM.dq/2-qfun.KK.qfun/2];
{MM,KK,Table[D[res,{qfun[[i]],1}],{i,1,len}],Table[D[res,{dq[[i]],1}],{i,1,len}]}
];

coeSmp[EXPR_,U_,LIST_]:=Module[{expr=Expand[EXPR],u=U,list=LIST,n=Length[U],clr,lst,sub,tmpVar},
clr=Table[u[[i]]->0,{i,1,n}];
sub=Table[{u[[i]],list[[i]],Coefficient[expr/.{u[[i]]->tmpVar},tmpVar]/.clr},{i,1,n}];
sub=Append[sub,{1,ToExpression[ToString[list[[1]]]<>"0"],Simplify[expr-Sum[sub[[i,1]]sub[[i,3]],{i,1,n}]]}];
{Sum[sub[[i,1]]sub[[i,2]]Boole[!PossibleZeroQ[sub[[i,3]]]],{i,1,n+1}],sub,Table[sub[[i,2]]->sub[[i,3]],{i,1,n+1}]}];

coeAppd[EXPR_,LISTN_,APPD_]:=EXPR/.Table[lst->ToExpression[ToString[lst]<>APPD],{lst,LISTN}];

End[]
EndPackage[]



