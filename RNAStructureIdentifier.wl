(* ::Package:: *)

BeginPackage["Bisect`"];


Structure2Pair::usage="Input a structure s in dot-bracket notation.\[IndentingNewLine] Output a pair table.";


DrawESB::usage="DrawESB[esb_]:=\[IndentingNewLine]Input an ensemble of secondary structures using pair tables. 
Draw a secondary structure given by base pairs and theirs probabilities p.";


TreeMEBP::usage="TreeMEBP[esb_,level_:10]: 
Input an ensemble of secondary structures using pair tables, \[IndentingNewLine]Return the set of all base pairs (i,j) for each structure in the ensemble, \[IndentingNewLine]  the vector of base pairs for each structure in the ensemble, \[IndentingNewLine]the MEBPSpace of all maximum entropy base pairs of all clusters,\[IndentingNewLine]the TreeSpace containing all clusters, each level of clusters representing the ones containing the maximum entropy base pairs or not.";


RanPdtMEBP::usage="RanPdtMEBP[TreeMEBP_,MFE_,esb_,lvl_:10,error_:0.01,error2_:0.05]: 
Input:\[IndentingNewLine] TreeMEBP, which is the output of TreeMEBP.\[IndentingNewLine]   MFE, target structure,\[IndentingNewLine]   esb, the ensemble of structures.\[IndentingNewLine]error: the error rate when the truthful answer is Yes.\[IndentingNewLine]error2: the error rate when the truthful answer is No.\[IndentingNewLine]\[IndentingNewLine]Return {str,ratio, leaveinfo, pathinfo}\[IndentingNewLine]
str denotes the structure predicted in the leave.\[IndentingNewLine]ratio denotes the ratio of the predicted structure in the leave.

leaveinfo={leave,leavecorrectness}
leave: the leave ensemble at the end of the path.\[IndentingNewLine]leavecorrectness: whether choose a correct leaf.\[IndentingNewLine]1 represents correct, 0 represents wrong.
\[IndentingNewLine]pathinfo={path, pathcorrectness,pathent,pathbpdist,bpset},\[IndentingNewLine]\[IndentingNewLine]path(result): the path to the cluster closest to MFE.\[IndentingNewLine]The path is a list of length 10.\[IndentingNewLine]For i-th position, 1 represents that the predicted structure contains the i-th MEBP.\[IndentingNewLine]2 represents that predicted does not contain the i-th MEBP.\[IndentingNewLine] 0 represnts that clustering stops.\[IndentingNewLine]\[IndentingNewLine]pathcorrectness: the path correctness to the cluster closest to MFE.\[IndentingNewLine]For i-th position, 1 represents correct, 0 represents wrong, -1 represents initial value.\[IndentingNewLine]\[IndentingNewLine]pathent: the entropy of the cluster on the path. -1 represents initial value.\[IndentingNewLine]\[IndentingNewLine]pathbpdist: bpdist between the cluster and MFE. -1 represents initial value.\[IndentingNewLine]\[IndentingNewLine]bpset: the set of bias base pairs selected on the path.";


DrawTree::usage="DrawTree[esb_,TreeMEBP_,lvl_]:=\[IndentingNewLine]Input an ensemble of secondary structures, and the output of TreeMEBP, and the max level of the tree (<=3).\[IndentingNewLine]Return the drawing of the tree of clusters, each level of clusters representing the ones containing the maximum entropy base pairs or not.";


Begin["`Private`"]


Structure2Pair[s_]:=
(* input a structure in dot-bracket notation
 output a pair table *)
Module[{n=StringLength[s],i,upstream={},ptable},
ptable=Table[0,{n}];
For[i=1,i<=n,i++,
If[StringTake[s,{i}]=="(",
AppendTo[upstream,i],
If[StringTake[s,{i}]==")",
ptable[[i]]=upstream[[-1]];
ptable[[upstream[[-1]]]]=i;
upstream=Drop[upstream,-1],
0
]
]
];
Return[ptable]
]


DrawESB[esb_]:=
(* Input an ensemble of secondary structures using pair tables. Draw a secondary structure given by base pairs and theirs probabilities p.  *)
Module[{tbp=TallyBPESB[esb],m=Length[esb],edge={},i,bp,n=Length[esb[[1]]]},
(* m is the number of structures in the ensemble, n is the length of each structure, bp is the number of base pairs.*)
bp=Length[tbp];
For [i=1,i<=bp,i++,
(*AppendTo[edge,UndirectedEdge[i,i+1]];*)
AppendTo[edge,Style[UndirectedEdge[tbp[[i,1,1]],tbp[[i,1,2]]] ,GrayLevel[1-tbp[[i,2]]/m]] ]
];
Return[Graph[Range[n],edge,VertexCoordinates->Table[{i,0},{i,1,n}],VertexLabels->None,(*EdgeShapeFunction\[Rule]Table[(i\[Rule]i+1)\[Rule]"DashedLine",{i,1,n=Length[x]-1}],*)
GraphLayout->"LinearEmbedding"]];
]


BPSet[x_]:=
(* Return all base pairs (i,j) in the secondary structure*)
Module[{ n=Length[x],BpSet={},i,j},
For[i=1,i<=n,i++,
 If[x[[i]]>i, 
AppendTo[BpSet,{i,x[[i]]}]
,0]
];
Return[BpSet]]


BPESB[esb_]:=
Map[BPSet,esb]


TallyBPESB[esb_]:=
Sort[Tally[Flatten[Map[BPSet,esb],1]],#1[[2]]>#2[[2]]&]


BPSetESB[esb_]:=
Module[{tally=TallyBPESB[esb],i},
Return[Table[tally[[i,1]],{i,1,Length[tally]}]]
]


Entf[flist_,n_]:=
(* Return the entropy of a random variable given by its frequency in an ensemble.  *)
Module[{sum=0,i},
For[i=1,i<=Length[flist],i++,
If[flist[[i]]<n&&flist[[i]]>0,sum=sum-flist[[i]]/n*Log[2,flist[[i]]/n],0]
];
Return[N[sum]]
]


StrEnt[esb_]:=
(* Return the structural entropy for the ensemble.*)
Module[{t=Tally[esb],s,i,m=Length[esb],ent=0},
s=Length[t];
For[i=1,i<=s,i++,
ent=ent-t[[i,2]]/m*Log[2,N[t[[i,2]]/m]];
];
Return[ent]
]


MaxEntfBP[esb_]:=
(*
Input an ensemble. Return  the maximum entropy base pair.
*)
Module[{bpEnt=BPEntf[esb],tbp=TallyBPESB[esb],idx},
idx=Position [bpEnt, Max[bpEnt]][[1,1]];
Return[tbp[[idx]]]
]


VectorEnt[vector_]:=
(* Return entropies of a vector ensemble. *)
Module[{m=Length[vector],b=Length[vector[[1]]],i,n,list},
list=Table[Count[vector,n_/;n[[i]]==1],{i,1,b}];
Return[Table[Entf[{list[[i]],m-list[[i]]},m],{i,1,b}]]
]


BPVector[esb_]:=
(* Input an ensemble.
Return the base-pair vectors of structures in the ensemble.*)
Module[{bp=BPSetESB[esb],bpe=BPESB[esb],vc,n,m,i,j},
m=Length[esb];
n=Length[bp];
vc=Table[Table[0,{j,1,n}],{i,1,m}];
For[i=1,i<=m,i++,
For[j=1,j<=Length[bpe[[i]]],j++,
If[MemberQ[bp,bpe[[i,j]]],vc[[i,Position[bp,bpe[[i,j]]][[1,1]]]]=1,0]
];
];
Return[vc]
]


BP[esb_]:=
(* Return the set of all base pairs (i,j) for each structure in the ensemble 
and the vector of base pairs for each structure in the ensemble.*)
Module[{bpset,TallyBp,bpe,vc,n,m,i,j},
bpe=BPESB[esb];

TallyBp=Sort[Tally[Flatten[bpe,1]],#1[[2]]>#2[[2]]&];

bpset=Table[TallyBp[[i,1]],{i,1,Length[TallyBp]}];

m=Length[esb];
n=Length[bpset];
vc=Table[Table[0,{j,1,n}],{i,1,m}];
For[i=1,i<=m,i++,
For[j=1,j<=Length[bpe[[i]]],j++,
If[MemberQ[bpset,bpe[[i,j]]],vc[[i,Position[bpset,bpe[[i,j]]][[1,1]]]]=1,0]
];
];
Return[{bpset,vc}]
]


ClusterMeV[vector_]:=
(*
Input the vector form for an ensemble. 
Return {idx,{esb1,esb2}}
idx is the index of the MEBP.
esb1 is the clusters containing the MEBP.
esb2 is the cluseter not conttaining the MEBP.*)
Module[{bpent,esb1={},esb2={},idx,j,m=Length[vector]},
bpent=VectorEnt[vector];
\:3000If[Max[bpent]==0.0,
esb1=Range[m],

idx=Position [bpent, Max[bpent]][[1,1]];
For[j=1,j<=m,j++,
If[vector[[j,idx]]==1,AppendTo[esb1,j],AppendTo[esb2,j]]
]
];
If[Length[esb2]==0||Length[esb1]==0,
Return[{0,{esb1,esb2}}],
Return[{idx,{esb1,esb2}}]]
(*If esb2 is empty, it means all structures are the same, MEBP is not valid.*)
]


TreeMEBP[esb_,level_:10]:=
(*
Input an ensemble. 
Return the set of all base pairs (i,j) for each structure in the ensemble, 
  the vector of base pairs for each structure in the ensemble, 
the MEBPSpace of all maximum entropy base pairs of all clusters,
the TreeSpace containing all clusters, each level of clusters representing the ones containing the maximum entropy base pairs or not.
*)
Module[{TreeSpace,MEBPSpace,bp=BP[esb],vector,
Cluster,idx={},idx1,idx2,i,j,m= Length[esb],lvl, result},
vector=bp[[2]];

lvl=level;

For[i=0,i<=lvl,i++,
idx=Union[idx,Tuples[{1,2},i]]];

TreeSpace=Table[{},{i,1,Length[idx]}];

MEBPSpace=Table[0,{i,1,Length[idx]}];


TreeSpace[[1]]=Range[m];




For[j=1,j<=Length[idx],j++,

If[Length[idx[[j]]]<lvl,
idx1=Position[idx,Append[idx[[j]],1]][[1,1]] ;
  
idx2=Position[idx,Append[idx[[j]],2]][[1,1]] ;
 
If[TreeSpace[[j]]=={}
,
MEBPSpace[[j]]=0;
TreeSpace[[idx1]]={};
TreeSpace[[idx2]]={}
,
Cluster=ClusterMeV[vector[[TreeSpace[[j]]]]];


MEBPSpace[[j]]=Cluster[[1]];
TreeSpace[[idx1]]=TreeSpace[[j]][[Cluster[[2,1]]]];

TreeSpace[[idx2]]=TreeSpace[[j]][[Cluster[[2,2]]]]

]
,0]
];
result={bp[[1]],bp[[2]],MEBPSpace,TreeSpace};
Return[result]
]



TreeMaxEntfBP10d[TreeSpace_]:=
(*
Input an TreeSpace_ of secondary structures using pair tables, Return the information about the level 10 clusters, each level of clusters representing the ones containing the maximum entropy base pairs or not.

*)
Module[{LevelEnt,ClusterEnt,ClusterWt,maxbp,idx={},tem,i,j,m= Length[TreeSpace[[1]]],str7=TreeSpace[[1,7]],str7info={},LeavesEntMean,LeavesEntVar,LeavesEty,pos},



ClusterEnt={};

ClusterWt={};


For[i=1,i<=1024,i++,

tem=TreeSpace[[1023+i]];

ClusterEnt=Append[ClusterEnt,StrEnt[tem]];


ClusterWt=Append[ClusterWt,N[ Length[tem]/m]];



If[MemberQ[tem,str7],
str7info={Length[tem],Count[tem,str7],StrEnt[tem]}
,0]

];



LeavesEty=Count[ClusterWt,0.];



pos=Position[ClusterWt,n_/;n!=0.];


ClusterWt=Extract[ClusterWt,pos];

ClusterEnt=Extract[ClusterEnt,pos];


LeavesEntMean=ClusterEnt.ClusterWt;

LeavesEntVar=((ClusterEnt-LeavesEntMean)^2).ClusterWt;


Return[{{LeavesEntMean,LeavesEntVar,LeavesEty,Max[ClusterEnt]},str7info}]
]



TreeMaxEntfBP10e[TreeSpace_,lvl_:10]:=
(*
Input an TreeSpace_ of secondary structures using pair tables, Return the information about the clusters of each level, each level of clusters representing the ones containing the maximum entropy base pairs or not.
*)
Module[{LevelEnt,ClusterEnt,ClusterWt,maxbp,idx={},tem,i,j,k,r,initial, m= Length[TreeSpace[[1]]],str7=TreeSpace[[1,7]],Levelstr7info=Range[10],LevelEntMean=Range[10],LevelEntVar=Range[10],LevelEty=Range[10],LevelMaxEnt=Range[10]},


ClusterEnt={};


ClusterWt={};


For[k=1,k<=lvl,k++,

initial=Sum[2^r,{r,0,k-1}];

For[i=1,i<=2^k,i++,

tem=TreeSpace[[initial+i]];

ClusterEnt=Append[ClusterEnt,StrEnt[tem]];

ClusterWt=Append[ClusterWt,N[ Length[tem]/m]];



If[MemberQ[tem,str7],
Levelstr7info[[k]]={Length[tem],Count[tem,str7],StrEnt[tem]}
,0]

];



LevelEty[[k]]=Count[ClusterWt,0.];


LevelEntMean[[k]]=ClusterEnt.ClusterWt;


LevelEntVar[[k]]=((ClusterEnt-LevelEntMean[[k]])^2).ClusterWt;


LevelMaxEnt[[k]]=Max[ClusterEnt];

ClusterEnt={};
ClusterWt={}


];


Return[{LevelEntMean,LevelEntVar,LevelEty,LevelMaxEnt,Levelstr7info}]
]



BpDistanceE[s1_,esb_]:=
(* The average base pair distance of one structure and a ensemble. *)
Module[{tally=Tally[esb],m=Length[esb],i,dis=0,bpset1=BPSet[s1],bpset2},
bpset2=Table[BPSet[tally[[i,1]]],{i,1,Length[tally]}];
Do[
dis=dis+N[tally[[i,2]]/m]*(Length[Complement[bpset1,bpset2[[i]]]]+Length[Complement[bpset2[[i]],bpset1]])
,{i,1,Length[tally]}];
Return[dis]
]


DomESB[esb_]:=
(* Input an ensemble. 
Return the dominating structure.*)
Module[{t},
t=Sort[Tally[esb],#1[[2]]>#2[[2]]&];
Return[{t[[1,1]],N[t[[1,2]]/Length[esb]]}]
]


ClusterMeVbias[vector_,indicator_,delta_:0.9]:=
(*
Input the vector for an ensemble, 
Return {idx,{esb1,esb2},{biasRank,ratio}}
idx is the index of the MEBP.
esb1 is the clusters containing the MEBP.
esb2 is the cluseter not conttaining the MEBP.
the clusters containing the maximum entropy base pairs or not.

biasRank stores all rankings of MEBP bias (the bias base pair of the cluster)
-2 initial value, clusters without computing bias MEBP
-1 cluster is empty
-3 cluster is non-empty but all structures are the same.
-4 cluster is non-empty but bias MEBP with indicator 0, but entropy ratio to MEBP without bias is < delta=0.9, we choose MEBP has indicator 1
0 all base pairs has indicator 1, meaning MFE contains all remaining base pairs, we choose MEBP has indicator 1.


ratio stores the ratio of the entropy of bias base pair compared with max entropy base pair
 of the cluster

-2 initial value
-1 cluster is empty
-3 cluster is non-empty but all structures are the same.
always non zero
1 means either we pick MEBP having indicator 0
or all base pairs has indicator 1, meaning MFE contains all remaining base pairs, we choose MEBP has indicator 1.
 or bias MEBP with indicator 0, but entropy ratio to MEBP without bias is < delta=0.9 thus we choose MEBP has indicator 1
*)
Module[{bpent,esb1={},esb2={},idx,i,j,m=Length[vector],order,biasRank=0,biasIndex,ratio},
bpent=VectorEnt[vector];
\:3000If[Max[bpent]==0.0,

esb1=Range[m];
biasRank=-3;
ratio=-3
,

order=Ordering[bpent,All,Greater];

Do[
If[indicator[[order[[i]]]]==0,

biasRank=i;

Break[],
0]
,
{i,1,Length[order]}];

If[biasRank==0,

idx=order[[1]];

ratio=1,

biasIndex=order[[biasRank]];

ratio=bpent[[biasIndex]]/Max[bpent];


If[ratio>=delta,


idx=biasIndex
,

biasRank=-4;
ratio=1;

idx=order[[1]]

]
];

For[j=1,j<=m,j++,
If[vector[[j,idx]]==1,AppendTo[esb1,j],AppendTo[esb2,j]]
]
];
If[Length[esb2]==0||Length[esb1]==0,
Return[{0,{esb1,esb2},{biasRank,ratio}}],
Return[{idx,{esb1,esb2},{biasRank,ratio}}]]
(*If esb2 is empty, it means all structures are the same, MEBP is not valid.*)
]


ClusterMeVbiasProbe[vector_,probeindicator_,delta_:0.9]:=
(*
Input the vector for an ensemble of secondary structures using pair tables, 
Return {idx,{esb1,esb2},{biasRank,ratio}}
idx is the index of the MEBP.
esb1 is the clusters containing the MEBP.
esb2 is the cluseter not conttaining the MEBP.
the clusters containing the maximum entropy base pairs or not.

biasRank stores all rankings of MEBP bias (the bias base pair of the cluster)
-2 initial value, clusters without computing bias MEBP
-1 cluster is empty
-3 cluster is non-empty but all structures are the same.
-4 cluster is non-empty but bias MEBP with indicator 0, but entropy ratio to MEBP without bias is < delta=0.9, we choose MEBP has indicator 1
0 all base pairs has indicator 1, meaning MFE contains all remaining base pairs, we choose MEBP has indicator 1.


ratio stores the ratio of the entropy of bias base pair compared with max entropy base pair
 of the cluster

-2 initial value
-1 cluster is empty
-3 cluster is non-empty but all structures are the same.
always non zero
1 means either we pick MEBP having indicator 0
or all base pairs has indicator 1, meaning MFE contains all remaining base pairs, we choose MEBP has indicator 1.
 or bias MEBP with indicator 0, but entropy ratio to MEBP without bias is < delta=0.9 thus we choose MEBP has indicator 1
*)
Module[{bpent,esb1={},esb2={},idx,i,j,m=Length[vector],order,biasRank=0,biasIndex,ratio},
bpent=VectorEnt[vector];
\:3000If[Max[bpent]==0.0,

esb1=Range[m];
biasRank=-3;
ratio=-3

,

order=Ordering[bpent,All,Greater];

Do[
If[probeindicator[[order[[i]]]]==0,

biasRank=i;

Break[],
0]
,
{i,1,Length[order]}];

If[biasRank==0,

idx=order[[1]];

ratio=1
,

biasIndex=order[[biasRank]];

ratio=bpent[[biasIndex]]/Max[bpent];


If[ratio>=delta,

idx=biasIndex
,

biasRank=-4;
ratio=1;

idx=order[[1]]
]
];

For[j=1,j<=m,j++,
If[vector[[j,idx]]==1,AppendTo[esb1,j],AppendTo[esb2,j]]
]
];
If[Length[esb2]==0||Length[esb1]==0,
Return[{0,{esb1,esb2},{biasRank,ratio}}],
Return[{idx,{esb1,esb2},{biasRank,ratio}}]]
(*If esb2 is empty, it means all structures are the same, MEBP is not valid.*)
]


TreeMEBPbias[esb_,MFE_,level_:10,delta_:0.9]:=
(*
Input an ensemble of secondary structures using pair tables, 
Return {bp[[1]],bp[[2]],MEBPSpace,TreeSpace,biasRankSpace,ratioSpace}
Return the set of all base pairs (i,j) for each structure in the ensemble, 
  the vector of base pairs for each structure in the ensemble, 
 MEBPSpace of all bias maximum entropy base pairs of all clusters,

TreeSpace containing all clusters, each level of clusters representing the ones containing the bias maximum entropy base pairs or not.

biasRankSpace stores all rankings of MEBP bias (the bias base pair of the cluster)
-2 initial value, clusters without computing bias MEBP
-1 cluster is empty
-3 cluster is non-empty but all structures are the same.
-4 cluster is non-empty but bias MEBP with indicator 0, but entropy ratio to MEBP without bias is < delta=0.9, we choose MEBP has indicator 1
0 all base pairs has indicator 1, meaning MFE contains all remaining base pairs, we choose MEBP has indicator 1.


ratioSpace stores the ratio of the entropy of bias base pair compared with max entropy base pair
 of the cluster

-2 initial value
-1 cluster is empty
-3 cluster is non-empty but all structures are the same.
always non zero
1 means either we pick MEBP having indicator 0
or all base pairs has indicator 1, meaning MFE contains all remaining base pairs, we choose MEBP has indicator 1.
 or bias MEBP with indicator 0, but entropy ratio to MEBP without bias is < delta=0.9 thus we choose MEBP has indicator 1

*)
Module[{TreeSpace,MEBPSpace,biasRankSpace,ratioSpace,bp=BP[esb],vector,mfe=BPSet[MFE],indicator,
Cluster,idx={},idx1,idx2,i,j,m= Length[esb],lvl, result},
vector=bp[[2]];

indicator=Table[Count[mfe,bp[[1,i]]],{i,1,Length[bp[[1]]]}];

lvl=level;

For[i=0,i<=lvl,i++,
idx=Union[idx,Tuples[{1,2},i]]];


TreeSpace=Table[{},{i,1,Length[idx]}];


MEBPSpace=Table[0,{i,1,Length[idx]}];

biasRankSpace=Table[-2,{i,1,Length[idx]}];


ratioSpace=Table[-2,{i,1,Length[idx]}];


TreeSpace[[1]]=Range[m];




For[j=1,j<=Length[idx],j++,

If[Length[idx[[j]]]<lvl,
idx1=Position[idx,Append[idx[[j]],1]][[1,1]] ;
 
idx2=Position[idx,Append[idx[[j]],2]][[1,1]] ;

If[TreeSpace[[j]]=={}
,
MEBPSpace[[j]]=0;
TreeSpace[[idx1]]={};
TreeSpace[[idx2]]={};

biasRankSpace[[j]]=-1;
ratioSpace[[j]]=-1

,
Cluster=ClusterMeVbias[vector[[TreeSpace[[j]]]],indicator,delta];

MEBPSpace[[j]]=Cluster[[1]];
TreeSpace[[idx1]]=TreeSpace[[j]][[Cluster[[2,1]]]];

TreeSpace[[idx2]]=TreeSpace[[j]][[Cluster[[2,2]]]];

biasRankSpace[[j]]=Cluster[[3,1]];
ratioSpace[[j]]=Cluster[[3,2]]
]
,0]
];
result={bp[[1]],bp[[2]],MEBPSpace,TreeSpace,biasRankSpace,ratioSpace};
Return[result]
]



RanPdtMEBP[TreeMEBP_,MFE_,esb_,lvl_:10,error_:0.01,error2_:0.05]:=
(*
Input:
 TreeMEBP, which is the output of TreeMEBP.
   MFE, target structure,
   esb, the ensemble of structures.
error: the error rate when the truthful answer is "Yes".
error2: the error rate when the truthful answer is "No".

Return {str,ratio, leaveinfo, pathinfo}

str denotes the structure predicted in the leave.
ratio denotes the ratio of the predicted structure in the leave.

leaveinfo={leave,leavecorrectness}
leave: the leave ensemble at the end of the path.
leavecorrectness: whether choose a correct leaf.
1 represents correct, 0 represents wrong.

pathinfo={path, pathcorrectness,pathent,pathbpdist,bpset},

path(result): the path to the cluster closest to MFE.
The path is a list of length 10.
For i-th position, 1 represents that the predicted structure contains the i-th MEBP.
2 represents that predicted does not contain the i-th MEBP.
 0 represnts that clustering stops.

pathcorrectness: the path correctness to the cluster closest to MFE.
For i-th position, 1 represents correct, 0 represents wrong, -1 represents initial value.

pathent: the entropy of the cluster on the path. -1 represents initial value.

pathbpdist: bpdist between the cluster and MFE. -1 represents initial value.

bpset: the set of bias base pairs selected on the path.";
*)
Module[{i,idx={},index={},bp=BPSet[MFE],pos,result,MEBPSpace,str,ratio,leave,leavecorrectness=1,
pathcorrectness,pathent,pathbpdist,bpset={},currentcluster},


MEBPSpace=TreeMEBP[[1]][[TreeMEBP[[3]]]];


result=Table[0,{i,1,lvl}];

pathcorrectness=Table[-1,{i,1,lvl}];


pathent=Table[-1,{i,1,lvl+1}];


pathbpdist=Table[-1,{i,1,lvl+1}];


For[i=0,i<=lvl,i++,
idx=Union[idx,Tuples[{1,2},i]]];


pos=Position[idx,index][[1,1]] ;


For[i=1,i<=lvl,i++,

currentcluster=esb[[TreeMEBP[[4,pos]]]];
pathent[[i]]=StrEnt[currentcluster];
pathbpdist[[i]]=BpDistanceE[MFE,currentcluster];

If[TrueQ[MEBPSpace[[pos]]==List],

Break[],

bpset=Append[bpset,MEBPSpace[[pos]]];



If[MemberQ[bp,MEBPSpace[[pos]]]==True,

If[RandomInteger[{1,100}]>100*error,


AppendTo[index,1];
pathcorrectness[[i]]=1;
result[[i]]=1,


AppendTo[index,2];
pathcorrectness[[i]]=0;
leavecorrectness=0;
result[[i]]=2
],


If[RandomInteger[{1,100}]>100*error2,



AppendTo[index,2];
pathcorrectness[[i]]=1;
result[[i]]=2,


AppendTo[index,1];
pathcorrectness[[i]]=0;
leavecorrectness=0;
result[[i]]=1
]

]
];
pos=Position[idx,index][[1,1]] 

];
leave=esb[[TreeMEBP[[4,pos]]]];

pathent[[lvl+1]]=StrEnt[leave];
pathbpdist[[lvl+1]]=BpDistanceE[MFE,leave];

{str,ratio}=DomESB[leave];


Return[{str,ratio,{leave,leavecorrectness},{result, pathcorrectness,pathent,pathbpdist,bpset}}]
]


DrawTree[esb_,TreeMEBP_,lvl_]:=
(*
Input an ensemble of secondary structures using pair tables, and the tree space of indices of structures and the max height of the tree.
Return the drawing of the tree of clusters, each level of clusters representing the ones containing the maximum entropy base pairs or not.
*)
Module[{TreeSpace,idx={},idx1={},i,j,vertexset,edgeset,vertexshape,edgelabel={},TreeSpaceindex,MEBPSpace},
For[i=0,i<=lvl,i++,
idx=Union[idx,Tuples[{1,2},i]]];


For[i=0,i<=lvl-1,i++,
idx1=Union[idx1,Tuples[{1,2},i]]];

TreeSpaceindex=TreeMEBP[[4]][[1;;Length[idx]]];
TreeSpace=Table[esb[[TreeSpaceindex[[i]]]],{i,1,Length[idx]}];


vertexset=idx;

edgeset=Flatten[Table[{DirectedEdge[idx[[i]],Append[idx[[i]],1]],DirectedEdge[idx[[i]],Append[idx[[i]],2]]},{i,1,Length[idx1]}]];

vertexshape=Table[idx[[j]]->DrawESB[TreeSpace[[j]]],{j,1,Length[idx]}] ;

MEBPSpace=TreeMEBP[[1]][[TreeMEBP[[3]][[1;;Length[idx1]]]]];

If[Length[MEBPSpace]>0,
edgelabel=Flatten[Table[{edgeset[[2j-1]]->ToString[MEBPSpace[[j]]]<>" Yes",edgeset[[2j]]->"No"},{j,1,Length[idx1]}]]
,0
];
Return[TreeGraph[vertexset,edgeset,EdgeLabels->edgelabel,VertexShape->vertexshape,VertexSize->0.9]]
]



End[]


EndPackage[];
