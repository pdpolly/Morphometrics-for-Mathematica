(* ::Package:: *)

(*  This function prints the version number for the installed verison of this package.  *)

MorphometricsVersion[]:=Print["PollyMorphometrics 12.4\n(c) P. David Polly, 12 February 2022\n"];
MorphometricsVersion[];

Quiet[<<MultivariateStatistics`];


(* This functions does a Procrustes superimposition of landmark coordinates followed by orthogonal projection into tangent space.  The algorithm is the one 
    presented by Rohlf & Slice, 1990, Syst. Zool. 39: 40-59 and uses the tangent space projection from Rohlf, 1999, J. Class. 16: 197-223.

   Usage:    Procrustes[data, n, k] 
             where data are the landmark coodiantes to be aligned, n is the number of landmarks, and k is the number of dimensions of each landmark (2 or 3).
             Aligned coordinates are returned with each shape in a single row, with the columns being the x,y (,z) coordinates of the n landmarks.

	Updated 21 February 2010 to fix problem with 3D rotations.
	Updated 28 January 2012 to change output so that shapes retain original orientation.


*)

Procrustes[data_,nlandmarks_,ndims_]:=Block[{l,II,PP,SS,x,y,u,w,v,hh,ResidSS,NewResidSS,shape1},
l=Partition[Partition[Flatten[data],ndims],nlandmarks];
II=IdentityMatrix[nlandmarks];
PP=Table[N[1/nlandmarks],{nlandmarks},{nlandmarks}];
SS=Table[N[Sqrt[Tr[(II-PP) . l[[x]] . Transpose[l[[x]]] . (II-PP)]]],{x,Length[l]}];
l=Table[((II-PP) . l[[x]])/SS[[x]],{x,Length[l]}];
y=Mean[l];
ResidSS=Plus@@Flatten[Table[(Flatten[l[[x]]]-Flatten[y])^2,{x,Length[l]}]];
While[True,
	For[x=1,x<=Length[l],x++,
		{u,w,v}=SingularValueDecomposition[Transpose[l[[x]]] . y];
		hh=u . (w*Inverse[Abs[w]]) . Transpose[v];
		l[[x]]=l[[x]] . hh;

		];
	y=Mean[l];
    NewResidSS=Plus@@Flatten[Table[(Flatten[l[[x]]]-Flatten[y])^2,{x,Length[l]}]];
	If[Abs[NewResidSS-ResidSS]<0.0001,Break[]];
	ResidSS=NewResidSS;
];
l=Partition[PrincipalComponents[Partition[Flatten[l],ndims]],nlandmarks];
If[ndims==2,
	shape1=Partition[Flatten[data[[1]]],ndims];
	shape1=#-Mean[shape1]&/@shape1;
	shape1=shape1/Sqrt[Plus@@Plus@@(shape1^2)]//N;
	{u,w,v}=SingularValueDecomposition[Transpose[PrincipalComponents[shape1]] . shape1];
	hh=u . (w*Inverse[Abs[w]]) . Transpose[v];
	l=Partition[#,2] . hh&/@Partition[Flatten[l],ndims*nlandmarks];
];
l=Partition[Flatten[l],ndims* nlandmarks];
l=l . (IdentityMatrix[ndims *nlandmarks]-Mean[l] . Mean[l]);
Return[l]]


(* This functions does a Procrustes superimposition of landmark coordinates followed by orthogonal projection into tangent space but returns the aligned shapes 
   rescaled to their original sizes instead of to unit size.  The algorithm is modified from the one presented by Rohlf & Slice, 1990, Syst. Zool. 39: 40-59  
   and uses the tangent space projection from Rohlf, 1999, J. Class. 16: 197-223.

   Usage:    Procrustes[data, scaling, n, k] 
             where data are the landmark coodiantes to be aligned with each shape in a single row, scaling is the scaling factor, 
             n is the number of landmarks, and k is the number of dimensions of each landmark (2 or 3).   Aligned coordinates are returned 
             with each shape in a single row, with the columns being the x,y (,z) coordinates of the n landmarks.

	Updated 21 February 2010 to fix problem with 3D rotations.
	Updated 28 January 2012 to change output so that shapes retain original orientation.

*)

ProcrustesSized[data_,scalingfactors_,nlandmarks_,ndims_]:=Block[{l,II,PP,SS,x,y,u,w,v,hh,ResidSS,NewResidSS,shape1},
l=Partition[Flatten[data],ndims*nlandmarks];
l=l*scalingfactors;
l=Partition[Partition[Flatten[l],ndims],nlandmarks];
l=Table[Transpose[Transpose[l[[x]]]-Mean[l[[x]]]],{x,Length[l]}];

y=Mean[l];
ResidSS=Plus@@Flatten[Table[(Flatten[l[[x]]]-Flatten[y])^2,{x,Length[l]}]];
While[True,
	For[x=1,x<=Length[l],x++,
		{u,w,v}=SingularValueDecomposition[Transpose[y] . l[[x]]];
		hh=Transpose[v] . (w*Inverse[Abs[w]]) . u;
		l[[x]]=l[[x]] . hh;
		];
	y=Partition[Mean[Partition[Flatten[l],ndims *nlandmarks]],ndims];
    NewResidSS=Plus@@Flatten[Table[(Flatten[l[[x]]]-Flatten[y])^2,{x,Length[l]}]];
	If[Abs[NewResidSS-ResidSS]<0.01,Break[]];
	ResidSS=NewResidSS;
];
l=Partition[PrincipalComponents[Partition[Flatten[l],ndims]],nlandmarks];
If[ndims==2,
	shape1=Partition[Flatten[data[[1]]],ndims];
	shape1=#-Mean[shape1]&/@shape1;
	shape1=shape1/Sqrt[Plus@@Plus@@(shape1^2)]//N;
	{u,w,v}=SingularValueDecomposition[Transpose[PrincipalComponents[shape1]] . shape1];
	hh=u . (w*Inverse[Abs[w]]) . Transpose[v];
	l=Partition[#,2] . hh&/@Partition[Flatten[l],ndims*nlandmarks];
];
l=Partition[Flatten[l],ndims* nlandmarks];
l=l . (IdentityMatrix[ndims nlandmarks]-Mean[l] . Mean[l]);
Return[l]]


(* This functions returns centroid sizes of the several objects.  

   Usage:    CentroidSizes[data] 
             where data are the landmark coodiantes to be aligned, with each shape in a single row with the columns being the x,y (,z) 
             coordinates of the n landmarks and with the first column containing the scaling factor.
*)

CentroidSizes[data_]:=Block[{sized,centroids,sizes},
sized=data[[1;;,2;;]]*data[[1;;,1]];
centroids={Mean[#[[1;;;;2]]],Mean[#[[2;;;;2]]]}&/@sized;
sizes=Table[Sqrt[Plus@@Flatten[(#-centroids[[x]]&/@Partition[sized[[x]],2])^2]],{x,Length[sized]}];
sizes
]


(* This function creates a thin plate spline graphic for the deformation of 2D landmark shapes.  The algorithm is from Dryden & Mardia 1998, and Hammer & Harper, 2006. 

   Usage:   tpSpline[source, target (, label, gridscale)]
            where source and target are sets of aligned landmark coordinates.  The results shows the deformation of target from source.  Colour and size of primitives in 
            the plot can be changed by modifying the final Return[] line options.

	Updates:   

		13 Aug 2014: fixed problem in equation for calculating A matrix in which it was not transposed in previous versions.
*)
tpSpline[source_,target_,label_:"none",grdscl_:20]:=Block[{T,Y,st,SigmaH,i,j,InvBigG,W,c,A,LastTerm,yyyy,xmax,xmin,ymax,ymin,xint,yint,plotsize,gridpts,grid,dots},
T=Partition[Flatten[source],2];
Y=Partition[Flatten[target],2];
SigmaH[vec1_,vec2_]:=Module[{h},If[(h=Norm[vec1-vec2])>0,Return[(h^2 )*Log[h]],Return[0]]];
st=Table[SigmaH[T[[i]],T[[j]]],{i,Length[T]},{j,Length[T]}];
InvBigG=Inverse[ArrayFlatten[{{st,1,T},{1,0,0},{Transpose[T],0,0}}]];
W=InvBigG[[1;;Length[T],1;;Length[T]]] . Y;
{c,A}={(InvBigG[[Length[T]+1;;Length[T]+3,1;;Length[T]]] . Y)[[1,1;;2]], Transpose[(InvBigG[[Length[T]+1;;Length[T]+3,1;;Length[T]]] . Y)[[2;;,1;;2]]]};
LastTerm[pt_]:=Transpose[W] . Table[SigmaH[pt,T[[i]]],{i,Length[T]}];
yyyy[pt_]:=Flatten[c]+A . pt+LastTerm[pt];
xmax=2 Max[Transpose[T][[1]]];
xmin=2 Min[Transpose[T][[1]]];
ymax=2 Max[Transpose[T][[2]]];
ymin=2 Min[Transpose[T][[2]]];
xint=((xmax-xmin) ymax)/(grdscl xmax);
yint=(ymax-ymin)/grdscl;
plotsize=Mean[{xmax-xmin,ymax-ymin}];
gridpts=Table[yyyy[{x,y}],{x,xmin,xmax,xint},{y,ymin,ymax,yint}];
grid=Flatten[{Table[Line[{gridpts[[i,j]],gridpts[[i+1,j]]}],{i,Length[gridpts]-1},{j,Length[gridpts[[1]]]}],Table[Line[{gridpts[[i,j]],gridpts[[i,j+1]]}],{i,Length[gridpts]},{j,Length[gridpts[[1]]]-1}]}];
dots=Table[Point[Flatten[{yyyy[T[[x]]]}]],{x,Length[T]}];
If[label!="none",
Return[Show[Graphics[{Gray,grid}],Graphics[{PointSize[0.05],dots}],AspectRatio->Automatic,Frame->False,PlotLabel->label]],
Return[Show[Graphics[{Gray,grid}],Graphics[{PointSize[0.05],dots}],AspectRatio->Automatic,Frame->False]]
];

]


(*  This function rotates one set of Procrustes aligned coordinates to the same orientation as another to make them comparable in case they have been aligned separately.
    
    Usage:  CommonOrientation[data1, data2, n, k]
            where data1 are aligned coordinates with the base orientation and data2 are another set of aligned coordinates that will be rotated into the base alignment.
            n is the number of landmarks and k is the number of dimensions for each coordinate (2 or 3). The function returns data2 rotated into the base alignment with 
            each specimen on a row.
*)

CommonOrientation[baseline_,data_,landmarks_,dims_]:=Block[{tempdata,Consensus1,Consensus2,u,w,v,h, hh},
Consensus1=Partition[Mean[baseline],dims];
Consensus2=Partition[Mean[data],dims];
{u,w,v}=SingularValueDecomposition[Transpose[Consensus1] . Consensus2];
hh=Transpose[v] . DiagonalMatrix[Diagonal[w]/Abs[Diagonal[w]]] . u;
tempdata=Partition[Partition[Flatten[data],dims],landmarks];
tempdata=Table[tempdata[[x]] . h,{x,Length[tempdata]}];
Return[Partition[Flatten[tempdata],dims landmarks]];]



(*  This function reads in a TPS file such as the ones used by the tps-series of programs by F. James Rohlf.  

    Usage:  tpsImport[filename (, labelsource)]
            where filename is the name of the file to be imported and labelsource is an optional parameter specifying what identifier should be used to lable the objects
			(the default is the ID tag, option is the IMAGE tag or a list of lables).  If the file is not in the working directory a full path may need to be included with the 
            filename.  If the TPS file contains landmarks but no outlines then the function returns landmark coordinates with each specimen in a single row.  The first column contains 
			a specimen label that is based on the file name in the ID tag of the tps file unless "IMAGE" is specified as an option.  If a scale has been set, as with tpsDIG, the second 
			column contains the scaling factor.  Subsequent columns contain the x,y (,z) coordinates of the landmarks.  If the TPS file contains outlines, then these are returned in a
			second data block (multiple blocks if there is more than one outline or curve per object).  Each outline block has object labels in the first column and scale factors in the 
			second column if appropriate.  

	Created: 11 April 2014
	Updated: 22 April 2014 to fix problem with transposing outlines when there are no outlines.
*)


tpsImport[filename_,lab_:"ID"]:=Block[{coords,LMlabels,LMnumber,OUTLINESlabels,OUTs,ObjectBreakPoints,objects,LMs,ThisObject,OUTLINESnumber,IMAGElabels,SCALElabels,IDlabels,landmarks,outlines,POINTSlabels,POINTSnumber,OutlineStart,OutlineBlock,obj,labellist},

(* First read in the file and determine the number of objects and whether there are outlines. Complain to user if the number of landmarks or outlines differ between objects. *)
coords=Import[filename,"Table"];
LMlabels=Select[coords,If[StringQ[#[[1]]],StringMatchQ[#[[1]],"LM=*"]]&];
LMnumber=ToExpression[StringReplace[Flatten[LMlabels],"LM="->""]];
If[Length[Union[LMnumber]]!=1,Return["Objects have different number of landmarks."]];
OUTLINESlabels=Select[coords,If[StringQ[#[[1]]],StringMatchQ[#[[1]],"OUTLINES=*"]]&];
If[Length[Union[Flatten[OUTLINESlabels]]]>1,Return["Objects have different number of outlines."]];
If[Length[Union[Flatten[OUTLINESlabels]]]>0,OUTs={}];

(* break the objects into blocks and being a landmark output block. *)
ObjectBreakPoints=Append[Flatten[Position[coords,LMlabels[[1]]]],Length[coords]+1];
objects=coords[[#]]&/@Table[ObjectBreakPoints[[x]];;ObjectBreakPoints[[x+1]]-1,{x,Length[ObjectBreakPoints]-1}];
LMs={};

(* create the label list *)
IDlabels=StringReplace[Flatten[Select[coords,If[StringQ[#[[1]]],StringMatchQ[#[[1]],"ID=*"]]&]],"ID="->""];
IMAGElabels=Select[coords,If[StringQ[#[[1]]],StringMatchQ[#[[1]],"IMAGE=*"]]&];
IMAGElabels=StringReplace[Flatten[IMAGElabels],"IMAGE="->""];
If[lab=="IMAGE",labellist=IMAGElabels];
If[Length[lab]==Length[objects],labellist=lab];
If[Length[labellist]!=Length[objects],labellist=IDlabels];

(* process each object *)
Do[
ThisObject=objects[[obj]];
(* get interior labels *)
OUTLINESlabels=Flatten[Select[ThisObject,If[StringQ[#[[1]]],StringMatchQ[#[[1]],"OUTLINES=*"]]&]];
OUTLINESnumber=ToExpression[StringReplace[Flatten[OUTLINESlabels],"OUTLINES="->""]];
SCALElabels=Select[ThisObject,If[StringQ[#[[1]]],StringMatchQ[#[[1]],"SCALE=*"]]&];
SCALElabels=ToExpression[StringReplace[Flatten[SCALElabels],"SCALE="->""]];
(* process the landmarks *)
If[LMnumber[[obj]]!=0,landmarks=Flatten[ThisObject[[2;;LMnumber[[obj]]+1]]],landmarks="No landmarks"];
If[LMnumber[[obj]]!=0&&Length[SCALElabels]>0,landmarks=Prepend[landmarks,SCALElabels[[1]]]];
landmarks=Prepend[landmarks,labellist[[obj]]];
LMs=Append[LMs,landmarks];

(* process the outlines, if any. *) 
If[Length[OUTLINESlabels]>0,
POINTSlabels=Select[ThisObject,If[StringQ[#[[1]]],StringMatchQ[#[[1]],"POINTS=*"]]&];
POINTSnumber=ToExpression[StringReplace[Flatten[POINTSlabels],"POINTS="->""]];
OutlineStart=Flatten[Position[ThisObject,OUTLINESlabels[[1]]]];
OutlineBlock=ThisObject[[OutlineStart[[1]]+1;;(Plus@@POINTSnumber+Length[POINTSnumber]+OutlineStart[[1]])]];
outlines=OutlineBlock[[#]]&/@Table[(Plus@@(POINTSnumber[[1;;x]]+1))-POINTSnumber[[x]]+1;;Plus@@(POINTSnumber[[1;;x]]+1),{x,Length[POINTSnumber]}];
If[Length[SCALElabels]>0,outlines=Prepend[outlines,SCALElabels[[1]]]];
outlines=Prepend[outlines,labellist[[obj]]];
OUTs=Append[OUTs,Flatten[outlines]];
];
,{obj,Length[objects]}];

If[Length[OUTs]==1,OUTs=OUTs[[1]]];
(* return landmarks and outlines in two blocks if appropriate, else just return landmarks *)
If[Length[OUTs]>0, Return[{LMs,OUTs}],Return[LMs]];
]





(*  This function reads in a basic TPS file such as the ones used by the tps-series of programs by F. James Rohlf.  

    Usage:  tpsImportOrig[filename, dimensions]
            where filename is the name of the file to be imported and dimensions is an optional parameter specifying whether the landmarks are two or three dimensional
			(the default is two).  If the file is not in the working directory a full path may need to be included with the 
            filename.  The function returns data with each specimen in a single row.  The first column contains a specimen label that is based on the file name in the
            IMAGE tag of the tps file, or the ID tag if there is no IMAGE tag.  If a scale has been set, as with tpsDIG, the second column contains the scaling factor.  
			Subsequent columns contain the x,y (,z) coordinates of the landmarks.  

			The function is simple and doesn't expect any lines in the tps file other than LM=, IMAGE=, SCALE=, and ID=.  If your file doesn't read in correctly, check 
			whether it has additional lines. 

	Updated 28 January 2012 to fix incompatibilities between PC and Mac line breaks, to allow for comma or tab separtion of landmarks,
                             and to take labels from the ID tag if there is no IMAGE tag.
*)

tpsImportOrig[filename_,ndims_:2]:=Block[{tps,x,i},
tps=ReadList[filename,Word,RecordLists->True,WordSeparators->{" ",",","\n","\t", "\r"},RecordSeparators->{"LM="}];

Do[
tps[[x]]=StringReplace[tps[[x]],"E"~~i:DigitCharacter->"*10^"<>ToString[i]];

If[Length[Select[tps[[x]],StringMatchQ[#,"IMAGE=*"]&]]>0,
If[Length[Select[tps[[x]],StringMatchQ[#,"SCALE=*"]&]]==0,
(* For TPS with both IMAGE but no SCALE tags *)
tps[[x]]=Flatten[{StringReplace[StringReplace[Select[tps[[x]],StringMatchQ[#,"IMAGE=*"]&],"IMAGE="->""],".JPG"->"",IgnoreCase->True],ToExpression[Take[tps[[x]],{2,1+ndims ToExpression[tps[[x,1]]]}]]}],

(* For TPS with both IMAGE and SCALE tags *)

tps[[x]]=Flatten[{StringReplace[StringReplace[Select[tps[[x]],StringMatchQ[#,"IMAGE=*"]&],"IMAGE="->""],".JPG"->"",IgnoreCase->True],ToExpression[StringReplace[Select[tps[[x]],StringMatchQ[#,"SCALE=*"]&],"SCALE="->""]],ToExpression[Take[tps[[x]],{2,1+ndims ToExpression[tps[[x,1]]]}]]}]],

(* For TPS with neither IMAGE nor SCALE tags. *)

tps[[x]]=Flatten[{StringReplace[Select[tps[[x]],StringMatchQ[#,"ID=*"]&],"ID="->""],ToExpression[Take[tps[[x]],{2,1+ndims *ToExpression[tps[[x,1]]]}]]}]
],



{x,Length[tps]}];

Return[tps];

]



(*  This function returns the Procrustes distance between two shapes after it places them in optimal Procrustes alignment.  It requires the Procrustes[] function from 
    this module to be available.  

    Usage:  ProcrustesDistance[shape1, shape2, k]
            where shape1 and shape2 are the two shapes whose distance is to be calculated and k is the number of dimensions of each coordinate (2 or 3).
*)

ProcrustesDistance[shape1_,shape2_,dims_]:=Block[{superd},
superd=Procrustes[{shape1,shape2},Length[shape1]/dims,dims];
Sqrt[Plus@@((superd[[1]]-superd[[2]])^2)]
]


(*  This function performs a Mantel Test on two matrices.  It returns the matrix correlation and the probability (p) that the two have the same structure.
.  The function incorporates the recommendations of Manly, 1986, Res. Populat. Ecol., 28:201-281.

    Usage:  Mantel[matrix1, matrix2, n]
            where matrix1 and matrix2 are the two symmetric matrices whose Mantel correlation is to be calculated and n is the number of randomizations to be performed to 
            compute the test statistics.
*)

Mantel[matrix1_,matrix2_,n_]:=Block[{M1,M2,i,j,z,znull,R1,R2,rz,p},
M1=Flatten[Table[Table[matrix1[[i,j]],{j,i+1,Length[matrix1]}],{i,Length[matrix1-1]}]];
M1=(#-Mean[M1])/StandardDeviation[M1]&/@M1;
M2=Flatten[Table[Table[matrix2[[i,j]],{j,i+1,Length[matrix2]}],{i,Length[matrix2-1]}]];
M2=(#-Mean[M2])/StandardDeviation[M2]&/@M2;
z=Plus@@(M1*M2)/(Length[M1]-1);
znull=Table[R1=RandomSample[M1];R2=RandomSample[M2];
rz=(Plus@@(R1*R2))/(Length[R1]-1),{n}];
p=Length[Select[znull,#>z&]]/n//N;
Return[{z,p}]
]


(*  This function performs a Mantel Test on two matrices.  It returns the test statistics with helpful labels.  The function incorporates the recommendations 
    of Manly, 1986, Res. Populat. Ecol., 28:201-281. As is appropriate for geometric morphometrics, the diagonal elements are included
    in the matrix correlation (Klingenberg and McIntyre, 1998). 

    Usage:  MantelForMorphometrics[matrix1, matrix2, n]
            where matrix1 and matrix2 are the two symmetric matrices whose Mantel correlation is to be calculated and n is the number of randomizations to be performed to 
            compute the test statistics.
            
    Created: 5 March 2017.
*)

MantelForMorphometrics[matrix1_,matrix2_,n_]:=Block[{M1,M2,i,j,z,znull,R1,R2,rz,p},
M1=Flatten[Table[Table[matrix1[[i,j]],{j,i,Length[matrix1]}],{i,Length[matrix1]}]];
M1=(#-Mean[M1])/StandardDeviation[M1]&/@M1;
M2=Flatten[Table[Table[matrix2[[i,j]],{j,i,Length[matrix2]}],{i,Length[matrix2]}]];
M2=(#-Mean[M2])/StandardDeviation[M2]&/@M2;
z=Plus@@(M1*M2)/(Length[M1]-1);
znull=Table[R1=RandomSample[M1];R2=RandomSample[M2];
rz=(Plus@@(R1*R2))/(Length[R1]-1),{n}];
p=Length[Select[znull,#>z&]]/n//N;
Return[{z,p}]
]


(* This function does a complete Principal Components analysis of shape on a set of procrustes superimposed landmarks.  

   Usage:    PrincipalComponentsOfShape[data,{PC1,PC2}, labels] 
             where data are Procrustes aligned shapes, PC1 and PC2 are two integers specifying which two principal components should be displayed in results, and labels is a list of specimen labels whose order corresponds to the Procrustes aligned shapes.
			 A Principal Components plot is returned, followed by a summary of the variance explained by the PCs in the plot, followed by a graphic of the mean shape with its landmarks numbered, followed by a morphospace diagram in which representative thin-plate spline diagrams have been placed to show the variation of shape through the PC space.             

	Created 28 January 2012.

*)


PrincipalComponentsOfShape[proc_,PCs_,labels_]:=Block[{marginfactor,MShape,Resids,P,u,w,v,EValues,PropEValues,scores,xrange,yrange,pts,labs,InsScale},
marginfactor=1.5;
MShape=Mean[proc];
Resids=#-MShape&/@proc;
P=Covariance[Resids];
{u,w,v}=SingularValueDecomposition[P];
EValues=Tr[w,List];
PropEValues=EValues/Plus@@EValues;
scores=Resids . u;
u=Transpose[u];
xrange=marginfactor*{Min[scores[[1;;,PCs[[1]]]]],Max[scores[[1;;,PCs[[1]]]]]};
yrange=marginfactor*{Min[scores[[1;;,PCs[[2]]]]],Max[scores[[1;;,PCs[[2]]]]]};
InsScale=.4*Mean[{Plus@@Abs[xrange],Plus@@Abs[yrange]}];
pts={PointSize[0.015],Point[#]}&/@scores[[1;;,PCs]];
labs=Table[{FontFamily-> "Helvetica",Text[labels[[x]],scores[[x,PCs]],{0,1}]},{x,Length[labels]}];

Print[Graphics[{pts,labs},Frame->True,AspectRatio->Automatic,FrameLabel->{"PC "<>ToString[PCs[[1]]], "PC "<>ToString[PCs[[2]]]},BaseStyle-> {FontFamily->"Helvetica"},PlotRange->{xrange,yrange}]];

Print["PC "<>ToString[PCs[[1]]]<>" explains "<>ToString[PropEValues[[PCs[[1]]]]]<>" of total shape variance"];
Print["PC "<>ToString[PCs[[2]]]<>" explains "<>ToString[PropEValues[[PCs[[2]]]]]<>" of total shape variance"];

Print[PlanarGraphPlot[Partition[Mean[proc],2],ConvexHull[Partition[Mean[proc],2]],BaseStyle->{FontFamily->"Helvetica",FontSize->14},PlotLabel-> "Mean Shape",AspectRatio->Automatic,ImageSize->250]];

Print[Graphics[{},Frame->True,AspectRatio->Automatic,PlotLabel->"Morphospace",FrameLabel->{"PC "<>ToString[PCs[[1]]], "PC "<>ToString[PCs[[2]]]},BaseStyle-> {FontFamily->"Helvetica",FontSize->12},PlotRange->{1.5*xrange,1.5*yrange},

Epilog-> {
Inset[tpSpline[MShape,xrange[[1]]*u[[PCs[[1]]]]+yrange[[1]]*u[[PCs[[2]]]]+MShape],{xrange[[1]],yrange[[1]]},Center,InsScale*{1,1}],
Inset[tpSpline[MShape,yrange[[1]]*u[[PCs[[2]]]]+MShape],{0,yrange[[1]]},Center,InsScale*{1,1}],
Inset[tpSpline[MShape,xrange[[2]]*u[[PCs[[1]]]]+yrange[[1]]*u[[PCs[[2]]]]+MShape],{xrange[[2]],yrange[[1]]},Center,InsScale*{1,1}],

Inset[tpSpline[MShape,xrange[[1]]*u[[PCs[[1]]]]+yrange[[2]]*u[[PCs[[2]]]]+MShape],{xrange[[1]],yrange[[2]]},Center,InsScale*{1,1}],
Inset[tpSpline[MShape,yrange[[2]]*u[[PCs[[2]]]]+MShape],{0,yrange[[2]]},Center,InsScale*{1,1}],
Inset[tpSpline[MShape,xrange[[2]]*u[[PCs[[1]]]]+yrange[[2]]*u[[PCs[[2]]]]+MShape],{xrange[[2]],yrange[[2]]},Center,InsScale*{1,1}],

Inset[tpSpline[MShape,xrange[[1]]*u[[PCs[[1]]]]+MShape],{xrange[[1]],0},Center,InsScale*{1,1}],
Inset[tpSpline[MShape,MShape],{0,0},Center,InsScale*{1,1}],
Inset[tpSpline[MShape,xrange[[2]]*u[[PCs[[1]]]]+MShape],{xrange[[2]],0},Center,InsScale*{1,1}]

}]];

Return[];

]




(* This function provides a diagnostic plot of Procrustes superimposed shapes. 

   Usage:    ProcrustesPlot[data] 
			 where data are Procrustes aligned shapes.
             A plot of the superimposed landmarks is returned in which each landmark has its own color and the landmark number is placed at the landmark's centroid.

	Created 28 January 2012.

*)

ProcrustesPlot[proc_]:=Block[{PlotPts,MShape,Labels,x},
PlotPts=Partition[#,2]&/@proc;
MShape=Mean[PlotPts];
Labels=Table[x,{x,Length[MShape]}];
Return[
Show[{
ListPlot[Transpose[PlotPts],Axes->False,AspectRatio->Automatic,ImageSize->500],
Graphics[{FontFamily->"Helvetica",FontSize->14,Table[Text[Labels[[x]],MShape[[x]]],{x,Length[MShape]}]},AspectRatio->Automatic]
}]
];
]


(* This function flips the y coordinates of a file that was created in ImageJ so that shapes will be right-side up.  
   
   Usage:  FlipImageJ[data]
           where data is a rectangular matrix with the first column containing labels and the second through final 
           columns containing the x and y coordinates of several 2-d landmarks.

	Created 28 January 2012.
*)

FlipImageJ[data_]:=Block[{},
Return[Partition[Flatten[Riffle[data[[1;;,1]],Partition[Flatten[({1,-1}*#&/@Partition[#,2])&/@data[[1;;,2;;]]],Length[data[[1,2;;]]]]]],Length[data[[1]]]]]
]


(*  This function rotates a set of 2D landmark shapes by a specified number of degrees.

    Usage:  RotateShape[data, deg]
            where data is a matrix of 2D shapes with coordinates in the columns and objects in the rows and deg is the number
            of degrees counterclockwise to rotate the shape.


	Created 5 February 2012.
*)

RotateShape[data_,deg_]:=Block[{},
Return[Partition[Flatten[Transpose[RotationMatrix[deg Degree] . Transpose[Partition[#,2]]]&/@data],Length[data[[1]]]]]
]


(*  This function creates a diagnostic plot of Procrustes superimposed 3D landmarks.  Each landmark is plotted in a different color 
    with the landmark number at its centroid.  Solid grey lines connect landmarks to their landmark centroid, and dotted lines connect
    each landmark centroid to the object's centroid.  

    Usage:   ProcrustesPlot3D[data]
             where data is a matrix of 3D shapes with coordinates in the columns and objects in the rows.


	Created 5 February 2012.
*)

ProcrustesPlot3D[proc_]:=Block[{lands,pts,x,labs,lns,lns2},
lands=Transpose[Partition[Partition[Flatten[proc],3],Length[proc[[1]]]/3]];
pts=Table[{PointSize[0.015],ColorData["Rainbow",x/Length[lands]],Point[#]}&/@lands[[x]],{x,Length[lands]}];
labs=Table[{FontFamily->"Helvetica",FontSize->14,Text[x,Mean[lands[[x]]]]},{x,Length[lands]}];
lns={Thickness[0.005],Dotted,LightGray,Line[{{0,0,0},#}]}&/@Mean[Transpose[lands]];
lns2=Table[{Thickness[0.001],LightGray,Line[{Mean[lands[[x]]],#}]}&/@lands[[x]],{x,Length[lands]}];
Return[Graphics3D[{pts,labs,lns,lns2},Boxed->False]];
]


(* This function does a complete Principal Components analysis of shape on a set of 3D procrustes superimposed landmarks.  

   Usage:    PrincipalComponentsOfShape3D[data,{PC1,PC2}, labels] 
             where data are Procrustes aligned 3D shapes, PC1 and PC2 are two integers specifying which two principal 
             components should be displayed in results, and labels is a list of specimen labels whose order corresponds 
             to the Procrustes aligned shapes.
			 A Principal Components plot is returned, followed by a summary of the variance explained by the PCs in the 
             plot, followed by a graphic of the mean shape with its landmarks numbered, followed by two shape models, one
             for each PC.  Each shape model has landmarks numbered at their centroids with vectors passing through them. 
             The tails of the vectors indicate the positions of the landmarks at the negative end of the PC and the heads
             indicate the positions of the landmarks at the positive end of the PC.

	Created 28 January 2012.

*)

PrincipalComponentsOfShape3D[proc_,PCs_,labels_]:=Block[{MShape,Resids,P,u,w,v,EValues,PropEValues,scores,xrange,yrange,InsScale,pts,labs,x,y,lands,lns,PCxPts,PCxLns,EndPtsPCx,MnTxtPCx,PCyPts,PCyLns,EndPtsPCy,MnTxtPCy,marginfactor},
marginfactor=1.5;

MShape=Mean[proc];
Resids=#-MShape&/@proc;
P=Covariance[Resids];
{u,w,v}=SingularValueDecomposition[P];
EValues=Tr[w,List];
PropEValues=EValues/Plus@@EValues;
scores=Resids . u;
u=Transpose[u];
xrange=marginfactor*{Min[scores[[1;;,PCs[[1]]]]],Max[scores[[1;;,PCs[[1]]]]]};
yrange=marginfactor*{Min[scores[[1;;,PCs[[2]]]]],Max[scores[[1;;,PCs[[2]]]]]};
InsScale=.4*Mean[{Plus@@Abs[xrange],Plus@@Abs[yrange]}];

pts={PointSize[0.015],Point[#]}&/@scores[[1;;,PCs]];
labs=Table[{FontFamily-> "Helvetica",Text[labels[[x]],scores[[x,PCs]],{0,1}]},{x,Length[labels]}];

Print[Graphics[{pts,labs},Frame->True,AspectRatio->Automatic,FrameLabel->{"PC "<>ToString[PCs[[1]]], "PC "<>ToString[PCs[[2]]]},BaseStyle-> {FontFamily->"Helvetica"},PlotRange->{xrange,yrange}]];


Print["PC "<>ToString[PCs[[1]]]<>" explains "<>ToString[PropEValues[[PCs[[1]]]]]<>" of total shape variance"];
Print["PC "<>ToString[PCs[[2]]]<>" explains "<>ToString[PropEValues[[PCs[[2]]]]]<>" of total shape variance"];

lands=Transpose[Partition[Partition[Flatten[proc],3],Length[proc[[1]]]/3]];
labs=Table[{FontFamily->"Helvetica",FontSize->14,Text[x,Mean[lands[[x]]]]},{x,Length[lands]}];
lns={Thickness[0.005],Dotted,LightGray,Line[{{0,0,0},#}]}&/@Partition[MShape,3];
Print[Graphics3D[{labs,lns},PlotLabel->"Mean Shape",BaseStyle->{FontFamily->"Helvetica"}]];

PCxPts=Table[Partition[x*u[[PCs[[1]]]]+MShape,3],{x,xrange[[1]],xrange[[2]],(xrange[[2]]-xrange[[1]])/10}];
PCxLns={LightGray,Thickness[0.002],Arrow[{#[[1]],#[[-1]]}]}&/@Transpose[PCxPts];
EndPtsPCx=Table[{PointSize[0.02],ColorData["Rainbow",x/(Length[Transpose[PCxPts]])],Point[Transpose[PCxPts][[x,1]]],Point[Transpose[PCxPts][[x,-1]]]},{x,Length[Transpose[PCxPts]]}];
MnTxtPCx=Table[{FontFamily->"Helvetica",FontSize->14,Text[x,Partition[MShape,3][[x]]]},{x,Length[MShape]/3}];
Print[Graphics3D[{PCxLns,MnTxtPCx,EndPtsPCx},Boxed->False,BaseStyle->{FontFamily->"Helvetica"},PlotLabel->"PC "<>ToString[PCs[[1]]]<>" Shape Model"]];

PCyPts=Table[Partition[y*u[[PCs[[2]]]]+MShape,3],{y,yrange[[1]],yrange[[2]],(yrange[[2]]-yrange[[1]])/10}];
PCyLns={LightGray,Thickness[0.002],Arrow[{#[[1]],#[[-1]]}]}&/@Transpose[PCyPts];
EndPtsPCy=Table[{PointSize[0.02],ColorData["Rainbow",y/(Length[Transpose[PCyPts]])],Point[Transpose[PCyPts][[y,1]]],Point[Transpose[PCyPts][[y,-1]]]},{y,Length[Transpose[PCyPts]]}];
MnTxtPCy=Table[{FontFamily->"Helvetica",FontSize->14,Text[y,Partition[MShape,3][[y]]]},{y,Length[MShape]/3}];
Print[Graphics3D[{PCyLns,MnTxtPCy,EndPtsPCy},Boxed->False,BaseStyle->{FontFamily->"Helvetica"},PlotLabel->"PC "<>ToString[PCs[[2]]]<>" Shape Model"]];

]


(*  This function uses polynomial interpolation to find n equally spaced points from 
	an input curve of x - y coordinates.  Note that it only works on one object at a 
	time.  Point coordinates are returned as a simple list.

    Usage:     newpts = EqualSpace[oldpoints, 100]
               where oldpoints is the array containing x-y coordinates and 100 is 
               the number of new equally spaced points to construct.  Optional 

    Created:   1 March 2012
	Updated:   1 April 2014
			  12 April 2014

*)

EqualSpace[curve_,n_,dim_:2,interporder_:1]:=Block[{newcurve,dists,pos,f},
newcurve=Partition[Flatten[curve],dim];
dists=Prepend[Table[Sqrt[Plus@@((newcurve[[x]]-newcurve[[x-1]])^2)]//N,{x,2,Length[newcurve]}],0];
pos=Accumulate[dists];
f=Interpolation[Transpose[{pos,#}],InterpolationOrder->interporder]&/@Transpose[newcurve];
newcurve=Table[#[x]&/@f,{x,0,pos[[-1]],pos[[-1]]/(n-1)}];
dists=Prepend[Table[Sqrt[Plus@@((newcurve[[x]]-newcurve[[x-1]])^2)]//N,{x,2,Length[newcurve]}],0];
pos=Accumulate[dists];
f=Interpolation[Transpose[{pos,#}],InterpolationOrder->interporder]&/@Transpose[newcurve];
newcurve=Table[#[x]&/@f,{x,0,pos[[-1]],pos[[-1]]/(n-1)}];
If[(Length[Flatten[newcurve]])!=(n*dim),Print["Warning: "<>ToString[Length[Flatten[newcurve]]/dim]<>" points returned instead of "<>ToString[n]<>"."]];
Return[Flatten[newcurve]];
]




(*  This function imports outline coordinates created using the polygon tool and 
    save XY coordinates function of ImageJ.  It reads in all files that match the 
    filter, assuming that each one is a plain text file.  The file name is prepended
    to each row of coordinates as a label.

    Usage:     outlinecoords = ImportImageJPoints[path, filter]
               where path is the path to the directory where the x-y 
               coordinates are stored and filter is of the form "*.txt".

    Created:  1 March 2012.

*)

ImportImageJPoints[path_,filter_]:=Block[{files,coordinates,x},
SetDirectory[path];
files=FileNames[filter];
coordinates=Import[#,"Table"]&/@files;
ResetDirectory[];
coordinates=Table[Prepend[coordinates[[x, 2;;, {2,3}]],files[[x]]],{x,Length[coordinates]}];
coordinates=Flatten[#]&/@coordinates;
Return[coordinates];
]



(*  This function imports outline coordinates created using the polygon tool and 
    save XY coordinates function of ImageJ.  It reads in all files that match the 
    filter, assuming that each one is a plain text file.  The file name is prepended
    to each row of coordinates as a label.

    Usage:     outlinecoords = ImportImageJOutlines[path, filter]
               where path is the path to the directory where the x-y 
               coordinates are stored and filter is of the form "*.txt".

    Created:  1 March 2012.

*)

ImportImageJOutlines[path_,filter_]:=Block[{files,coordinates,x},
SetDirectory[path];
files=FileNames[filter];
coordinates=Import[#,"Table"]&/@files;
ResetDirectory[];
coordinates=Table[Prepend[coordinates[[x]],files[[x]]],{x,Length[coordinates]}];
coordinates=Flatten[#]&/@coordinates;
Return[coordinates];
]



(*  This function does a Euclidean Distance Matrix Analysis (EDMA) on two samples of 
    two-dimensional landmark coordinates.  Normally the coordinates are scaled in 
     real units and have not been Procrustes superimposed.  The function returns two 
    graphs, the first showing the values of the Form Distance Matrix of sample 2 compared
    to sample 1, and the second showing the mean shape of sample 2 with the lengthened
    interlandmark distances colored deep red and the shortened distance colored blue.

    Usage:     EDMA[landmarks1, landmarks2, landmarklabels]
               where landmarks1 and landmarks2 are two samples of 2D landmarks and
               landmarklabels is a list of labels for the landmarks.


    Created:  24 October 2012.

*)


EDMA[Sample1_,Sample2_,LandmarkLabels_]:=Block[{k,E1,E2,Ebar1,Ebar2,Var1,Var2,Delta1,Delta2,FDM,SortedFDM,SortedLabelPairs,
SD,lowlines,highlines,LabelPairs,consensus2,points,lines,j,i,x,samp1,samp2, boot1, boot2, bootE1,bootE2,
bootEbar1,bootEbar2,bootVar1,bootVar2,bootDelta1,bootDelta2,bootFDM,bootFDMmax,bootFDMmin,SortedFDMmin,SortedFDMmax},
k=Length[Flatten[Sample1[[1]]]]/Length[LandmarkLabels];
samp1=Partition[Partition[Flatten[Sample1],k],Length[LandmarkLabels]];
samp2=Partition[Partition[Flatten[Sample2],k],Length[LandmarkLabels]];
E1=Flatten[Table[Table[Sqrt[Plus@@((#[[i]]-#[[j]])^2)],{i,j+1,Length[#],1}],{j,Length[#]-1}]]&/@samp1;
E2=Flatten[Table[Table[Sqrt[Plus@@((#[[i]]-#[[j]])^2)],{i,j+1,Length[#],1}],{j,Length[#]-1}]]&/@samp2;
Ebar1=Mean[E1];
Ebar2=Mean[E2];
Var1=Variance[E1];
Var2=Variance[E2];
If[k==2,Delta1=(Ebar1-Var1)^0.25;Delta2=(Ebar2-Var2)^0.25, Delta1=(Ebar1-(1.5*Var1))^0.25;Delta2=(Ebar2-(1.5*Var2))^0.25];
FDM=Delta2/Delta1//N;
LabelPairs=Flatten[Table[Table[LandmarkLabels[[i]]<>"-"<>LandmarkLabels[[j]],{i,j+1,Length[LandmarkLabels]}],{j,Length[LandmarkLabels]-1}]];
SD=StandardDeviation[FDM];

(* Next section bootstraps the samples for 95% confidence intervals *)
bootFDM=Table[Null,{1000}];
Do[
boot1=samp1[[RandomInteger[{1,Length[samp1]},Length[samp1]]]];
boot2=samp2[[RandomInteger[{1,Length[samp2]},Length[samp2]]]];
bootE1=Flatten[Table[Table[Sqrt[Plus@@((#[[i]]-#[[j]])^2)],{i,j+1,Length[#],1}],{j,Length[#]-1}]]&/@boot1;
bootE2=Flatten[Table[Table[Sqrt[Plus@@((#[[i]]-#[[j]])^2)],{i,j+1,Length[#],1}],{j,Length[#]-1}]]&/@boot2;
bootEbar1=Mean[bootE1];
bootEbar2=Mean[bootE2];
bootVar1=Variance[bootE1];
bootVar2=Variance[bootE2];
bootDelta1=(Clip[bootEbar1-bootVar1,{0.0001,Infinity}])^0.25;
bootDelta2=(Clip[bootEbar2-bootVar2,{0.0001,Infinity}])^0.25;
bootFDM[[x]]=bootDelta2/bootDelta1//N,
{x,1000}];

bootFDMmax=Sort[#][[950]]&/@Transpose[bootFDM];
bootFDMmin=Sort[#][[50]]&/@Transpose[bootFDM];

SortedFDM=Sort[Transpose[{FDM,bootFDMmax,bootFDMmin,LabelPairs}]][[1;;,1]];
SortedFDMmin=Sort[Transpose[{FDM,bootFDMmax,bootFDMmin,LabelPairs}]][[1;;,3]];
SortedFDMmax=Sort[Transpose[{FDM,bootFDMmax,bootFDMmin,LabelPairs}]][[1;;,2]];
SortedLabelPairs=Sort[Transpose[{FDM,bootFDMmax,bootFDMmin,LabelPairs}]][[1;;,4]];

lowlines=Flatten[Position[bootFDMmax,#]&/@Select[bootFDMmax,#<1&]];
highlines=Flatten[Position[bootFDMmin,#]&/@Select[bootFDMmin,#>1&]];


Print[
Graphics[
{
{Gray,Thick,Line[{{0,1},{Length[SortedFDM],1}}]},
{RGBColor[0.138811, 0.31458, 0.424048],Line[Table[{x,SortedFDMmin[[x]]},{x,Length[SortedFDMmin]}]]},
{RGBColor[0.138811, 0.31458, 0.424048],Line[Table[{x,SortedFDMmax[[x]]},{x,Length[SortedFDMmax]}]]},
{PointSize[0.015],Point[#]}&/@Table[{x,SortedFDM[[x]]},{x,Length[SortedFDM]}],
{FontFamily->"Helvetica",Table[Text[SortedLabelPairs[[x]],{x,SortedFDMmin[[1]]-0.1},{-1,0},{0,1}],{x,Length[SortedLabelPairs]}]},
{FontFamily->"Helvetica",Text[ToString[PaddedForm[#,{3,2}]],{-2,#}]&/@{SortedFDMmin[[1]],SortedFDM[[1]],1.00,SortedFDM[[-1]],SortedFDMmax[[-1]]}}
}

,AspectRatio->1/GoldenRatio, PlotLabel->"FDM Plot (Relative Difference)",BaseStyle->{FontFamily->"Helvetica"}]
];


consensus2=Partition[Mean[Procrustes[Partition[Flatten[samp2],Length[samp2[[1]]]*2],Length[samp2[[1]]],2]],2];
points={Gray,PointSize[0.05],Point[#]}&/@consensus2;
lines=Flatten[Table[Table[Line[{consensus2[[i]],consensus2[[j]]}],{i,j+1,Length[consensus2],1}],{j,Length[consensus2]-1}]];

Print[
Graphics[{{RGBColor[0.138811, 0.31458, 0.424048],lines[[lowlines]]},{Thick,RGBColor[0.409247, 0.0267643, 0.0470741],lines[[highlines]]},points}]
];

]


(*  This function projects a phylogenetic tree into a GMM morphospace using PGLS
    to estimate the most likely ancestral shapes under a Brownian motion model of
    evolution, then projecting the node shapes into the principal components space
    of the shape data.

    Usage:     TreeToMorphospace[proc, labels, PCs, tree]
               where proc is a matrix of procrustes superimposed 2D landmark data, labels
               is a list of string labels for each row in the procrustes matrix, PCs is a
               list of digits specifying which principal components to plot, and tree
               is a string describing the tree and its branch lengths in Newick format.
               The labels and the tip names in the tree must match exactly.


    Created:  20 May 2012.
	Updated:  11 April 2014 to exclude PCs with less than 1% of the total variance from calculations.  
*)

TreeToMorphospace[ProcrustesCoords_,TaxonLabels_,PCs_,Tree_]:=Block[{Consensus, Resids, P, u, w, v, Scores, Tips, TipLabs,PltRnge,TreeTable,TipEntries,varY,varAY,varA,Rates,Nodes,SEs,NodePts,NodeEntries,NdLabs,Branches,Pts,TipPositions,NodeLabels,NodePositions,CEllipses,CIs},

TreeTable=Sort[TreeToTable[Tree][[2;;]]];
TipPositions=Flatten[Position[TreeTable[[1;;,1]],#]&/@TaxonLabels];
If[TreeTable[[TipPositions,1]]==TaxonLabels,"WunderBar",Return["Labels do not match taxon names in tree"]];
TipEntries=TreeTable[[TipPositions]];
NodeLabels=Sort[Union[TreeTable[[1;;,2]]]];
NodePositions=Flatten[Position[TreeTable[[1;;,1]],#]&/@NodeLabels];
NodeEntries=TreeTable[[NodePositions]];
TreeTable=TreeTable[[Flatten[{TipPositions,NodePositions}]]];
{varY,varAY,varA}=PhylogeneticMatrices[TreeTable];

Consensus=Mean[ProcrustesCoords];
Resids=#-Consensus&/@ProcrustesCoords;
P=Covariance[Resids];
{u,w,v}=SingularValueDecomposition[P];
w=Tr[w,List];
Scores=Resids . u;
Scores=Scores[[1;;,1;;Length[Select[w,#>=0.01*(Plus@@w)&]]]];


{Rates,Nodes,SEs}=Transpose[ReconstructNodes[TreeTable,Transpose[{TaxonLabels,#}]]&/@Transpose[Scores]];
Rates=#[[2]]&/@Rates;Nodes=Transpose[#[[1;;,2]]&/@Nodes];
SEs=Transpose[#[[1;;,2]]&/@SEs];

NodePts={RGBColor[0.621027, 0.491814, 0.34374],PointSize[0.015],Point[#]}&/@Nodes[[1;;,PCs]];
NdLabs={FontFamily->"Arial",Table[Text[NodeLabels[[x]],Nodes[[1;;,PCs]][[x]],{-1.3,0}],{x,Length[NodeLabels]}]};
Pts=Partition[Flatten[Append[Partition[Flatten[Transpose[{TaxonLabels,Scores[[1;;,PCs]]}]],3],
Transpose[{NodeLabels,Nodes[[1;;,PCs]]}]]],3];
Branches=Line[{Flatten[Cases[Pts,{#[[1]],___,___}]][[2;;3]],Flatten[Cases[Pts,{#[[2]],___,___}]][[2;;3]]}]&/@TreeTable;

CEllipses={
{Opacity[.1],RGBColor[0.728496, 0.427222, 0.162768],Table[Disk[Nodes[[x,PCs]],1.96*SEs[[x,PCs]]],{x,Length[Nodes]}]}
};

CIs={Opacity[.25],Dashed,RGBColor[0.728496, 0.427222, 0.162768],{
Table[Line[{{Nodes[[x,PCs[[1]]]]-SEs[[x,PCs[[1]]]],Nodes[[x,PCs[[2]]]]},{Nodes[[x,PCs[[1]]]]+SEs[[x,PCs[[1]]]],Nodes[[x,PCs[[2]]]]}}],{x,Length[Nodes]}],
Table[Line[{{Nodes[[y,PCs[[1]]]],Nodes[[y,PCs[[2]]]]-SEs[[y,PCs[[2]]]]},{Nodes[[y,PCs[[1]]]],Nodes[[y,PCs[[2]]]]+SEs[[y,PCs[[2]]]]}}],{y,Length[Nodes]}]
}};

Tips={PointSize[0.02],RGBColor[0.409247, 0.0267643, 0.0470741],Table[Point[Scores[[1;;,PCs]][[x]]],{x,Length[Scores]}]};
TipLabs={FontFamily->"Arial",Table[Text[TaxonLabels[[x]],Scores[[1;;,PCs]][[x]],{-1.5,0}],{x,Length[Scores]}]};PltRnge=1.3*{{Min[Scores[[1;;,PCs[[1]]]]],Max[Scores[[1;;,PCs[[1]]]]]},{Min[Scores[[1;;,PCs[[2]]]]],Max[Scores[[1;;,PCs[[2]]]]]}};

Return[Graphics[{CEllipses,CIs,Branches,Tips,TipLabs,NodePts,NdLabs},BaseStyle->{PointSize->10,FontFamily->"Arial"},Frame->True,PlotRange->All,FrameLabel->{"PC "<>ToString[PCs[[1]]],"PC "<>ToString[PCs[[2]]]}]];
];



(* This function does a multivariate least squares regression of shape onto a single predictor variable. The function calculates the PC scores
   of the Procrustes superimposed shape, then regresses them onto the input variable.  The function calculates an r-squared value for the
   regression and assesses its significance by randomization.  The scores are randomized with respect to the input variable 10,000 times
   and the observed r-squared value compared to the randomized distribution.  The function also reports the intercept, regression, and 
   univariate r-square for each of the non-zero PCs.  A graph showing the scores and regression of one PC is returned.  The first PC is graphed
   unless a different one is specified.

	Usage:  ShapeRegress[proc, X, 3]
	where proc is a set of Procrustes superimposed coordinates, X is a continuous variable vector with an observation for each
	object in proc, and "3" indicates that the regression of PC3 onto X should be plotted in the output.

	Created:  19 April 2012

*)

ShapeRegress[proc_,Var_,PC_:1]:=Block[{Y,X,lm,P,EHat,Rsquare,MinX,ln,MaxX,BootR2,BootX,Bootlm,BootP,BootYHat,BootEHat,Pvalue,
BootIter,UniVarRsquare,scores},
If[Length[Dimensions[Var]]>1, Return["Only one X variable is allowed.  Consider using TwoBlockPartialLeastSquares[] instead."]];
scores=PrincipalComponents[proc];
Y=scores[[1;;,1;;MatrixRank[scores]]];
X=Transpose[{Table[1,{Length[Var]}],Var}];
lm=Inverse[Transpose[X] . X] . Transpose[X] . Y;
P=X . Inverse[Transpose[X] . X] . Transpose[X];
EHat=Y-P . Y;
Rsquare=(Plus@@(Variance[Y]-Variance[EHat]))/(Plus@@Variance[Y]);
UniVarRsquare=NumberForm[#,{3,2}]&/@(((Variance[Y]-Variance[EHat]))/(Variance[Y]));
BootR2={};
BootIter=10000;
Do[
BootX=Transpose[{Table[1,{Length[Var]}],Var[[RandomSample[Range[Length[Var]]]]]}];
Bootlm=Inverse[Transpose[BootX] . BootX] . Transpose[BootX] . Y;
BootP=BootX . Inverse[Transpose[BootX] . BootX] . Transpose[BootX];
BootYHat=BootP . Y;
BootEHat=Y-BootP . Y;
BootR2=Append[BootR2,(Plus@@(Variance[Y]-Variance[BootEHat]))/(Plus@@Variance[Y])];
,{BootIter}];
Pvalue=Length[Select[BootR2,#>=Rsquare&]]/BootIter//N;

MinX=Min[X[[1;;,2]]];MaxX=Max[X[[1;;,2]]];
ln=Line[{{MinX,MinX*lm[[2,PC]]+lm[[1,PC]]},{MaxX,MaxX*lm[[2,PC]]+lm[[1,PC]]}}];
 
Return[
{
Graphics[{{PointSize[0.02],Point[#]}&/@Transpose[{X[[1;;,2]],Y[[1;;,PC]]}],ln},Frame->True, ImageSize->500,BaseStyle->{FontFamily->"Arial", FontSize->10},FrameLabel->{"Var","PC "<>ToString[PC]},AspectRatio->1/GoldenRatio],
"R-square (all PCs) = "<>ToString[NumberForm[Rsquare,{3,2}]],
"P[R-square is random] = "<>ToString[NumberForm[Pvalue,{3,2}]],
TableForm[Partition[Flatten[{lm,UniVarRsquare}],Length[Y[[1]]]],TableHeadings->{{"Intercept","Slope","Univariate R-square"},Table["PC"<>ToString[x],{x,Length[Y[[1]]]}]}]
}//TableForm
]; 
];


(* This function performs a two-block partial least squares analysis following the
   methodology of Rohlf and Corti (2000).  Two blocks of data are given to the function,
   along with a list of two strings indicating the type of data.  Allowable types are
   "Shape" (Procrustes superimposed coordinates), "Standardized" (independent variables 
   with different units of measurement that need to be standardized), and "Unstandardized"
   (independent variables with the same unit of measurement that do not need to be 
   standardized).  An optional argument is the number of the PLS axis to plot in the 
   output graph.  By default PLS 1 is plotted.

   Usage:   TwoBlockPartialLeastSquares[proc, X, {"Shape","Standardized"}]
   where proc is a table of Procrustes superimposed coordinates, X is a table
   of independent variables with labels in the first row.

   Created:  20 April 2012

*)

TwoBlockPartialLeastSquares[Data1_,Data2_,Mode_,PLS_:1]:=Block[{Y1,Y2,N1,N2,R12,F1,Dmat,F2,Z1,Z2,con1,con2,corrs,output,labels1, labels2},
Which[
Mode=={"Shape","Shape"},con1=Mean[Data1];Y1=#-con1&/@Data1;con2=Mean[Data2];Y2=#-con2&/@Data2;,
Mode=={"Shape","Standardized"},con1=Mean[Data1];Y1=#-con1&/@Data1;Y2=Standardize[Data2[[2;;]]];labels2=Data2[[1]],
Mode=={"Shape","Unstandardized"},con1=Mean[Data1];Y1=#-con1&/@Data1;Y2=Data2[[2;;]]; labels2=Data2[[1]],
Mode=={"Standardized","Standardized"},Y1=Standardize[Data1[[2;;]]]; labels1=Data1[[1]]; Y2=Standardize[Data2[[2;;]]];labels2=Data2[[1]],
Mode=={"Unstandardized","Standardized"},Y1=Data1[[2;;]]; labels1=Data1[[1]];Y2=Standardize[Data2[[2;;]]]; labels2=Data2[[1]],
Mode=={"Unstandardized","Unstandardized"},Y1=Data1[[2;;]]; labels1=Data1[[1]];Y2=Data2[[2;;]]; labels2=Data2[[1]],
Mode==Mode,Return["Mode must be {\"Shape\",\"Shape\"},{\"Shape\",\"Standardized\"}, {\"Shape\",\"Unstandardized\"}, {\"Standardized\",\"Standardized\"}, {\"UnStandardized\",\"Standardized\"}, or  {\"Unstandardized\",\"Unstandardized\"}."];
];


N1=Plus@@Variance[Y1];
N2=Plus@@Variance[Y2];
R12=Covariance[Y1,Y2];
{F1,Dmat,F2}=SingularValueDecomposition[R12];
Dmat=Tr[Dmat,List];
Z1=Y1 . F1;
Z2=Y2 . F2;
corrs=Table[Correlation[Z1[[1;;,i]],Z2[[1;;,i]]],{i,Length[Dmat]}];

Which[
Mode=={"Shape","Shape"},output=GraphicsGrid[{{tpSpline[con1,con1+F1[[1;;,PLS]]*Max[Z1[[1;;,PLS]]]],tpSpline[con2,con2+F2[[1;;,PLS]]*Max[Z2[[1;;,PLS]]]]}},ImageSize->500,BaseStyle->{FontFamily->"Arial",FontSize->10},PlotLabel->"Positive End of PLS "<>ToString[PLS]],
Mode=={"Shape","Standardized"},
output={tpSpline[con1,con1+F1[[1;;,PLS]]*Max[Z1[[1;;,PLS]]]],
TableForm[Append[NumberForm[#,{4,4}]&/@F2[[1;;,PLS]],NumberForm[corrs[[PLS]],{3,2}]],TableHeadings->{Append[labels2,"Corr"],{"PLS"}}]},
Mode=={"Shape","Unstandardized"},
output={tpSpline[con1,con1+F1[[1;;,PLS]]*Max[Z1[[1;;,PLS]]]],
TableForm[Append[NumberForm[#,{4,4}]&/@F2[[1;;,PLS]],NumberForm[corrs[[PLS]],{3,2}]],TableHeadings->{Append[labels2,"Corr"],{"PLS"}}]},
Mode=={"Standardized","Standardized"},
output={TableForm[Append[NumberForm[#,{4,4}]&/@F1[[1;;,PLS]],NumberForm[corrs[[PLS]],{3,2}]],TableHeadings->{Append[labels1,"Corr"],{"PLS"}}],
TableForm[Append[NumberForm[#,{4,4}]&/@F2[[1;;,PLS]],NumberForm[corrs[[PLS]],{3,2}]],TableHeadings->{Append[labels2,"Corr"],{"PLS"}}]},
Mode=={"Unstandardized","Standardized"},
output={TableForm[Append[NumberForm[#,{4,4}]&/@F1[[1;;,PLS]],NumberForm[corrs[[PLS]],{3,2}]],TableHeadings->{Append[labels1,"Corr"],{"PLS"}}],
TableForm[Append[NumberForm[#,{4,4}]&/@F2[[1;;,PLS]],NumberForm[corrs[[PLS]],{3,2}]],TableHeadings->{Append[labels2,"Corr"],{"PLS"}}]},
Mode=={"Unstandardized","Unstandardized"},
output={TableForm[Append[NumberForm[#,{4,4}]&/@F1[[1;;,PLS]],NumberForm[corrs[[PLS]],{3,2}]],TableHeadings->{Append[labels1,"Corr"],{"PLS"}}],
TableForm[Append[NumberForm[#,{4,4}]&/@F2[[1;;,PLS]],NumberForm[corrs[[PLS]],{3,2}]],TableHeadings->{Append[labels2,"Corr"],{"PLS"}}]}
];


Return[{
"Proportion of Actual Squared Covariance Explained by PLS "<>ToString[PLS]<>":  "<>ToString[NumberForm[((Dmat^2)/(Plus@@(Dmat^2)))[[PLS]],{3,2}]],
"Proportion of Total Possible Squared Covariance Explained by All PLS Axes:  "<>ToString[NumberForm[(Plus@@(Dmat^2))/(N1*N2),{3,2}]],
Graphics[{PointSize[0.02],Point[#]}&/@Transpose[{Z1[[1;;,PLS]],Z2[[1;;,PLS]]}],AspectRatio->1/GoldenRatio,Axes->False,Frame->True,BaseStyle->{FontFamily->"Arial",FontSize->10},ImageSize->400,FrameLabel->{"Data 1, PLS "<>ToString[PLS],"Data 2, PLS "<>ToString[PLS]}],
output
}//TableForm];

];



(* 
   This function projects a phylogenetic tree and fossils into a GMM morphospace. 
   The function uses PGLS on the tip taxa in the tree to estimate the most likely 
   ancestral shapes under a Brownian motion model of evolution, then projecting the 
   node shapes into the principal components space of the shape data. The function also 
   projects the remaining non-tip shapes, which are presumed to be fossil candidates
   for ancestors, into the morphospace.

    Usage:     TreeAndFossilsToMorphospace[proc, labels, PCs, tree]
               where proc is a matrix of procrustes superimposed 2D landmark data, labels
               is a list of string labels for each row in the procrustes matrix, PCs is a
               list of digits specifying which principal components to plot, and tree
               is a string describing the tree and its branch lengths in Newick format.
               All tip names in the tree string must be in the list of labels.  All labels
               that are not in the tree string are presumed to be fossils.


    Created:  20 May 2012.
	Upadated: 11 April 2014 to exclude PCs with less than 1% of the total variance from the calculations.

*)

TreeAndFossilsToMorphospace[ProcrustesCoords_,TaxonLabels_,PCs_,Tree_]:=Block[{Consensus, Resids, P, u, w, v, Scores, 
Tips, TipLabs,PltRnge,TreeTable,TipEntries,Rates,Nodes,SEs,NodePts,NodeEntries,
NdLabs,Branches,Pts,TipPositions,NodeLabels,NodePositions,CEllipses,CIs,tiplabels,
TreePositions,FossilLabels,FossilPositions,TipScores,FossilResids,FossilScores,Fossils,FossilLabs},

TreeTable=Sort[TreeToTable[Tree][[2;;]]];

tiplabels=TreeTable[[Flatten[Position[TreeTable[[1;;,4]],1]],1]];
TreePositions=Flatten[Position[TaxonLabels,#]&/@tiplabels];
FossilLabels=DeleteCases[Flatten[If[MemberQ[tiplabels,#],Null,#]&/@TaxonLabels],Null];
FossilPositions=Flatten[Position[TaxonLabels,#]&/@FossilLabels];

TipPositions=Flatten[Position[TreeTable[[1;;,1]],#]&/@TaxonLabels[[TreePositions]]];
If[Length[TreePositions]==Length[tiplabels],"WunderBar",Return["Not all tree taxa are present in the Procrustes coordinate matrix."]];
If[Length[TreePositions]==Length[TaxonLabels],Return["There are no non-tip shapes in the data set."]];
TipEntries=TreeTable[[TipPositions]];
NodeLabels=Sort[Union[TreeTable[[1;;,2]]]];
NodePositions=Flatten[Position[TreeTable[[1;;,1]],#]&/@NodeLabels];
NodeEntries=TreeTable[[NodePositions]];
TreeTable=TreeTable[[Flatten[{TipPositions,NodePositions}]]];

Consensus=Mean[ProcrustesCoords[[TreePositions]]];
Resids=#-Consensus&/@ProcrustesCoords[[TreePositions]];
P=Covariance[Resids];
{u,w,v}=SingularValueDecomposition[P];
w=Tr[w,List];
TipScores=Resids . u;
TipScores=TipScores[[1;;,1;;Length[Select[w,#>=0.01*(Plus@@w)&]]]];

FossilResids=#-Consensus&/@ProcrustesCoords[[FossilPositions]];
FossilScores=FossilResids . u;
FossilScores=FossilScores[[1;;,1;;Length[Select[w,#>=0.01*(Plus@@w)&]]]];


Scores=Partition[Flatten[Append[TipScores[[1;;,PCs]],FossilScores[[1;;,PCs]]]],2];

{Rates,Nodes,SEs}=Transpose[ReconstructNodes[TreeTable,Transpose[{TaxonLabels[[TreePositions]],#}]]&/@Transpose[TipScores]];
Rates=#[[2]]&/@Rates;Nodes=Transpose[#[[1;;,2]]&/@Nodes];
SEs=Transpose[#[[1;;,2]]&/@SEs];

NodePts={RGBColor[0.621027, 0.491814, 0.34374],PointSize[0.015],Point[#]}&/@Nodes[[1;;,PCs]];
NdLabs={FontFamily->"Arial",Table[Text[NodeLabels[[x]],Nodes[[1;;,PCs]][[x]],{-1.3,0}],{x,Length[NodeLabels]}]};
Pts=Partition[Flatten[Append[Partition[Flatten[Transpose[{TaxonLabels[[TreePositions]],TipScores[[1;;,PCs]]}]],3],
Transpose[{NodeLabels,Nodes[[1;;,PCs]]}]]],3];
Branches=Line[{Flatten[Cases[Pts,{#[[1]],___,___}]][[2;;3]],Flatten[Cases[Pts,{#[[2]],___,___}]][[2;;3]]}]&/@TreeTable;

CEllipses={Opacity[.1],
RGBColor[0.728496, 0.427222, 0.162768],
Table[
Disk[EllipsoidQuantile[MultinormalDistribution[Nodes[[x,PCs]],DiagonalMatrix[SEs[[x,PCs]]^2]],0.95][[1]],
EllipsoidQuantile[MultinormalDistribution[Nodes[[x,PCs]],DiagonalMatrix[SEs[[x,PCs]]^2]],0.95][[2]]]
,{x,Length[Nodes]}]
};

CIs={Opacity[.25],Dashed,RGBColor[0.728496, 0.427222, 0.162768],
{
Table[Line[{{Nodes[[x,PCs[[1]]]]-SEs[[x,PCs[[1]]]],Nodes[[x,PCs[[2]]]]},{Nodes[[x,PCs[[1]]]]+SEs[[x,PCs[[1]]]],Nodes[[x,PCs[[2]]]]}}],{x,Length[Nodes]}],
Table[Line[{{Nodes[[y,PCs[[1]]]],Nodes[[y,PCs[[2]]]]-SEs[[y,PCs[[2]]]]},{Nodes[[y,PCs[[1]]]],Nodes[[y,PCs[[2]]]]+SEs[[y,PCs[[2]]]]}}],{y,Length[Nodes]}]
}
};

Tips={PointSize[0.02],RGBColor[0.409247, 0.0267643, 0.0470741],Table[Point[TipScores[[1;;,PCs]][[x]]],{x,Length[TipScores]}]};
TipLabs={FontFamily->"Arial",Table[Text[TaxonLabels[[TreePositions]][[x]],TipScores[[1;;,PCs]][[x]],{-1.5,0}],{x,Length[TipScores]}]};

Fossils={PointSize[0.02],RGBColor[0.426566, 0.474784, 0.257816],Table[Point[FossilScores[[1;;,PCs]][[x]]],{x,Length[FossilScores]}]};
FossilLabs={FontFamily->"Arial",Table[Text[TaxonLabels[[FossilPositions]][[x]],FossilScores[[1;;,PCs]][[x]],{-1.5,0}],{x,Length[FossilScores]}]};

PltRnge=2*{{Min[TipScores[[1;;,PCs[[1]]]]],Max[TipScores[[1;;,PCs[[1]]]]]},{Min[TipScores[[1;;,PCs[[2]]]]],Max[TipScores[[1;;,PCs[[2]]]]]}};

Return[Graphics[{CEllipses,CIs,Branches,Tips,TipLabs,NodePts,NdLabs,Fossils,FossilLabs},BaseStyle->{PointSize->10,FontFamily->"Helvetica"},Frame->True,PlotRange->All,FrameLabel->{"PC "<>ToString[PCs[[1]]],"PC "<>ToString[PCs[[2]]]}]];

];





(* 
   This function calculates the probability that a shape, usually a fossil, falls within
   the expected distribution of an ancestral node under the assumption of a Brownian motion 
   model of evolution.  The node values and their covariances are calculated from the tree and 
   the taxa whose names appear in its tip labels using the ReconstructAncestors[] function from 
   the Phylogenetics package. P-values are then calculated for the remaining taxa in the data
   set.  P-values are multivariate, based on all dimensions of the morphospace, unless the optional
   value PCs is smaller than the total number of morphospace dimensions.  For each fossils, the 
   squared Mahalanobis distance between it and the node is calculated using the scores of the fossil, 
   the scores of the node, and the covariance matrix of the node (a diagonal matrix of the squared
   standard deviations associated with the evolutionary process).  P is taken from the Chi-Squared
   distribution whose degrees of freedom equals the number of morphospace dimensions being considered. 
   The P-value is interpreted as the probability that the fossil could represent the ancestral population
   at that node given a Brownian motion process of evolution, the rate of evolution estimated from the
   tree, the topology of the tree, and the shapes of the tip taxa.  Remember that the ancestral 
   estimates are population means rather than individuals, so the interpretation is based on the fossil
   being representative of the mean of the population from which it came.  Ideally, each tip and
   fossil taxon will be represented in this analysis by a sample mean shape (consensus shape) of a larger
   sample.  P-values greater than 0.05 are highlighted in the results as being significantly compatible
   with the distribution of a node.  Remember that a fossil that did not live at the time of the 
   ancestor cannot logically be an ancestor, regardless of how similar its morphology is to the 
   expected ancestral shape.

   This function is a companion to the TreeAndFossilsToMorphospace[] function in this package 
   and to the ReconstructNodes[] function in the Phylogenetics package.


    Usage:     ProbabilitiesOfShapesAsAncestors[proc , labels, tree]
               where proc is a matrix of procrustes superimposed 2D landmark data, labels
               is a list of string labels for each row in the procrustes matrix, and tree
               is a string describing the tree and its branch lengths in Newick format.
               All tip names in the tree string must be in the list of labels.  All labels
               that are not in the tree string are presumed to be fossils.


    Created:  20 May 2012.

*)

ProbabilitiesOfShapesAsAncestors[ProcrustesCoords_,TaxonLabels_,Tree_,PCs_:5000]:=Block[{TreeTable,tiplabels,TreePositions,FossilLabels,FossilPositions,TipPositions,TipEntries,NodeLabels,NodePositions,NodeEntries,varY,varAY,varA,Consensus,Resids,P,Rank,u,w,v,TipScores,FossilResids,FossilScores,Rates,Nodes,SEs,m,p,x,i,j},

If[StringQ[Tree], TreeTable=TreeToTable[Tree][[2;;]],If[Tree[[1]]=={"Descendant","Ancestor","Branch Length","Tip?"},
TreeTable=Tree[[2;;]],
If[Length[Tree[[1]]]==4,TreeTable=Tree,Print["Tree variable must be a Newick tree string or a tree table."];];];
];

If[PCs<=0,Return["Number of PCs must be greater than 0."]];


tiplabels=TreeTable[[Flatten[Position[TreeTable[[1;;,4]],1]],1]];
TreePositions=Flatten[Position[TaxonLabels,#]&/@tiplabels];
FossilLabels=DeleteCases[Flatten[If[MemberQ[tiplabels,#],Null,#]&/@TaxonLabels],Null];
FossilPositions=Flatten[Position[TaxonLabels,#]&/@FossilLabels];



TipPositions=Flatten[Position[TreeTable[[1;;,1]],#]&/@TaxonLabels[[TreePositions]]];
If[Length[TreePositions]==Length[tiplabels],"WunderBar",Return["Not all tree taxa are present in the Procrustes coordinate matrix."]];
TipEntries=TreeTable[[TipPositions]];
NodeLabels=Sort[Union[TreeTable[[1;;,2]]]];
NodePositions=Flatten[Position[TreeTable[[1;;,1]],#]&/@NodeLabels];
NodeEntries=TreeTable[[NodePositions]];
TreeTable=TreeTable[[Flatten[{TipPositions,NodePositions}]]];
{varY,varAY,varA}=PhylogeneticMatrices[TreeTable];

Consensus=Mean[ProcrustesCoords[[TreePositions]]];
Resids=#-Consensus&/@ProcrustesCoords[[TreePositions]];
P=Covariance[Resids];
Rank=MatrixRank[P];
If[PCs<=Rank,Rank=PCs,Rank=Rank];
{u,w,v}=SingularValueDecomposition[P];
w=w[[1;;Rank,1;;Rank]];
TipScores=Resids . u[[1;;,1;;Rank]];

FossilResids=#-Consensus&/@ProcrustesCoords[[FossilPositions]];
FossilScores=FossilResids . u[[1;;,1;;Rank]];

{Rates,Nodes,SEs}=Transpose[ReconstructNodes[TreeTable,Transpose[{TaxonLabels[[TreePositions]],#}]]&/@Transpose[TipScores]];
Rates=#[[2]]&/@Rates;Nodes=Transpose[#[[1;;,2]]&/@Nodes];
SEs=Transpose[#[[1;;,2]]&/@SEs];

Return[Text[Style[TableForm[Table[Table[
m=(FossilScores[[i]]-Nodes[[j]]) . Inverse[DiagonalMatrix[SEs[[j]]^2]] . (FossilScores[[i]]-Nodes[[j]]);
AccountingForm[NumberForm[Style[p=(1-CDF[ChiSquareDistribution[Rank],m]),If[p>=0.05,Bold,GrayLevel[0.7]]],{3,4}],4,NumberPadding->{" ",""}],{j,Length[Nodes]}],{i,Length[FossilScores]}],TableHeadings->{FossilLabels,NodeLabels}],FontFamily->"Arial"]]];

];




(* 

This function creates thin-plate spline representations of the shapes of taxa at the tips and nodes of a phylognetic tree.  The node shapes are
reconstructed using the PGLS method assuming a Brownian motion mode of evolution (see ReconstructNodesSimple[] function in Phylogenetics package).  The
method essentially follows Rohlf (2002).   The thin-plate spline representations are shown as the deformation of each tip and node from the shape 
at the root of the tree, thus showing derived changes.  Note that the shape reconstructions themselves are done from a shape space centered at 
the mean (consensus) of the tip taxa as recommended by Rohlf (1998).

   This function is a companion to the TreeToMorphospace[] function in this package 
   and to the ReconstructNodes[] function in the Phylogenetics package.


    Usage:     ReconstructAncestorShapes[proc , labels, tree]
               where proc is a matrix of procrustes superimposed 2D landmark data, labels
               is a list of string labels for each row in the procrustes matrix, and tree
               is a string describing the tree and its branch lengths in Newick format.
               All tip names in the tree string must be in the list of labels.  All labels
               that are not in the tree string are presumed to be fossils.


    Created:  20 May 2012.
	Updated:  11 April 2014
	Updated:   5 Mar 2017 to fix bug in the number of PCs used.  Now the number of PCs equals the rank of the covariance matrix. 
			   Also added simplified function for ReconstructNodes to speed up calculation.


*)

ReconstructAncestorShapes[proc_,labels_,tree_]:=Block[{consensus,resids,P,rank,u,w,v,scores,trait,x,rates,nodes,SDs,nodelabels,nodeshapes,rootposition},
consensus=Mean[proc];
resids=#-consensus&/@proc;
P=Covariance[resids];
rank=MatrixRank[P];
{u,w,v}=SingularValueDecomposition[P];
u=u[[1;;,1;;rank]];
w=Tr[w,List][[1;;rank]];
scores=resids . u;
nodes=ReconstructNodesSimple[tree,Transpose[{labels,#}]]&/@Transpose[scores];
nodelabels=nodes[[1,1;;,1]];
nodes=Transpose[#[[1;;,2]]&/@nodes];
rootposition=Flatten[Position[nodelabels, "Node 0"]];
nodeshapes=Partition[#+consensus,2]&/@(nodes . Transpose[u]);
Return[GraphicsGrid[{Table[tpSpline[nodeshapes[[rootposition]],proc[[x]],labels[[x]]],{x,Length[proc]}],Table[tpSpline[nodeshapes[[rootposition]],nodeshapes[[x]],nodelabels[[x]]],{x,Length[nodeshapes]}]}]];
];




(* 

	This function returns the estimated ancestral shapes for geometric morphometric
    landmark data that have been Procrustes superimposed.  It uses the Browninan motion 
    ancestral node reconstruction method described by Martins and Hansen (1997).  

   This function is a companion to the TreeToMorphospace[] function in this package 
   and to the ReconstructNodesSimple[] function in the Phylogenetics package, and it requires the latter package to 
	be installed.


    Usage:     AncestorShapes[proc, labels, tree]
               where proc is a matrix of procrustes superimposed 2D landmark data, labels
               is a list of string labels for each row in the procrustes matrix, PCs is a
               list of digits specifying which principal components to plot, and tree
               is a string describing the tree and its branch lengths in Newick format.
               All tip names in the tree string must be in the list of labels.  All labels
               that are not in the tree string are presumed to be fossils.


    Created:  23 July 2012.
	Updated:  11 April 2014 to stop calculation for PCs whose variance is less than 1% of the total. 
	Updated:   5 Mar 2017 to fix bug in the number of PCs used.  Now the number of PCs equals the rank of the covariance matrix. 
			   Also added simplified function for ReconstructNodes to speed up calculation.


*)

AncestorShapes[proc_,labels_,tree_]:=Block[{consensus,resids,P,rank,u,w,v,scores,trait,x,rates,nodes,SDs,nodelabels,nodeshapes,rootposition},
consensus=Mean[proc];
resids=#-consensus&/@proc;
P=Covariance[resids];
rank=MatrixRank[P];
{u,w,v}=SingularValueDecomposition[P];
u=u[[1;;,1;;rank]];
w=Tr[w,List][[1;;rank]];
scores=resids . u;
trait=Table[Prepend[scores[[x]],labels[[x]]],{x,Length[scores]}];
nodes=ReconstructNodesSimple[tree,Transpose[{labels,#}]]&/@Transpose[scores];
nodelabels=nodes[[1,1;;,1]];
nodes=Transpose[#[[1;;,2]]&/@nodes];
nodeshapes=#+consensus&/@(nodes . Transpose[u]);
nodeshapes=Table[Prepend[nodeshapes[[x]],nodelabels[[x]]],{x,Length[nodeshapes]}];
Return[nodeshapes];
];



(* 

	This function does a phylogenetic principal components of shape using the techniques described by Revell (2009).  
	The function also projects the phylogenetic tree into the pPCA space using the techniques descrbed by Rohlf 
	(2002), making use of the Browninan motion ancestral node reconstruction method 
	described by Martins and Hansen (1997).  

   This function is a companion to the TreeToMorphospace[] function in this package 
   and to the ReconstructNodes[] function in the Phylogenetics package, and it requires the latter package to 
	be installed.


    Usage:     PhylogenticPrincipalComponentsOfShape[proc, labels, PCs, tree]
               where proc is a matrix of procrustes superimposed 2D landmark data, labels
               is a list of string labels for each row in the procrustes matrix, PCs is a
               list of digits specifying which principal components to plot, and tree
               is a string describing the tree and its branch lengths in Newick format.
               All tip names in the tree string must be in the list of labels.  All labels
               that are not in the tree string are presumed to be fossils.


    Created:  23 July 2012.


*)

PhylogeneticPrincipalComponentsOfShape[proc_,labels_,PCs_,tree_]:=Block[{ancs,VarY,VarAY,VarA,phyloresids,pP,u,v,w,Scores,Nodes,NodeLabels,TaxonLabels,TreeTable,NodePts,NdLabs,x,Pts,Branches,Tips,TipLabs,PltRnge,rootposition},
{VarY,VarAY,VarA}=PhylogeneticMatrices[tree];
VarY=VarY[[2;;,2;;]];
ancs=AncestorShapes[proc,labels,tree];
rootposition=Flatten[Position[ancs[[1;;,1]], "Node 0"]];
phyloresids=#-Flatten[ancs[[rootposition,2;;]]]&/@proc;
pP=(Transpose[phyloresids] . Inverse[VarY] . phyloresids)/(Length[proc]-1);
{u,v,w}=SingularValueDecomposition[pP];
Scores=Chop[phyloresids . u];
Nodes=Chop[(#-Flatten[ancs[[rootposition,2;;]]]&/@ancs[[1;;,2;;]]) . u];
NodeLabels=ancs[[1;;,1]];
TaxonLabels=labels;
TreeTable=TreeToTable[tree][[2;;]];


NodePts={RGBColor[0.621027, 0.491814, 0.34374],PointSize[0.015],Point[#]}&/@Nodes[[1;;,PCs]];
NdLabs={FontFamily->"Arial",Table[Text[NodeLabels[[x]],Nodes[[1;;,PCs]][[x]],{-1.3,0}],{x,Length[NodeLabels]}]};
Pts=Partition[Flatten[Append[Partition[Flatten[Transpose[{TaxonLabels,Scores[[1;;,PCs]]}]],3],Transpose[{NodeLabels,Nodes[[1;;,PCs]]}]]],3];
Branches=Line[{Flatten[Cases[Pts,{#[[1]],___,___}]][[2;;3]],Flatten[Cases[Pts,{#[[2]],___,___}]][[2;;3]]}]&/@TreeTable;
Tips={PointSize[0.02],RGBColor[0.409247, 0.0267643, 0.0470741],Table[Point[Scores[[1;;,PCs]][[x]]],{x,Length[Scores]}]};TipLabs={FontFamily->"Arial",Table[Text[TaxonLabels[[x]],Scores[[1;;,PCs]][[x]],{-1.5,0}],{x,Length[Scores]}]};
PltRnge=1.3*{{Min[Scores[[1;;,PCs[[1]]]]],Max[Scores[[1;;,PCs[[1]]]]]},{Min[Scores[[1;;,PCs[[2]]]]],Max[Scores[[1;;,PCs[[2]]]]]}};

Return[Graphics[{Branches,Tips,TipLabs,NodePts,NdLabs},BaseStyle->{PointSize->10,FontFamily->"Arial"},Frame->True,PlotRange->All,FrameLabel->{"pPC "<>ToString[PCs[[1]]],"pPC "<>ToString[PCs[[2]]]}]];
];




(* 

This function replaces missing landmarks by flipping their bilateral counterpart. It 
works on 2D data from bilaterally symmetric objects where there are at least two
landmarks on the midline.  Missing landmarks are coded with {"NA", "NA"}.  

Usage:     ReplaceMissingWithBilateral[object, midline, bilaterals]

           where object is a table of landmark coordinates with landmarsk in the
           rows and x and y coordinates in the columns.  The x and y columns of
           missing landmarks are filled with the string "NA".  Midline is a list
           containing the number of two midline landmarks.  For example, if the 
           first and fifth landmarks are on the midline, then midline is {1,5}.
           Bilaterals is a list containing the numbers of pairs of bilaterally 
           symmetric landmarks.  For example, if landmarks 3&4 and 2&6 are 
           bilateral pairs, then this parameter is {{3,4},{2,6}}.
*)

ReplaceMissingWithBilateral[Object_,midline_,bilaterals_]:=Block[{MissingLandmarks,NewObject,MissingPosition,CorrespondingPosition,CenteredShape,RotationAngle,RotatedShape},

MissingLandmarks=Flatten[Position[Object,{"NA","NA"}]];
NewObject=Object;
NewObject[[MissingLandmarks]]={0,0};
MissingPosition=
Flatten[Position[bilaterals,#]&/@MissingLandmarks,1];
CorrespondingPosition=MissingPosition;

Do[If[CorrespondingPosition[[x,2]]==1,CorrespondingPosition[[x,2]]=2,CorrespondingPosition[[x,2]]=1],{x,Length[CorrespondingPosition]}];

CenteredShape=#-Object[[midline[[2]]]]&/@NewObject//N;
RotationAngle=VectorAngle[CenteredShape[[midline[[1]]]],{0,1}];

RotatedShape=Transpose[RotationMatrix[RotationAngle] . Transpose[CenteredShape]];
Do[
RotatedShape[[MissingLandmarks[[x]]]]={-1,1}*RotatedShape[[bilaterals[[CorrespondingPosition[[x,1]],CorrespondingPosition[[x,2]]]]]],{x,Length[MissingLandmarks]}];
CenteredShape=Transpose[RotationMatrix[-RotationAngle] . Transpose[RotatedShape]];
NewObject=#+Object[[midline[[2]]]]&/@CenteredShape//N;
Return[NewObject];
]





<<HypothesisTesting`

(* This function performs a one-way multivariate analysis of variance (MANOVA) of Procrustes superimposed shapes onto a grouping factor.  The 
   function calculates the PC scores of the Procrustes superimposed shape, then partitions sums of squares into within and between components. 
   Significance is assessed by non-parametric randomization, in which the real between-group sums of squares is compared to a distribution of between-group 
	sums of squares calculated by randomizing shapes among the grouping variables.  Post-hoc pairwise tests are also performed using randomization.  By default, 
	10,000 randomizations are used to produce this distribution, but the user can optionally change this with the third input parameter. The function 
	also calculates an F-ratio and uses the F-distribution to report a parametric p-value.  Note however, that the parametric p-value depends on the 
	assumption of normality, equal variance, and equal sampling in the groups, which is seldom the case for geometric morphometric shape data.   
	R-squared reports the proportion of the total variance explained by the grouping factor.

	Usage:  ShapeMANOVA[proc, labels, PC, iterations]
	where proc is a set of Procrustes superimposed coordinates, labels is a list of group label, PC indicates which principal component should be graphed
	(the choice does not affect the dimensions that are used in the test), and iterations is the number of randomizations performed in the non-parametric 
	test of significance (by default this is 10,000).

	Created:  28 February 2014
	Updated:   1 March 2014
	Updated:  11 April 2014 to prevent data with group N < 3 from being analyzed.
	Updated:  12 February 2022 to fix misreporting of total sum of squares in ANOVA table and F-ratio in MANOVA results table

*)

ShapeMANOVA[Proc_,GroupLabels_,PC_:1,Iterations_:10000]:=Block[{Consensus,Resids,x,CM,u,w,v,Scores,Groups,EachGroup,GroupScores,GroupMeans,GrandMean,TotalSS,RealBSS,RealWSS,RandomBSS,RandomizedScores,RandomizedMeans,P,Bdf,Wdf,Tdf,ANOVAtable,i,j,RandomResults,PairwiseBSS,RandomPairwiseBSS,PairwiseLabels,PairwiseP,MANOVAResults,PairwiseResults,ANOVAGraph,GroupNs},
Consensus = Mean[Proc];
Resids =Table[Proc[[x]]-Consensus, {x, Length[Proc]}];
CM = Covariance[Resids];
{u,w,v}=SingularValueDecomposition[CM];
Scores=Resids . u;
Scores=Scores[[1;;,1;;(MatrixRank[CM]-1)]];
Groups=Union[GroupLabels];
Bdf=Length[Groups]-1;
Wdf=Length[Scores]-Length[Groups];
Tdf=Length[Scores]-1;
EachGroup=Flatten[Position[GroupLabels,#]]&/@Groups;
GroupNs=Length[#]&/@EachGroup;
If[Min[GroupNs]<=2,Return["One or more groups have fewer than 3 individuals. A meaningful MANOVA cannot be performed."]];
GroupScores=Scores[[#]]&/@EachGroup;
GroupMeans={Length[#],Mean[#]}&/@GroupScores;
GrandMean=Mean[Scores];
RealBSS=Plus@@Plus@@((#[[1]]*((GrandMean-#[[2]])^2))&/@GroupMeans);
TotalSS=Plus@@Plus@@(((GrandMean-#)^2)&/@Scores);
PairwiseBSS=Flatten[Table[Table[Plus@@Plus@@((GroupMeans[[i,2]]-GroupMeans[[j,2]])^2),{j,i+1,Length[Groups]}],{i,Length[Groups]-1}]];
PairwiseLabels=Flatten[Table[Table[StringJoin[Groups[[i]],"-",Groups[[j]]],{j,i+1,Length[Groups]}],{i,Length[Groups]-1}]];
RealWSS=Plus@@Table[Plus@@(Plus@@(((#-GroupMeans[[x,2]])^2)&/@GroupScores[[x]])),{x,Length[GroupScores]}];
RandomResults=Table[RandomizedScores=RandomSample[Scores]; RandomizedMeans={Length[#],Mean[RandomizedScores[[#]]]}&/@EachGroup;{Plus@@Plus@@((#[[1]]*((GrandMean-#[[2]])^2))&/@RandomizedMeans),Flatten[Table[Table[Plus@@((RandomizedMeans[[i,2]]-RandomizedMeans[[j,2]])^2),{j,i+1,Length[Groups]}],{i,Length[Groups]-1}]]},{Iterations}];
RandomBSS=RandomResults[[1;;,1]];
RandomPairwiseBSS=RandomResults[[1;;,2]];
P=NumberForm[Length[Select[RandomBSS,#>=RealBSS&]]/10000//N,{3,4}];
PairwiseP=Transpose[{PairwiseLabels,Table[NumberForm[Length[Select[RandomPairwiseBSS[[1;;,x]],#>=PairwiseBSS[[x]]&]]/10000//N,{3,4}],{x,Length[PairwiseBSS]}]}];
ANOVAtable=TableForm[{{NumberForm[RealBSS,{3,4}],Bdf,NumberForm[RealBSS/Bdf,{3,4}]},{NumberForm[RealWSS,{3,4}],Wdf,NumberForm[RealWSS/Wdf,{3,4}]},{NumberForm[RealBSS+RealWSS,{3,4}], Length[Scores]-1,NumberForm[(RealBSS+RealWSS)/( Length[Scores]-1),{3,4}]}},TableHeadings->{{"Between", "Within","Total"},{"SS","df","MS"}}];
MANOVAResults=TableForm[{{"P (non-parametric)",P},
{"F-ratio",NumberForm[(RealBSS/Bdf)/(RealWSS/Wdf),{4,2}]},{"P (parametric)",NumberForm[Chop[FRatioPValue[(RealBSS/Bdf)/(RealWSS/Wdf),Bdf,Wdf][[2]]],{3,4}]},{"R-squared",NumberForm[1-(RealWSS/(RealBSS+RealWSS)),{3,4}]}},TableHeadings->{None,{"MANOVA results",None}}];
PairwiseResults=TableForm[PairwiseP,TableHeadings->{None,{"Post-hoc tests","P (non-parametric)"}}];
ANOVAGraph=Graphics[Table[{Style[Text[Groups[[x]],{x,Min[GroupScores[[1;;,1;;,PC]]]*1.25}],Medium],{RGBColor[0.629267, 0.637751, 0.646723],Thick,Line[{{x,Mean[GroupScores[[x,1;;,PC]]]+StandardDeviation[GroupScores[[x,1;;,PC]]]},{x,Mean[GroupScores[[x,1;;,PC]]]-StandardDeviation[GroupScores[[x,1;;,PC]]]}}]&/@GroupScores[[x]]},{PointSize[0.01],Point[{x,#}]&/@GroupScores[[x,1;;,PC]]},{RGBColor[0.428687, 0, 0.0558022],PointSize[0.02],Point[{x,Mean[GroupScores[[x,1;;,PC]]]}]}},{x,Length[GroupScores]}],AspectRatio->1,PlotRange->{{0,Length[GroupScores]+1},{Max[GroupScores[[1;;,1;;,PC]]]*1.25,Min[GroupScores[[1;;,1;;,PC]]]*1.5}},Axes->{False,True},AxesLabel->{None,"PC"<>ToString[PC]},BaseStyle->{FontFamily->"Helvetica"}];
Return[{ANOVAtable,MANOVAResults,PairwiseResults,ANOVAGraph}];
]



(* This function breaks an outline into segments based on landmark positions and places
   a specified number of equally spaced points along each segment. 

   Usage:    BreakOutline[outline, landmarks, {breakpoints}] 
			 where outline is the x,y coordinates of the outline curve, partitioned into blocks
			 of two, landmarks are the x,y coordinates of the points where the curve should be
             broken (also partitioned into blocks of two) and breakpoints is a list of integers
			 indicating the number of equally spaced points to be placed on each segment.
			 Note that there will be one more segment than there are landmarks.

	Created 15 January 2014.

	Updated 1 September 2014 to deal with cases where Nearest[] finds more than one point.
    Updated 27 April 2015 so that the points between the breaks are also spaced equally.

*)


BreakOutline[thiscurve_,theselandmarks_,ptnums_:Null]:=Block[{breakpoints,brokencurves,x},
If[ptnums==Null,Return["You must specify the number of semilandmarks for each segement of the broken outline."]];
If[Length[ptnums]!=Length[theselandmarks]+1,Return["You have not specified the number of semilandmarks for each segment of the broken outline."]];
If[Length[Dimensions[thiscurve]]!=2,Return["The outline coordinates must be partitioned into groups of two for 2D data or groups of three for 3D data."]];
If[Length[Dimensions[theselandmarks]]!=2,Return["The outline coordinates must be partitioned into groups of two for 2D data or groups of three for 3D data."]];

breakpoints=Sort[Prepend[Append[Flatten[{Nearest[thiscurve->Automatic,#][[1]]}&/@theselandmarks],Length[thiscurve]],0]];
brokencurves=thiscurve[[#]]&/@Table[breakpoints[[x-1]]+1;;breakpoints[[x]],{x,2,Length[breakpoints]}];
Do[
brokencurves[[1]]=Append[brokencurves[[1]],brokencurves[[2,1]]];
brokencurves=RotateRight[brokencurves],
{Length[brokencurves]}];

brokencurves=Table[EqualSpace[brokencurves[[x]],ptnums[[x]]+1],{x,Length[brokencurves]}];
brokencurves=Drop[#,-2]&/@brokencurves;
Return[brokencurves];
]





(* This function returns principal coordinates scores using either Euclidean or Gower distances. 

   Usage:    PrincipalCoordinates[data, (Gower)] 
			 where data is a matrix with specimens in the rows and variables in the columns.
			 and Gower is an optional parameter, which if set to True gives scores using Gower
			 distances.

	Created 21 March 2015.


*)

PrincipalCoordinates[data_,Gower_:False]:=Block[{dists,Q,uq,vq,wq,PCOscores,newdata},
If[Gower==False,newdata=data,newdata=Transpose[Table[(data[[1;;,i]]-Min[data[[1;;,i]]])/(Max[data[[1;;,i]]]-Min[data[[1;;,i]]]),{i,Length[data[[1]]]}]]];
dists=Table[Table[-0.5*Plus@@((newdata[[i]]-newdata[[j]])^2),{i,Length[newdata]}],{j,Length[newdata]}];
Q=Table[Table[dists[[i,j]]-Mean[dists[[1;;,j]]]-Mean[dists[[i,1;;]]]+Mean[Flatten[dists]],{i,Length[dists]}],{j,Length[dists]}];
{uq,vq,wq}=SingularValueDecomposition[Q];
PCOscores=Transpose[Transpose[uq]*Sqrt[Tr[vq,List]]];
Return[PCOscores];

]


(* This function returns the RV coefficient of Klingenberg (2009) and Escoufier (1973) for two blocks of data. 

   Usage:    RVCoefficient[data1, data2] 
			 where data1 and data2 are two blocks of columns from a matrix with specimens in the rows and variables in the columns.

	Created 23 March 2019.

*)

RVCoefficient[data1_,data2_]:=Block[{cm1,cm2,cm12,cm21},
cm1=Covariance[data1];
cm2=Covariance[data2];
cm12=Covariance[data1,data2];
cm21=Covariance[data2,data1];
Return[Tr[cm12 . cm21]/(Sqrt[Tr[cm1 . cm1]*Tr[cm2 . cm2]])];
];



(* This function returns the CR coefficient of Adams (2016) for two blocks of data. 

   Usage:    CRCoefficient[data1, data2] 
			 where data1 and data2 are two blocks of columns from a matrix with specimens in the rows and variables in the columns.

	Created 23 March 2019.

*)


CRCoefficient[data1_,data2_]:=Block[{cm1,cm2,cm12,cm21},
cm1=Covariance[data1];
cm2=Covariance[data2];
Do[cm1[[i,i]]=0,{i,Length[cm1]}];
Do[cm2[[i,i]]=0,{i,Length[cm2]}];
cm12=Covariance[data1,data2];
cm21=Covariance[data2,data1];
Return[Sqrt[Tr[cm12 . cm21]/(Sqrt[Tr[cm1 . cm1]*Tr[cm2 . cm2]])]];
];




(* This function performs a permutation test using the CR coefficient as described by Adams (2016) to test whether two blocks of 
columns have significantly more covariance within them than between them (i.e., whether they are "modular") and whether the same two 
blocks have more covariance between than within (i.e., whether they are integrated).  

   Usage:    AdamsCRTest[proc, landmarks1, landmarks2 (, dims, permuteby)] 
			 where proc is a block of procrustes superimposed landmarks (or any other data set) and landmarks1 and landmarks 2
			 are lists of landmarks by number, each of which is assumed to have two coordinates, x & y, unless specified in an 
			 option.  Options:  dims specifies the number of dimensions in each landmark (by default this is set to 2) and permuteby 
			 is a string "Landmarks" if all dimensions of each landmark are to be permuted togethr or "Variables" if each column is
			 to be permuted independently.  

	Created 23 March 2019.

*)


AdamsCRTest[proc_,lands1_,lands2_,dims_:2, permuteby_:"Landmarks"]:=Block[{RealCR,modulecolumns,randomcolumns,randCRs,pModular,pIntegrated,RandLands,permutemessage},
If[Length[proc[[1]]]!=(Length[lands1]+Length[lands2])*dims,Return["Number of landmarks does not match data.  Data have "<>ToString[Length[proc[[1]]]/dims]<>" "<>ToString[dims]<>"D landmarks and module specification has total of "<>ToString[Length[lands1]+Length[lands2]]<>".  Dims is set to "<>ToString[dims]<>"; Is this correct?"]];
modulecolumns=Flatten[Partition[Table[x,{x,Length[proc[[1]]]}],dims][[#]]]&/@{lands1,lands2} ;
RealCR=CRCoefficient[proc[[1;;,modulecolumns[[1]]]],proc[[1;;,modulecolumns[[2]]]]];
randCRs=Table[
If[permuteby=="Landmarks",
RandLands=Flatten[RandomSample[Partition[Table[x,{x,Length[proc[[1]]]}],dims]]],
RandLands=RandomSample[Table[x,{x,Length[proc[[1]]]}]]];
randomcolumns={RandLands[[1;;Length[lands1]*dims]],RandLands[[Length[lands1]*dims+1;;]]};
CRCoefficient[proc[[1;;,randomcolumns[[1]]]],proc[[1;;,randomcolumns[[2]]]]]
,{1000}];
pModular=Length[Select[randCRs,#>=RealCR&]]/1000;
pIntegrated=Length[Select[randCRs,#<=RealCR&]]/1000;
If[permuteby=="Landmarks",permutemessage="Permuted by landmarks.", permutemessage="Permuted by columns."];
Return[{"CR"->RealCR,"P(modular)"->pModular//N,"P(integrated)"->pIntegrated//N,Show[Histogram[randCRs,PlotRange->{{0,1.5},All}],Graphics[Line[{{RealCR,0},{RealCR,50}}]]],permutemessage}];
];

