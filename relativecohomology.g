LoadPackage("Sonata");
LoadPackage("hap");
#####################################################################
#
# MachedPairs(G)
#
# Calcula los pares cruzados de un grupo salvo conjugación y orden.
#
MachedPairs:=function(G)
    local i, j, L, Pares;
    Pares:=[];
    L:=List(ConjugacyClassesSubgroups(G),Representative);
        for i in [1..Length(L)-1] do 
              for j in [i+1..Length(L)] do
                  if Size(L[i])>1 and Size(Intersection(L[i],L[j]))=1 and Size(DoubleCosets(G,L[i],L[j]))=1 
                      then Add(Pares, [L[i],L[j]]);
                  fi;
             od;
        od;
    return(Pares);
end;

#####################################################################
#
# RelCoh(G,R,n)
#
# Calcula H^n(G,F;Z).
#
#
LoadPackage("Sonata");
LoadPackage("hap");


#####################################################################
#
# PermMod(G,R,S)
#
# Calcula el modulo de permutacion G-->GL(V) 
#  del G-conjunto G/RUG/S
#
PermMod:= function(G,Par)
local hom_R,hom_S, M_R,M_S,rep_S,rep_R,gensG,A,gensA,GhomA,R,S;
R:=Par[1];
S:=Par[2];
gensG:=GeneratorsOfGroup(G);

hom_R:=ActionHomomorphism(G, RightCosets(G, R), OnRight); # crea un grupo de permutaciones y el morfismo
hom_S:=ActionHomomorphism(G, RightCosets(G, S), OnRight); # crea un grupo de permutaciones y el morfismo
M_R:=PermToMatrixGroup(Image(hom_R),Size(G)/Size(R)); # accion por permutacion natural
M_S:=PermToMatrixGroup(Image(hom_S),Size(G)/Size(S)); # accion por permutacion natural
rep_R:=hom_R * M_R; # representacion por permutation # accion por permutacion del grupo original
rep_S:=hom_S * M_S; # representacion por permutation # accion por permutacion del grupo original

gensA:=List(gensG,x-> BlockMatrix([  [1,1,Image( rep_R,x)], [2,2, Image( rep_R,x)] ],2,2 ));

A:=Group(gensA);
#GhomA:=GroupHomomorphismByImagesNC(G,A,gensG,gensA);
GhomA:=GroupHomomorphismByImages(G,A,gensG,gensA);

return GhomA;


end;


#####################################################################
#
# Delta(P)
#
# Calcula la matriz  
#  P lista de matrices
#

aumentacion:=function(P)
local n,B,BaseV,M,MatrixPer,g,gensA,V ;
gensA:=[];# lista donde estaran las matrices
n:=Size(P[1]);
B:=IdentityMat(n);
BaseV:=List([1..n-1],x-> B[x]-B[n]);
V:= VectorSpace( Rationals, BaseV);
BaseV:=Basis( V, BaseV );
MatrixPer:=y-> List([1..n-1],x->Coefficients( BaseV,y* BaseV[x] ) ); # funcion que calcula la matriz
for g in P do
Append(gensA,TransposedMat(MatrixPer(g)));
od;

return(gensA);
end;


#####################################################################
#
# MatrixPermRe(G,R)
#
# Calcula la representacion sobre Delta   
#  
#
MatrixPermRe:=function(G,Par)

local gensG, gensA, rep, P,g,A,GhomA;

gensG:=GeneratorsOfGroup(G);
rep:=PermMod(G,Par);
#P:=List(gensG,x-> TransposedMat(Image(rep,x^-1));
#P:=List(gensG,x-> Image(rep,x^-1);
P:=List(gensG,x-> Image(rep,x));
gensA:=aumentacion(P);

A:=Group(gensA);
#GhomA:=GroupHomomorphismByImagesNC(G,A,gensG,gensA);
GhomA:=GroupHomomorphismByImages(G,A,gensG,gensA);

return GhomA;




end;

#####################################################################
#
# RelCoh(G,Par,n)
#
# Calcula H^n(G,ker(e) donde...   
#  
#
RelCoh:=function(G,Par,n)
local M,R;
M:=MatrixPermRe(G,Par);
R:=ResolutionGenericGroup(G,n+1);
return Cohomology(HomToIntegralModule(R,M),n);
end;

#####################################################################
#
# SumCoho(G,Par,n)
#
# Calcula H_n(S,Z)+H_n(R,Z) donde...   
#  
#
SumCoho:=function(G,Par,n)
local M,R;
M:=PermMod(G,Par);
R:=ResolutionGenericGroup(G,n+1);
return Cohomology(HomToIntegralModule(R,M),n);
end;
