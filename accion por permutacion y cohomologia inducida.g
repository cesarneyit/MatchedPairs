MachedPairsNoNormal:=function(G)
    local i, j, L,Pares,C,x;
    Pares:=[];
    C:=List(ConjugacyClassesSubgroups(G),Representative);
    L:=Filtered( C , x -> IsNormal(G,x)=false);
        for i in [1..Length(L)-1] do 
              for j in [i+1..Length(L)] do
                  if Size(L[i])>1  and Size(Intersection(L[i],L[j]))=1 and Size(DoubleCosets(G,L[i],L[j]))=1 
                      then Add(Pares, [L[i],L[j]]);
                  fi;
             od;
        od;
    return(Pares);
end;




LoadPackage("hap");
G:=AlternatingGroup(5);;
R:=MachedPairsNoNormal(G)[1][1]; # solo lo estoy usando para buscar un subgrupo interesante
hom:=ActionHomomorphism(G, RightCosets(G, R), OnRight); # crea un grupo de permutaciones y el morfismo
P:=Image(hom); 
LoadPackage("hap");
M:=PermToMatrixGroup(P,Size(G)/Size(R)); # accion por permutacion natural
rep:=hom * M; # representacion por permutation # accion por permutacion del grupo original
Res:=ResolutionGenericGroup(G,4); 
Cohomology(HomToIntegralModule(Res,rep),3);

#
# Ahora vamos a crear el Delta
#
B:=IdentityMat(12);
BaseW:=List([1..12-1],x-> B[x]-B[12]);
V:= VectorSpace( Rationals, BaseW);
BaseV:=Basis( V, BaseW );

#provemos
q:=Elements(G)[5];# un elemento al azar en G
Per:=Image(rep,q); # la matriz asociada
Coefficients( BaseV,Per* BaseW[11] );
MatrixPer:=y-> List([1..12-1],x->Coefficients( BaseV,y* BaseW[x] ) ); # matriz asociada a Per


# vamos a crear el mapa a GL(n,Z)

gensG:=GeneratorsOfGroup(G); 
gensA:=[]; # lista vacia donde estaran las imagenes de los generadores de G




gensG:=GeneratorsOfGroup(G);
gensA:=[];

for g in gensG do
Append(gensA,[MatrixPer(Image(rep,g))]);
od;

A:=Group(gensA);
GhomA:=GroupHomomorphismByImagesNC(G,A,gensG,gensA);

return GhomA;


end);

# matriz en bloque 
M1:=gensA[1];
M2:=gensA[2];
BlockMatrix([  [1,1,M1], [2,2,M2] ],2,2 );