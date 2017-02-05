#####################################################################
#
# MachedPairs(G)
#
# Calcula los pares cruzados de un grupo salvo conjugación y orden.
#
#
LoadPackage("Sonata");
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

MachedPairsNoCyclic:=function(G)
    local i, j, L,Pares,C,x;
    Pares:=[];
    C:=List(ConjugacyClassesSubgroups(G),Representative);
    L:=Filtered( C , x -> IsCyclic(x));
        for i in [1..Length(L)-1] do 
              for j in [i+1..Length(L)] do
                  if Size(L[i])>1  and Size(Intersection(L[i],L[j]))=1 and Size(DoubleCosets(G,L[i],L[j]))=1 
                      then Add(Pares, [L[i],L[j]]);
                  fi;
             od;
        od;
    return(Pares);
end;
MachedPairsNoCyclicNoNormal:=function(G)
    local i, j, L,Pares,C,x,R;
    Pares:=[];
    C:=List(ConjugacyClassesSubgroups(G),Representative);
    R:=Filtered( C , x -> IsCyclic(x));
    L:=Filtered(R,x -> IsNormal(G,x)=false);
        for i in [1..Length(L)-1] do 
              for j in [i+1..Length(L)] do
                  if Size(L[i])>1  and Size(Intersection(L[i],L[j]))=1 and Size(DoubleCosets(G,L[i],L[j]))=1 
                      then Add(Pares, [L[i],L[j]]);
                  fi;
             od;
        od;
    return(Pares);
end;

TestOneGroup := G -> Size(MachedPairsNoCyclicNoNormal(G))>1 ;
    
TestOneOrderEasy := function(n)
local i;
for i in [1..NrSmallGroups(n)] do
  if TestOneGroup( SmallGroup( n, i ) ) then
    return [n,i];
  fi;
od;
return fail;
end;

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