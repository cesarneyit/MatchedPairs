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
    L:=Subgroups( G );
        for i in [1..Length(L)-1] do 
              for j in [i+1..Length(L)] do
                  if  Size(Intersection(L[i],L[j]))=1 and Size(DoubleCosets(G,L[i],L[j]))=1 
                      then Add(Pares, [L[i],L[j]]);
                  fi;
             od;
        od;
    return(Pares);
end;