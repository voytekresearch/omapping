% For the Desikan-Killiany atlas
% Export scouts from Brainstorm using label DK_Scout

for i = 1 : 68
    if i<10
        varstr=strcat('DK_Scout0',int2str(i));
    else
        varstr=strcat('DK_Scout',int2str(i));
    end
    
    s=eval(varstr);
    Scout_Vertex.Label{i,1}=s.Label;
    Scout_Vertex.Vertices{i,1}=s.Vertices;
end
