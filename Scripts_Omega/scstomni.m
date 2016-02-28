function [] = scstomni(Vertex_co_scs)
%% Function scstomni converts vertex co-ordinates from SCS to MNI
%% Input 
%%      Vertex_co_scs - Co-ordinates of all vertices in SCS
%% 
%% Before calling scstomni export MRI data from brainstorm using label sMri
%% Also Compute MNI Transformation 

% Convert vertex co-ordinates from SCS to MNI
for i = 1:length(Vertex_co_scs)
    Vertex_co_mni(i,:) = cs_convert(sMri, 'scs', 'mni', Vertex_co_scs(i,:));
end   

end