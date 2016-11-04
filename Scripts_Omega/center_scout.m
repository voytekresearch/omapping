function [scout_cent] = center_scout(Label,Vertices,Vertex_co) 
%% Function scout_cent finds the center of each scout in an atlas.
%% Input 
%%      Label - Scout name
%%      Vertices - Set of vertices that make up each scout
%%      Vertex_co - Co-ordinates of all vertices
%% Output
%%      scout_cent - Array containing the vertex that is at the center of each scout

no_of_scouts = length(Label);
scout_cent= zeros(no_of_scouts,1);
x=1;y=2;z=3;

for scout_no = 1 : no_of_scouts
    fprintf('Finding center of scout no %d -',scout_no);
    disp(Label(scout_no));
    vertices = Vertices{scout_no,1};
    xsum = 0; ysum=0;zsum=0;
    no_of_vertices = length(vertices);
    
    %% Sum x, y, z co-ordinates
    for j = 1:no_of_vertices
        vno = vertices(j);
        xsum = xsum+Vertex_co(vno,x);
        ysum = ysum+Vertex_co(vno,y);
        zsum = zsum+Vertex_co(vno,z);
    end
    
    %% Find theoretical center
    xavg = xsum/no_of_vertices;
    yavg = ysum/no_of_vertices;
    zavg = zsum/no_of_vertices;
    
    %% Find vertex closest to theoretical center
    vno = vertices(1);
    mindist = sqrt((xavg-Vertex_co(vno,x))^2+(yavg-Vertex_co(vno,y))^2+(zavg-Vertex_co(vno,z))^2);
    scout_cent(scout_no,1) = vno;
    for j = 2:no_of_vertices
        vno = vertices(j);
        dist = sqrt((xavg-Vertex_co(vno,x))^2+(yavg-Vertex_co(vno,y))^2+(zavg-Vertex_co(vno,z))^2);
        if (dist<mindist)
            mindist = dist;
            scout_cent(scout_no,1) = vno;
        end
    end
    
end
end      