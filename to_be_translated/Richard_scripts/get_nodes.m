function [nodes]=get_nodes(nodexy,region)

nodes=find(nodexy(:,1)>=region(1) & nodexy(:,1)<=region(2) & nodexy(:,2)>=region(3) & nodexy(:,2)<=region(4));