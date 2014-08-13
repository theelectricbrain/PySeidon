function [elements]=get_elements(uvnode,region)

elements=find(uvnode(:,1)>=region(1) & uvnode(:,1)<=region(2) & uvnode(:,2)>=region(3) & uvnode(:,2)<=region(4));