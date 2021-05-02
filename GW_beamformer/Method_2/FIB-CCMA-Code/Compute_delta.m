function [ delta_matrix] = Compute_delta(r_vect, M_vect)

% % Compute the diffuse noise pseudo-coherence matrix with CCMA
% clc;clear all;close all;
% r_vect = [0.01 0.02]; M_vect = [4 8];

for index = 1 : length(r_vect);
    
phi_vector_2 = (0 : M_vect(index)-1)'*2*pi/M_vect(index);
coordinate_index = r_vect(index) * [cos(phi_vector_2) sin(phi_vector_2)];

if (index ==1)
    coordinate = coordinate_index;
else
    coordinate = [coordinate; coordinate_index];
end
end
delta_matrix = dist(coordinate');

% M_sum = length(coordinate(:,1));
% delta_matrix_2 = zeros(M_sum);
% for index1 = 1 :    M_sum
%     coordinate_index1 = coordinate(index1,:);
%     for index2 = 1 : M_sum
%         coordinate_index2 = coordinate(index2,:);
%         distance = dist(coordinate_index1,coordinate_index2');
%         delta_matrix_2(index1,index2) = distance;
%     end
% end
   
end

