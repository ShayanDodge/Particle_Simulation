function [fi] = create_linear_index_list(array_3d,i_list,j_list,k_list)
i_size=size(i_list,2);
j_size=size(j_list,2);
k_size=size(k_list,2);
number_of_indices=i_size*j_size*k_size;
I=zeros(number_of_indices,1);
J=zeros(number_of_indices,1);
K=zeros(number_of_indices,1);
ind=1;
for mk=k_list(1):k_list(k_size)
    for mj=j_list(1):j_list(j_size)
        for mi=i_list(1):i_list(i_size)
            I(ind)=mi;
            J(ind)=mj;
            K(ind)=mk;
            ind=ind+1;
        end
    end
end
fi=sub2ind(size(array_3d),I,J,K);
