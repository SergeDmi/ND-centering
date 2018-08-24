function [M,norms]=normalize_rows_ND(M,N)
% normalize_rows_ND(M,N)
% Normalizes a column of N-dimension row vectors
% N is given as an argument to save time... super duper useful
% Still making it two lines for clarity
norms=sqrt(sum(M.^2,2));
M(:,:)=M./(norms*ones(1,N));
end
