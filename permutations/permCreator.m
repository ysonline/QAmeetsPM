n = 9; %all permutations of ints between 0 and n (inclusive)
A = perms(0:n); %this call takes just 1 second for n = 9 :)

fid = fopen([int2str(n+1) '!.txt'],'wt'); %wt for writing in text mode
for row = 0:factorial(n+1)-1
    fprintf(fid,'%d\t', A(row+1, :)); %+1: matlab indices start from 1
    fprintf(fid, '\n');
end
fclose(fid);
%whole execution takes 110 seconds for n = 9 :(

%%%%%%% prolog code for fantasy purposes :) %%%%%%%
%insert(X,Y,[X|Y]).
%insert(X,[A|Y],[A|R]) :- insert(X,Y,R).
%permut(_,[],0).
%permut(X,[A|RestPerm],N) :- N>0,M is N-1,
%                insert(A,Rest,X),permut(Rest,RestPerm,M).
%?- permut([0,1,2,3,4,5,6,7,8,9],R,10).
