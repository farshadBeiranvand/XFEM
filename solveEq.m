function U = solveEq(Q,NXE,NYE,KIN,gdofIN,force)

F=force;

% Us=pinv(KIN)*F;
Us=lsqr(KIN,F,1e-7,1e6);
U=Us;
end