function F = pwmSwitchneg(x,modphase,phase,mf,Tac) %x(1) is either 1 or -1, x(2) is the phase of the 

F = abs(x*4*mf/Tac - 1 + modphase*cos(x + phase)); 


