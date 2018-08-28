function [biasmap, hit_bias, NoHit, NoFalse, NoMiss, SumMiss, SumFalse, sumSimhit] = ValidationFunction(simulated, observed,trr)
vector1 = simulated((simulated >=trr | observed >=trr) & (observed>=0 & simulated>=0));
vector2 = observed((simulated >=trr | observed >=trr) & (observed>=0 & simulated>=0));
simulated = vector1; observed = vector2; 
HitError = simulated(simulated > 0 & observed > 0) - observed(simulated > 0 & observed > 0);
hit_bias = sum(HitError); 
biasmap = sum(simulated)/sum(observed);
simhit = simulated(simulated > 0 & observed > 0);
NoHit = sum(simulated > 0 & observed > 0);
sumSimhit = sum(simhit);
obsmiss = observed(simulated == 0 & observed > 0);
NoMiss = sum(simulated == 0 & observed > 0);
SumMiss = sum(obsmiss);
simfalse = simulated(simulated > 0 & observed == 0);
NoFalse = sum(simulated > 0 & observed == 0);
SumFalse = sum(simfalse);




