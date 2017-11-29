include("LP_matrix.jl")

function create_SDLP(A,b,c)
  dimension = size(A);
  dv_num    = dimension[2];
  sv_num    = dimension[1];

  #add slack variable to the matrix
  Aug       = [A eye(sv_num)];
  cug       = [c; zeros(sv_num)];

  #initialize basic variable and nonbasic variable list
  nbv       = collect(1:1:dv_num);
  dv        = collect(1:1:dv_num);
  bv        = collect(dv_num+1:1:dv_num+sv_num);
  sv        = collect(dv_num+1:1:dv_num+sv_num);

  return LP_matrix(Aug,b,cug,0,bv,nbv,dv,sv)
end

function find_enter(LP)
  Ab = LP.A[:,LP.bv];
  An = LP.A[:,LP.nbv];
  cb = LP.c[LP.bv,:];
  cn = LP.c[LP.nbv,:];

  ceof    = cn' - cb'*inv(Ab)*An;
  enter_v = findfirst(ceof, maximum(ceof));
  return enter_v;
end

function get_constant(LP)
  Ab = LP.A[:,LP.bv];
  return inv(Ab)*LP.b;
end

function get_dict_coef(LP,index)
  Aj      = LP.A[:,index];
  coef    = inv(LP.A[:,LP.bv])*Aj
  return coef
end

function find_exit(LP,index)
  b     = LP.b;
  aj    = get_dict_coef(LP,index);
  bound = b./aj
  t_bound = minimum(bound);
  #assume all ratio >=0 for now i.e not degeneracy
  return [LP.bv[findfirst(bound, t_bound)],t_bound];
end

function next_dic(LP)
  #get the information for the next step
  #1. find the enter variable by computing
  #     CnT - CbT*Ab-1*An
  #     the largest coeffcients
  #2. find the exiting variable by computing
  #     b - t*Ab-1*Aj
  #3. Update b for the next step
  #     b* = b - t*d and replace the enter entries with t

  #get the enter and exit variable
  enter_v = find_enter(LP);
  exit_v  = find_exit(LP,enter_v);

  #compute new b
  d       = get_dict_coef(LP,enter_v);
  b       = LP.b - exit_v[2].*d

  #update constant for the entering entries
  exit_index    = findfirst(LP.bv,exit_v[1]);
  b[exit_index] = exit_v[2];
  LP.b          = b;

  #update basic varialb
  LP.bv[exit_index] = enter_v;

  #update non basic variable
  enter_index         = findfirst(LP.nbv,enter_v);
  LP.nbv[enter_index] = exit_v[1];
end

#we can just check the entering column at each step
#because eventually, we will get to the unbounded column
function detect_unboundness(enter_coe, dict_coe)

end

#detect cycle by keeping track with the objective value
#if the objective value hasn't change for 5 dictionary
#assume there is a cycle and swtich to bland's rule
function detect_cycle()

end

#this function will detect if the LP
#enter optimal form
function detect_optimality()
end

A = [1.0 4.0 0.0; 3 -1 1];
b_e = [1 2];
c = [4 1 2];
SDl = create_SDLP(A,b_e',c');
SDL_next = next_dic(SDl);
