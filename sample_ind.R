cppFunction(
'LogicalVector new_sample(NumericMatrix xy, int N, double buffer, LogicalVector use, IntegerVector new_order)
{
  int maxid = use.size();
  
	use[new_order[0]] = 1; // Automatically include the first point in the sample = 1st wolverine

	int ct = 1; // Number of individuals included so far
	int id = 1; // Current individual to check
	double mindist;
	double newdist;
  bool * trydist = new bool[maxid];
  
  for(int i = 0; i < maxid; i++)
  {
    newdist = sqrt(pow(xy(new_order[id],1) - xy(new_order[i],1), 2)+ pow(xy(new_order[id],0) - xy(new_order[i],2), 0))/1000.0;
    if(newdist <= buffer){ 
      trydist[new_order[i]] = 0;
    } else { 
      trydist[new_order[i]] = 1; 
    }
  } 

	while((ct < N) && (id < maxid))		// NOT <=*maxid because id starts at 0
	{ 
		if(trydist[new_order[id]] == 1){
      use[new_order[id]] = 1;
      ct++;
     for(int i=0; i < maxid; i++)
     {
      newdist = sqrt(pow(xy(new_order[id],1) - xy(new_order[i],1), 2)+ pow(xy(new_order[id],0) - xy(new_order[i],2), 0))/1000.0;
      if(newdist <= buffer) 
        trydist[new_order[i]] = 0;    
     }
    } 
		id++;
	}
  delete [] trydist;
  return use;
 }'
 )